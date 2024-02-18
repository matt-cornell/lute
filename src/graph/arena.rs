//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use crate::graph::bitfilter::BitFiltered;
use crate::graph::compact::GraphCompactor;
use crate::graph::isomorphism::*;
use crate::molecule::{Atom, Bond};
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::data::DataMap;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallbitvec::SmallBitVec;
use smallvec::{smallvec, SmallVec};
use std::cell::{Ref, RefCell, RefMut};
use std::collections::HashMap;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

type Ix = petgraph::graph::DefaultIx;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
struct InterFragBond {
    an: Ix,
    ai: Ix,
    bn: Ix,
    bi: Ix,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BrokenMol {
    frags: SmallVec<[Ix; 8]>,
    bonds: SmallVec<[InterFragBond; 8]>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MolRepr {
    Atomic(SmallBitVec),
    Broken(BrokenMol),
    Redirect(Ix),
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena {
    graph: StableUnGraph<Atom, Bond>,
    parts: SmallVec<[(MolRepr, Ix); 16]>,
}
impl Arena {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn graph(&self) -> &StableUnGraph<Atom, Bond> {
        &self.graph
    }

    fn contains_group_impl(&self, mol: usize, group: usize, seen: &mut SmallBitVec) -> bool {
        if mol == group {
            return true;
        }
        if seen[mol] {
            return false;
        }
        seen.set(mol, true);
        match self.parts.get(mol) {
            Some((MolRepr::Broken(b), _)) => b
                .frags
                .iter()
                .any(|f| self.contains_group_impl(*f as _, group, seen)),
            Some((MolRepr::Redirect(r), _)) => self.contains_group_impl(*r as _, group, seen),
            _ => false,
        }
    }

    /// Check if `mol` contains `group`
    pub fn contains_group(&self, mol: usize, mut group: usize) -> bool {
        while let Some((MolRepr::Redirect(r), _)) = self.parts.get(group) {
            group = *r as _;
        }
        let mut seen = SmallBitVec::from_elem(self.parts.len(), false);
        self.contains_group_impl(mol, group, &mut seen)
    }

    pub fn molecule(&self, mol: usize) -> Molecule<RefAcc> {
        Molecule::from_arena(self, mol)
    }

    pub fn insert_mol<G>(&mut self, mol: G) -> usize
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + DataMap
            + GraphProp<EdgeType = Undirected>
            + GraphRef
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
    {
        let compacted = (0..self.parts.len())
            .filter_map(|i| {
                if let (MolRepr::Atomic(a), _) = &self.parts[i] {
                    Some((
                        i,
                        GraphCompactor::<BitFiltered<&StableUnGraph<Atom, Bond>>>::new(
                            BitFiltered::new(unsafe {&*std::ptr::addr_of!(self.graph)}, a.clone()),
                        ),
                    ))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        // keep track of matched atoms so there's no overlap
        let mut matched = SmallBitVec::from_elem(mol.node_count(), false);
        // keep track of found isomorphisms, don't try to handle them in the search
        let mut found = SmallVec::<[_; 8]>::new();
        for (n, cmp) in &compacted {
            for ism in
                subgraph_isomorphisms_iter(&cmp, &mol, &mut Atom::eq_or_r, &mut PartialEq::eq)
            {
                if ism.iter().any(|&i| matched[i]) {
                    continue;
                }
                ism.iter().for_each(|&i| matched.set(i, true));
                found.push((*n, ism));
            }
        }
        let (ret, news): (_, SmallVec<[_; 8]>) = if found.len() == 1 {
            // simple case: no subgraph isomorhpisms found
            let mut map = vec![NodeIndex::end(); mol.node_count()];
            let mut bits = SmallBitVec::new();
            mol.node_references().for_each(|n| {
                let idx = self.graph.add_node(*n.weight());
                map[mol.to_index(n.id())] = idx;
                let idx = idx.index();
                if idx >= bits.len() {
                    bits.resize(idx, false);
                }
                bits.set(idx, true);
            });
            mol.edge_references().for_each(|e| {
                self.graph.add_edge(
                    map[mol.to_index(e.source())],
                    map[mol.to_index(e.target())],
                    *e.weight(),
                );
            });
            let out = self.parts.len() as Ix;
            self.parts
                .push((MolRepr::Atomic(bits), mol.node_count() as Ix));
            (out, smallvec![out])
        } else {
            let mut frags = SmallVec::with_capacity(found.len() + 1);
            frags.push(self.parts.len() as Ix);
            frags.extend(found.iter().map(|(i, _)| *i as Ix));
            let mut bonds = SmallVec::new();

            for (n, (i, ism)) in found.iter().enumerate() {
                let cmp = &compacted[*i as usize].1;
                for j in 0..ism.len() {
                    let idx = cmp.node_map[j];
                    let atom = self.graph[idx];
                    if atom.protons == 0 {
                        let an = 0;
                        let ai = ism[j] as Ix;
                        let bn = n as Ix + 1;
                        let bi = idx.index() as Ix;
                        bonds.push(InterFragBond { an, ai, bn, bi });
                    } else {
                       todo!()
                    }
                }
            }

            self.parts.push((
                MolRepr::Broken(BrokenMol { frags, bonds }),
                mol.node_count() as Ix,
            ));

            todo!()
        };
        for new in news {
            let Some((MolRepr::Atomic(a), _)) = self.parts.get(new as usize) else {
                continue;
            };
            for (n, cmp) in &compacted {
                if is_isomorphic_matching(
                    &cmp,
                    &GraphCompactor::<BitFiltered<&StableUnGraph<Atom, Bond>>>::new(
                        BitFiltered::new(&self.graph, a.clone()),
                    ),
                    &mut Atom::eq_match_r,
                    &mut PartialEq::eq,
                ) {
                    self.parts[new as usize].0 = MolRepr::Redirect(*n as _);
                    break;
                }
            }
        }
        ret as _
    }
}

/// A `Container` corresponds to a reaction vessel. It's at this layer that actual reactions are
/// handled.
#[derive(Debug, Clone)]
pub struct Container {
    pub quantities: HashMap<usize, f64>,
    pub temperature: f64,
}

/// A `Molecule` acts like a graph, and can have graph algorithms used on it. It's immutable, with all
/// mutations making (efficient) copies.
#[derive(Debug, Clone, Copy)]
pub struct Molecule<R> {
    arena: R,
    index: usize,
}
impl<R> Molecule<R> {
    pub fn from_arena<'a, 'b: 'a, A: ArenaAccessible<Access<'a> = R> + 'a>(
        arena: &'b A,
        index: usize,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor(),
            index,
        }
    }

    pub fn from_mut_arena<'a, 'b: 'a, A: ArenaAccessibleMut<AccessMut<'a> = R> + 'a>(
        arena: &'b A,
        index: usize,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor_mut(),
            index,
        }
    }

    pub fn arena(&self) -> R::Ref<'_>
    where
        R: ArenaAccessor,
    {
        self.arena.get_arena()
    }

    pub fn arena_mut(&self) -> R::RefMut<'_>
    where
        R: ArenaAccessorMut,
    {
        self.arena.get_arena_mut()
    }

    pub fn contains(&self, group: usize) -> bool
    where
        R: ArenaAccessor,
    {
        self.arena.get_arena().contains_group(self.index, group)
    }
}

/// This trait handles the access to the backing arena. Rather than just passing around references,
/// this allows for lock guards to be used while not forcing them to live for as long as the
/// accessor.
pub trait ArenaAccessor {
    type Ref<'a>: Deref<Target = Arena>
    where
        Self: 'a;

    fn get_arena<'a>(&'a self) -> Self::Ref<'a>
    where
        Self: 'a;
}

/// This trait provides mutable access to the underlying arena.
pub trait ArenaAccessorMut: ArenaAccessor {
    type RefMut<'a>: DerefMut<Target = Arena>
    where
        Self: 'a;

    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        Self: 'a;
}

/// This trait allows access to a backing arena.
pub trait ArenaAccessible {
    type Access<'a>: ArenaAccessor + 'a
    where
        Self: 'a;

    fn get_accessor(&self) -> Self::Access<'_>;
}

/// This trait also allows access to a graph through an internally mutable type.
pub trait ArenaAccessibleMut: ArenaAccessible {
    type AccessMut<'a>: ArenaAccessorMut + 'a
    where
        Self: 'a;

    fn get_accessor_mut(&self) -> Self::AccessMut<'_>;
}

impl<T: ArenaAccessible> ArenaAccessible for &T {
    type Access<'a> = T::Access<'a> where Self: 'a;

    fn get_accessor(&self) -> Self::Access<'_> {
        T::get_accessor(self)
    }
}
impl<T: ArenaAccessibleMut> ArenaAccessibleMut for &T {
    type AccessMut<'a> = T::AccessMut<'a> where Self: 'a;

    fn get_accessor_mut(&self) -> Self::AccessMut<'_> {
        T::get_accessor_mut(self)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct PtrAcc(*mut Arena);
impl ArenaAccessor for PtrAcc {
    type Ref<'a> = &'a Arena;
    fn get_arena<'a>(&'a self) -> Self::Ref<'a> where Self: 'a {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe {
            &*self.0
        }
    }
}
impl ArenaAccessorMut for PtrAcc {
    type RefMut<'a> = &'a mut Arena;
    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a> where Self: 'a {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe {
            &mut *self.0
        }
    }
}

/// Wrapper type around a `*mut Arena`. It has an `unsafe` constructor because the
/// `ArenaAccessible` implementation can't be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArenaPtr(*mut Arena);
impl ArenaPtr {
    pub const unsafe fn new(ptr: *mut Arena) -> Self {
        Self(ptr)
    }
    pub fn get_ptr(self) -> *mut Arena {
        self.0
    }
}
impl ArenaAccessible for ArenaPtr {
    type Access<'a> = PtrAcc;

    fn get_accessor(&self) -> PtrAcc {
        PtrAcc(self.0)
    }
}
impl ArenaAccessibleMut for ArenaPtr {
    type AccessMut<'a> = PtrAcc;

    fn get_accessor_mut(&self) -> PtrAcc {
        PtrAcc(self.0)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RefAcc<'a>(&'a Arena);
impl<'a> ArenaAccessor for RefAcc<'a> {
    type Ref<'b> = &'b Arena where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0
    }
}

impl ArenaAccessible for Arena {
    type Access<'a> = RefAcc<'a>;

    fn get_accessor(&self) -> RefAcc {
        RefAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RefCellAcc<'a>(&'a RefCell<Arena>);
impl<'a> ArenaAccessor for RefCellAcc<'a> {
    type Ref<'b> = Ref<'b, Arena> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.borrow()
    }
}
impl<'a> ArenaAccessorMut for RefCellAcc<'a> {
    type RefMut<'b> = RefMut<'b, Arena> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.borrow_mut()
    }
}

impl ArenaAccessible for RefCell<Arena> {
    type Access<'a> = RefCellAcc<'a>;

    fn get_accessor(&self) -> RefCellAcc {
        RefCellAcc(self)
    }
}
impl ArenaAccessibleMut for RefCell<Arena> {
    type AccessMut<'a> = RefCellAcc<'a>;

    fn get_accessor_mut(&self) -> RefCellAcc {
        RefCellAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RwLockAcc<'a>(&'a RwLock<Arena>);
impl<'a> ArenaAccessor for RwLockAcc<'a> {
    type Ref<'b> = RwLockReadGuard<'b, Arena> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.read()
    }
}
impl<'a> ArenaAccessorMut for RwLockAcc<'a> {
    type RefMut<'b> = RwLockWriteGuard<'b, Arena> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.write()
    }
}

impl ArenaAccessible for RwLock<Arena> {
    type Access<'a> = RwLockAcc<'a>;

    fn get_accessor(&self) -> RwLockAcc {
        RwLockAcc(self)
    }
}
impl ArenaAccessibleMut for RwLock<Arena> {
    type AccessMut<'a> = RwLockAcc<'a>;

    fn get_accessor_mut(&self) -> RwLockAcc {
        RwLockAcc(self)
    }
}

impl<T: ArenaAccessible> ArenaAccessible for Rc<T> {
    type Access<'a> = T::Access<'a> where T: 'a;

    fn get_accessor(&self) -> Self::Access<'_> {
        T::get_accessor(self)
    }
}
impl<T: ArenaAccessibleMut> ArenaAccessibleMut for Rc<T> {
    type AccessMut<'a> = T::AccessMut<'a> where T: 'a;

    fn get_accessor_mut(&self) -> Self::AccessMut<'_> {
        T::get_accessor_mut(&**self)
    }
}

impl<T: ArenaAccessible> ArenaAccessible for Arc<T> {
    type Access<'a> = T::Access<'a> where T: 'a;

    fn get_accessor(&self) -> Self::Access<'_> {
        T::get_accessor(self)
    }
}
impl<T: ArenaAccessibleMut> ArenaAccessibleMut for Arc<T> {
    type AccessMut<'a> = T::AccessMut<'a> where T: 'a;

    fn get_accessor_mut(&self) -> Self::AccessMut<'_> {
        T::get_accessor_mut(&**self)
    }
}
