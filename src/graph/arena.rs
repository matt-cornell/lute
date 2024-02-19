//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use crate::graph::bitfilter::BitFiltered;
use crate::graph::compact::GraphCompactor;
use crate::graph::isomorphism::*;
use crate::molecule::{Atom, Bond};
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::data::DataMap;
use petgraph::graph::{DefaultIx, IndexType};
use petgraph::prelude::*;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::cell::{Ref, RefCell, RefMut};
use std::collections::HashMap;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

const ATOM_BIT_STORAGE: usize = 8;

type Graph<Ix> = StableGraph<Atom, Bond, Undirected, Ix>;
type BSType = crate::utils::bitset::BitSet<usize, ATOM_BIT_STORAGE>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
struct InterFragBond<Ix> {
    an: Ix,
    ai: Ix,
    bn: Ix,
    bi: Ix,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BrokenMol<Ix> {
    frags: SmallVec<Ix, 8>,
    bonds: SmallVec<InterFragBond<Ix>, 8>,
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Eq)]
enum MolRepr<Ix> {
    Atomic(BSType),
    Broken(BrokenMol<Ix>),
    Redirect(Ix),
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena<Ix: IndexType = DefaultIx> {
    graph: Graph<Ix>,
    parts: SmallVec<(MolRepr<Ix>, Ix), 16>,
}
impl<Ix: IndexType> Arena<Ix> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
    }

    fn contains_group_impl(&self, mol: Ix, group: Ix, seen: &mut BSType) -> bool {
        if mol == group {
            return true;
        }
        if seen.get(mol.index()) {
            return false;
        }
        seen.set(mol.index(), true);
        match self.parts.get(mol.index()) {
            Some((MolRepr::Broken(b), _)) => b
                .frags
                .iter()
                .any(|f| self.contains_group_impl(*f as _, group, seen)),
            Some((MolRepr::Redirect(r), _)) => self.contains_group_impl(*r, group, seen),
            _ => false,
        }
    }

    /// Check if `mol` contains `group`
    pub fn contains_group(&self, mol: Ix, mut group: Ix) -> bool {
        while let Some((MolRepr::Redirect(r), _)) = self.parts.get(group.index()) {
            group = *r;
        }
        let mut seen = BSType::with_capacity(self.parts.len());
        self.contains_group_impl(mol, group, &mut seen)
    }

    /// Get a graph of the molecule at the given index. Note that `Molecule::from_arena` could give
    /// better results as it can borrow from `RefCell`s and `RwLock`s.
    pub fn molecule(&self, mol: Ix) -> Molecule<Ix, RefAcc<Ix>> {
        Molecule::from_arena(self, mol)
    }

    /// Insert a molecule into the arena, deduplicating common parts. May misbehave if multiple
    /// disjoint molecules are in the graph.
    #[allow(clippy::needless_range_loop)]
    pub fn insert_mol<G>(&mut self, mol: G) -> Ix
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
                        GraphCompactor::<BitFiltered<&Graph<Ix>, _, ATOM_BIT_STORAGE>>::new(
                            BitFiltered::new(
                                unsafe { &*std::ptr::addr_of!(self.graph) },
                                a.clone(),
                            ),
                        ),
                    ))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        // keep track of matched atoms so there's no overlap
        let mut matched: SmallVec<_, 16> = smallvec![(0, 0); mol.node_bound()];
        // keep track of found isomorphisms, don't try to handle them in the search
        let mut found = SmallVec::<_, 8>::new();
        for (k, (n, cmp)) in compacted.iter().enumerate() {
            for ism in
                subgraph_isomorphisms_iter(&cmp, &mol, &mut Atom::eq_or_r, &mut PartialEq::eq)
            {
                if ism.len() == mol.node_count() {
                    // simplest case, this molecule already exists
                    return Ix::new(*n);
                }
                if ism.iter().any(|&i| matched[i].0.index() != 0) {
                    continue;
                }
                ism.iter().enumerate().for_each(|(n, &i)| {
                    if self.graph[cmp.node_map[n]].protons != 0 {
                        matched[i] = (k + 1, n);
                    }
                });
                found.push((*n, ism));
            }
        }

        let (ret, rng) = if found.is_empty() {
            // simple case: no subgraph isomorhpisms found
            let mut map = vec![NodeIndex::end(); mol.node_bound()];
            let mut bits = BSType::new();
            mol.node_references().for_each(|n| {
                let idx = self.graph.add_node(*n.weight());
                map[mol.to_index(n.id())] = idx;
                bits.set(idx.index(), true);
            });
            mol.edge_references().for_each(|e| {
                self.graph.add_edge(
                    map[mol.to_index(e.source())],
                    map[mol.to_index(e.target())],
                    *e.weight(),
                );
            });
            let out = self.parts.len();
            self.parts
                .push((MolRepr::Atomic(bits), Ix::new(mol.node_count())));
            (Ix::new(out), out..(out + 1))
        } else {
            let mut nb = BSType::with_capacity(mol.node_bound());
            let mut needs_clear = false;
            let mut frags = SmallVec::with_capacity(found.len() + 1);
            frags.push(Ix::new(self.parts.len()));
            frags.extend(found.iter().map(|(i, _)| Ix::new(*i)));
            let mut bonds = SmallVec::new();

            for (n, (i, ism)) in found.iter().enumerate() {
                let cmp = &compacted[*i].1;
                for (j, &mi) in ism.iter().enumerate() {
                    let idx = cmp.node_map[j];
                    let atom = self.graph[idx];
                    let ma = mol.node_weight(mol.from_index(mi)).unwrap();
                    if atom.protons == 0 && ma.protons != 0 {
                        let an = Ix::new(0);
                        let ai = Ix::new(mi);
                        let bn = Ix::new(n + 1);
                        let bi = Ix::new(idx.index());
                        bonds.push(InterFragBond { an, ai, bn, bi });
                    }
                    if atom.data.unknown() > 0 {
                        if needs_clear {
                            nb.clear();
                        }
                        needs_clear = true;
                        for c in self.graph.neighbors(idx) {
                            nb.set(c.index(), true);
                        }
                        #[cfg(debug_assertions)]
                        let mut matched_count = 0;

                        // look for atoms bonded to this atom that don't have a corresponding node
                        // in the fragment, those *must* be an R-group
                        bonds.extend(mol.neighbors(mol.from_index(mi)).filter_map(|mi| {
                            let mix = mol.to_index(mi);
                            if nb.get(mix) {
                                None // already accounted for
                            } else {
                                #[cfg(debug_assertions)]
                                {
                                    matched_count += 1;
                                }
                                Some(InterFragBond {
                                    an: Ix::new(0),
                                    ai: Ix::new(mix),
                                    bn: Ix::new(n + 1),
                                    bi: Ix::new(idx.index()),
                                })
                            }
                        }));

                        // #[cfg(debug_assertions)]
                        // assert_eq!(matched_count, atom.data.unknown());
                    }
                }
            }

            let rng = if matched.iter().all(|m| m.0 != 0) {
                // all atoms make up other molecules
                // TODO: maybe use `swap_remove` here?
                frags.remove(0);
                for bond in &mut bonds {
                    let (mut an, ai) = matched[bond.ai.index()];
                    an -= 1;
                    if an <= bond.bn.index() {
                        bond.an = Ix::new(an);
                        bond.ai = Ix::new(ai);
                    } else {
                        bond.an = std::mem::replace(&mut bond.bn, Ix::new(an));
                        bond.ai = std::mem::replace(&mut bond.bi, Ix::new(ai));
                    }
                }

                0..0
            } else {
                let mut map = vec![NodeIndex::end(); mol.node_bound()];
                let mut bits = BSType::new();
                let mut count = 0;
                mol.node_references().for_each(|n| {
                    let ni = mol.to_index(n.id());
                    if matched[ni].0.index() == 0 {
                        count += 1;
                        let idx = self.graph.add_node(*n.weight());
                        map[ni] = idx;
                        bits.set(idx.index(), true);
                    }
                });
                mol.edge_references().for_each(|e| {
                    let si = mol.to_index(e.source());
                    let ti = mol.to_index(e.target());
                    if matched[ti].0.index() != 0 || matched[si].0.index() != 0 {
                        return;
                    }
                    self.graph.add_edge(
                        map[mol.to_index(e.source())],
                        map[mol.to_index(e.target())],
                        *e.weight(),
                    );
                });

                for bond in &mut bonds {
                    let ix = map[bond.ai.index()].index();
                    if ix == <Ix as IndexType>::max().index() {
                        let (an, ai) = matched[bond.ai.index()];
                        if an <= bond.bn.index() {
                            bond.an = Ix::new(an);
                            bond.ai = Ix::new(ai);
                        } else {
                            bond.an = std::mem::replace(&mut bond.bn, Ix::new(an));
                            bond.ai = std::mem::replace(&mut bond.bi, Ix::new(ai));
                        }
                    } else {
                        bond.ai = Ix::new(ix);
                    }
                }

                let start = self.parts.len();

                self.parts.push((MolRepr::Atomic(bits), Ix::new(count)));

                // TODO: split the fragment if disjoint!

                start..self.parts.len()
            };

            let out = Ix::new(self.parts.len());
            self.parts.push((
                MolRepr::Broken(BrokenMol { frags, bonds }),
                Ix::new(mol.node_count()),
            ));

            (out, rng)
        };

        for new in rng {
            let Some((MolRepr::Atomic(a), _)) = self.parts.get(new) else {
                continue;
            };
            let nc = GraphCompactor::<BitFiltered<&Graph<Ix>, _, ATOM_BIT_STORAGE>>::new(
                BitFiltered::new(&self.graph, a.clone()),
            );
            for (n, cmp) in &compacted {
                if is_isomorphic_matching(cmp, &nc, &mut Atom::eq_match_r, &mut PartialEq::eq) {
                    self.parts[new].0 = MolRepr::Redirect(Ix::new(*n));
                    break;
                }
            }
        }

        ret
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
pub struct Molecule<Ix, R> {
    arena: R,
    index: Ix,
}
impl<Ix: IndexType, R> Molecule<Ix, R> {
    pub fn from_arena<'a, 'b: 'a, A: ArenaAccessible<Ix = Ix, Access<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor(),
            index,
        }
    }

    pub fn from_mut_arena<'a, 'b: 'a, A: ArenaAccessibleMut<Ix = Ix, AccessMut<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
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
        R: ArenaAccessor<Ix = Ix>,
    {
        self.arena.get_arena()
    }

    pub fn arena_mut(&self) -> R::RefMut<'_>
    where
        R: ArenaAccessorMut<Ix = Ix>,
    {
        self.arena.get_arena_mut()
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Molecule<Ix, R> {
    pub fn contains(&self, group: Ix) -> bool
    where
        R: ArenaAccessor<Ix = Ix>,
    {
        self.arena.get_arena().contains_group(self.index, group)
    }
}

/// This trait handles the access to the backing arena. Rather than just passing around references,
/// this allows for lock guards to be used while not forcing them to live for as long as the
/// accessor.
pub trait ArenaAccessor {
    type Ix: IndexType;
    type Ref<'a>: Deref<Target = Arena<Self::Ix>>
    where
        Self: 'a;

    fn get_arena<'a>(&'a self) -> Self::Ref<'a>
    where
        Self: 'a;
}

/// This trait provides mutable access to the underlying arena.
pub trait ArenaAccessorMut: ArenaAccessor {
    type RefMut<'a>: DerefMut<Target = Arena<Self::Ix>>
    where
        Self: 'a;

    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        Self: 'a;
}

/// This trait allows access to a backing arena.
pub trait ArenaAccessible {
    type Ix: IndexType;
    type Access<'a>: ArenaAccessor<Ix = Self::Ix> + 'a
    where
        Self: 'a;

    fn get_accessor(&self) -> Self::Access<'_>;
}

/// This trait also allows access to a graph through an internally mutable type.
pub trait ArenaAccessibleMut: ArenaAccessible {
    type AccessMut<'a>: ArenaAccessorMut<Ix = Self::Ix> + 'a
    where
        Self: 'a;

    fn get_accessor_mut(&self) -> Self::AccessMut<'_>;
}

impl<T: ArenaAccessible> ArenaAccessible for &T {
    type Ix = T::Ix;
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
pub struct PtrAcc<Ix: IndexType>(*mut Arena<Ix>);
impl<Ix: IndexType> ArenaAccessor for PtrAcc<Ix> {
    type Ix = Ix;
    type Ref<'a> = &'a Arena<Ix>;

    fn get_arena<'a>(&'a self) -> Self::Ref<'a>
    where
        Self: 'a,
    {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe { &*self.0 }
    }
}
impl<Ix: IndexType> ArenaAccessorMut for PtrAcc<Ix> {
    type RefMut<'a> = &'a mut Arena<Ix>;

    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        Self: 'a,
    {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe { &mut *self.0 }
    }
}

/// Wrapper type around a `*mut Arena`. It has an `unsafe` constructor because the
/// `ArenaAccessible` implementation can't be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArenaPtr<Ix: IndexType>(*mut Arena<Ix>);
impl<Ix: IndexType> ArenaPtr<Ix> {
    /// # Safety
    /// This type basically wraps a pointer, but moves the `unsafe` to its construction. All
    /// pointer invariants must hold, as they won't be checked elsewhere.
    pub const unsafe fn new(ptr: *mut Arena<Ix>) -> Self {
        Self(ptr)
    }
    pub fn get_ptr(self) -> *mut Arena<Ix> {
        self.0
    }
}
impl<Ix: IndexType> ArenaAccessible for ArenaPtr<Ix> {
    type Ix = Ix;
    type Access<'a> = PtrAcc<Ix>;

    fn get_accessor(&self) -> PtrAcc<Ix> {
        PtrAcc(self.0)
    }
}
impl<Ix: IndexType> ArenaAccessibleMut for ArenaPtr<Ix> {
    type AccessMut<'a> = PtrAcc<Ix>;

    fn get_accessor_mut(&self) -> PtrAcc<Ix> {
        PtrAcc(self.0)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RefAcc<'a, Ix: IndexType>(&'a Arena<Ix>);
impl<'a, Ix: IndexType> ArenaAccessor for RefAcc<'a, Ix> {
    type Ix = Ix;
    type Ref<'b> = &'b Arena<Ix> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0
    }
}

impl<Ix: IndexType> ArenaAccessible for Arena<Ix> {
    type Ix = Ix;
    type Access<'a> = RefAcc<'a, Ix>;

    fn get_accessor(&self) -> RefAcc<Ix> {
        RefAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RefCellAcc<'a, Ix: IndexType>(&'a RefCell<Arena<Ix>>);
impl<'a, Ix: IndexType> ArenaAccessor for RefCellAcc<'a, Ix> {
    type Ix = Ix;
    type Ref<'b> = Ref<'b, Arena<Ix>> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.borrow()
    }
}
impl<'a, Ix: IndexType> ArenaAccessorMut for RefCellAcc<'a, Ix> {
    type RefMut<'b> = RefMut<'b, Arena<Ix>> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.borrow_mut()
    }
}

impl<Ix: IndexType> ArenaAccessible for RefCell<Arena<Ix>> {
    type Ix = Ix;
    type Access<'a> = RefCellAcc<'a, Ix>;

    fn get_accessor(&self) -> RefCellAcc<Ix> {
        RefCellAcc(self)
    }
}
impl<Ix: IndexType> ArenaAccessibleMut for RefCell<Arena<Ix>> {
    type AccessMut<'a> = RefCellAcc<'a, Ix>;

    fn get_accessor_mut(&self) -> RefCellAcc<Ix> {
        RefCellAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct RwLockAcc<'a, Ix: IndexType>(&'a RwLock<Arena<Ix>>);
impl<'a, Ix: IndexType> ArenaAccessor for RwLockAcc<'a, Ix> {
    type Ix = Ix;
    type Ref<'b> = RwLockReadGuard<'b, Arena<Ix>> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.read()
    }
}
impl<'a, Ix: IndexType> ArenaAccessorMut for RwLockAcc<'a, Ix> {
    type RefMut<'b> = RwLockWriteGuard<'b, Arena<Ix>> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.write()
    }
}

impl<Ix: IndexType> ArenaAccessible for RwLock<Arena<Ix>> {
    type Ix = Ix;
    type Access<'a> = RwLockAcc<'a, Ix>;

    fn get_accessor(&self) -> RwLockAcc<Ix> {
        RwLockAcc(self)
    }
}
impl<Ix: IndexType> ArenaAccessibleMut for RwLock<Arena<Ix>> {
    type AccessMut<'a> = RwLockAcc<'a, Ix>;

    fn get_accessor_mut(&self) -> RwLockAcc<Ix> {
        RwLockAcc(self)
    }
}

impl<T: ArenaAccessible> ArenaAccessible for Rc<T> {
    type Ix = T::Ix;
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
    type Ix = T::Ix;
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
