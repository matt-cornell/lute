//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use crate::graph::compact::GraphCompactor;
use crate::molecule::{Atom, Bond};
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::algo::subgraph_isomorphisms_iter;
use petgraph::data::DataMap;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallbitvec::{InternalStorage, SmallBitVec};
use smallvec::SmallVec;
use std::cell::{Ref, RefCell, RefMut};
use std::collections::HashMap;
use std::io::{self, Read, Write};
use std::mem::ManuallyDrop;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

type Ix = petgraph::graph::DefaultIx;
const IX_SIZE: usize = std::mem::size_of::<Ix>();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
struct InterFragBond {
    ai: Ix,
    ar: Ix,
    bi: Ix,
    br: Ix,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BrokenMol {
    frags: SmallVec<[Ix; 8]>,
    bonds: SmallVec<[InterFragBond; 8]>,
}
impl BrokenMol {
    pub fn comp_size(&self) -> usize {
        (2 + self.frags.len() + self.bonds.len() * 4) * IX_SIZE
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MolRepr {
    Atomic(SmallBitVec),
    Broken(BrokenMol),
}
impl MolRepr {
    pub fn comp_size(&self) -> usize {
        const US_SIZE: usize = std::mem::size_of::<usize>() * 8;
        match self {
            Self::Atomic(a) => IX_SIZE + (a.len() * (US_SIZE - 1)) / US_SIZE,
            Self::Broken(b) => b.comp_size(),
        }
    }
}

/// Write a slice of POD data to a writer
fn write_pod_slice<T: Copy, W: Write>(data: &[T], buf: &mut W) -> io::Result<()> {
    // Safety: we're reinterpreting as bytes here, which will always be valid
    unsafe {
        let ptr = data.as_ptr() as *const u8;
        let len = std::mem::size_of_val(data);
        buf.write_all(&len.to_ne_bytes())?;
        buf.write_all(std::slice::from_raw_parts(ptr, len))?;
        Ok(())
    }
}
/// Read a slice of data from a reader, reinterpreting the bytes
/// # Safety
/// `T` must be valid for any memory contents
unsafe fn read_pod_slice<T: Copy, R: Read>(buf: &mut R) -> io::Result<Vec<T>> {
    let mut arr = [0u8; std::mem::size_of::<usize>()];
    buf.read_exact(&mut arr)?;
    let len = usize::from_ne_bytes(arr);
    let mut mdv = ManuallyDrop::new(Vec::<T>::with_capacity(len));
    let ptr = mdv.as_mut_ptr();
    let cap = mdv.capacity();
    assert!(cap >= len);
    buf.read_exact(std::slice::from_raw_parts_mut(ptr as *mut _, len))?;
    Ok(Vec::from_raw_parts(ptr, len, cap))
}
/// Read a slice of data from a reader, reinterpreting the bytes
/// # Safety
/// `T` must be valid for any memory contents
unsafe fn read_pod_slice_with_buf<T: Copy, R: Read>(
    buf: &mut R,
    vec: &mut Vec<T>,
) -> io::Result<()> {
    let mut arr = [0u8; std::mem::size_of::<usize>()];
    buf.read_exact(&mut arr)?;
    let len = usize::from_ne_bytes(arr);
    vec.clear();
    vec.reserve(len);
    let mut mdv = ManuallyDrop::new(std::mem::take(vec));
    let ptr = mdv.as_mut_ptr();
    let cap = mdv.capacity();
    assert!(cap >= len);
    buf.read_exact(std::slice::from_raw_parts_mut(ptr as *mut _, len))?;
    *vec = Vec::from_raw_parts(ptr, len, cap);
    Ok(())
}

#[derive(Clone, Copy)]
#[repr(C)]
struct EdgeRepr {
    source: Ix,
    target: Ix,
    weight: Bond,
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena {
    graph: StableUnGraph<Atom, Bond>,
    parts: SmallVec<[MolRepr; 16]>,
}
impl Arena {
    pub fn new() -> Self {
        Self::default()
    }

    /// Load a previously compiled graph. Note that this doesn't perform additional verification,
    /// it assumes that the previously compiled graph was valid.
    pub fn from_compiled<R: Read>(mut buf: R) -> io::Result<Self> {
        let mut ibuf = [0u8; IX_SIZE];
        let mut abuf = Vec::with_capacity(1);
        let mut fbuf = Vec::new();
        let mut bbuf = Vec::new();
        unsafe {
            let nodes = read_pod_slice::<Atom, R>(&mut buf)?;
            let edges = read_pod_slice::<EdgeRepr, R>(&mut buf)?;
            let mut graph = StableGraph::with_capacity(nodes.len() as _, edges.len() as _);
            for node in nodes {
                graph.add_node(node);
            }
            for EdgeRepr {
                source,
                target,
                weight,
            } in edges
            {
                graph.add_edge(source.into(), target.into(), weight);
            }
            buf.read_exact(&mut ibuf)?;
            let part_count = Ix::from_ne_bytes(ibuf) as usize;
            let mut parts = SmallVec::with_capacity(part_count);
            let mut i = 0u8;
            for _ in 0..part_count {
                buf.read_exact(std::slice::from_mut(&mut i))?;
                if i == 0 {
                    read_pod_slice_with_buf(&mut buf, &mut abuf)?;
                    let is = if abuf.len() == 1 {
                        InternalStorage::Inline(abuf[0])
                    } else {
                        InternalStorage::Spilled(Box::from(abuf.as_slice()))
                    };
                    parts.push(MolRepr::Atomic(SmallBitVec::from_storage(is)));
                } else {
                    read_pod_slice_with_buf(&mut buf, &mut fbuf)?;
                    read_pod_slice_with_buf(&mut buf, &mut bbuf)?;
                    parts.push(MolRepr::Broken(BrokenMol {
                        frags: SmallVec::from_slice(&fbuf),
                        bonds: SmallVec::from_slice(&bbuf),
                    }));
                }
            }
            Ok(Self { graph, parts })
        }
    }

    /// Return the size of the binary output in bytes
    pub fn comp_size(&self) -> usize {
        const NODE_SIZE: usize = 5;
        const EDGE_SIZE: usize = 2 * IX_SIZE + 1;
        3 * IX_SIZE
            + self.graph.node_count() * NODE_SIZE
            + self.graph.edge_count() * EDGE_SIZE
            + self.parts.iter().map(MolRepr::comp_size).sum::<usize>()
    }

    /// Create a binary output from which this graph can be reconstructed
    pub fn compile<W: Write>(&self, mut buf: W) -> io::Result<()> {
        let mut table = vec![None; self.graph.node_count()];
        let mut c: Ix = 0;
        let nodes = self
            .graph
            .node_indices()
            .map(|i| {
                table[i.index()] = Some(c);
                c += 1;
                self.graph[i]
            })
            .collect::<Vec<_>>();
        let edges = self
            .graph
            .edge_references()
            .map(|e| EdgeRepr {
                source: table[e.source().index()].unwrap(),
                target: table[e.target().index()].unwrap(),
                weight: *e.weight(),
            })
            .collect::<Vec<_>>();
        write_pod_slice(&nodes, &mut buf)?;
        write_pod_slice(&edges, &mut buf)?;
        buf.write_all(&(self.parts.len() as Ix).to_ne_bytes())?;
        for part in &self.parts {
            match part {
                MolRepr::Atomic(a) => {
                    let i = a.clone().into_storage();
                    let s = match &i {
                        InternalStorage::Inline(s) => std::slice::from_ref(s),
                        InternalStorage::Spilled(s) => &**s,
                    };
                    buf.write_all(&[0])?;
                    write_pod_slice(s, &mut buf)?;
                }
                MolRepr::Broken(b) => {
                    buf.write_all(&[1])?;
                    write_pod_slice(&b.frags, &mut buf)?;
                    write_pod_slice(&b.bonds, &mut buf)?;
                }
            }
        }
        Ok(())
    }

    fn contains_group_impl(&self, mol: usize, group: usize, seen: &mut SmallBitVec) -> bool {
        if mol == group {
            return true;
        }
        if seen[mol] {
            return false;
        }
        seen.set(mol, true);
        if let Some(MolRepr::Broken(b)) = self.parts.get(mol) {
            b.frags
                .iter()
                .any(|f| self.contains_group_impl(*f as _, group, seen))
        } else {
            false
        }
    }

    /// Check if `mol` contains `group`
    pub fn contains_group(&self, mol: usize, group: usize) -> bool {
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
            + EdgeCount
            + IntoEdgesDirected,
    {
        let compacted = (0..self.parts.len())
            .filter_map(|i| {
                if let MolRepr::Atomic(a) = &self.parts[i] {
                    Some(GraphCompactor::new(NodeFiltered::from_fn(
                        &self.graph,
                        |i| a[i.index()],
                    )))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        let mut broken = smallvec::smallvec![self.parts.len() as Ix];
        let mut bonds = SmallVec::new();
        for (n, cmp) in compacted.iter().enumerate() {
            for ism in
                subgraph_isomorphisms_iter(&cmp, &mol, &mut Atom::eq_or_r, &mut PartialEq::eq)
                    .into_iter()
                    .flatten()
            {}
        }
        todo!()
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
