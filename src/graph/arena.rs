//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use crate::molecule::{Atom, Bond, Chirality};
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::prelude::*;
use petgraph::visit::*;
use std::cell::{Ref, RefCell, RefMut};
use std::io::{self, Read, Write};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

type Ix = petgraph::graph::DefaultIx;
const IX_SIZE: usize = std::mem::size_of::<Ix>();

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena {
    graph: StableUnGraph<Atom, Bond>,
}
impl Arena {
    pub fn new() -> Self {
        Self::default()
    }

    /// Load a previously compiled graph. Note that this doesn't perform additional verification,
    /// it assumes that the previously compiled graph was valid.
    pub fn from_compiled<R: Read>(mut buf: R) -> io::Result<Self> {
        let mut ibuf = [0u8; IX_SIZE];
        let mut abuf = [0u8; 5];
        buf.read_exact(&mut ibuf)?;
        let node_count = Ix::from_le_bytes(ibuf);
        buf.read_exact(&mut ibuf)?;
        let edge_count = Ix::from_le_bytes(ibuf);
        let mut graph = StableGraph::with_capacity(node_count as _, edge_count as _);
        for _ in 0..node_count {
            buf.read_exact(&mut abuf)?;
            let protons = abuf[0];
            let charge = unsafe { *(&abuf[1] as *const u8 as *const i8) };
            let chirality = match abuf[2] {
                0 => Chirality::None,
                1 => Chirality::Ccw,
                2 => Chirality::Cw,
                c => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected 0, 1, or 2 for chirality kind, found {c}"),
                ))?,
            };
            let mut isotope = Some(((abuf[4] as u16) << 8) | (abuf[3] as u16));
            if isotope == Some(0) && protons != 0 {
                isotope = None;
            }
            graph.add_node(Atom {
                protons,
                charge,
                chirality,
                isotope,
                aromatic: false,
            });
        }
        for _ in 0..edge_count {
            buf.read_exact(&mut ibuf)?;
            let source = NodeIndex::new(Ix::from_le_bytes(ibuf) as _);
            buf.read_exact(&mut ibuf)?;
            let target = NodeIndex::new(Ix::from_le_bytes(ibuf) as _);
            let mut c = 0u8;
            buf.read_exact(std::slice::from_mut(&mut c))?;
            let weight = match c {
                0 => Bond::Non,
                1 => Bond::Single,
                2 => Bond::Double,
                3 => Bond::Triple,
                4 => Bond::Quad,
                5 => Bond::Aromatic,
                6 => Bond::Left,
                7 => Bond::Right,
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected 0..=7 for bond kind, found {c}"),
                ))?,
            };
            graph.add_edge(source, target, weight);
        }
        Ok(Self { graph })
    }
    /// Return the size of the binary output in bytes
    pub fn comp_size(&self) -> usize {
        const NODE_SIZE: usize = 5;
        const EDGE_SIZE: usize = 2 * IX_SIZE + 1;
        2 * IX_SIZE + self.graph.node_count() * NODE_SIZE + self.graph.edge_count() * EDGE_SIZE
    }
    /// Create a binary output from which this graph can be reconstructed
    pub fn compile<W: Write>(&self, mut buf: W) -> io::Result<()> {
        let mut ebuf = [0u8; 2 * IX_SIZE + 1];
        ebuf[0..IX_SIZE].copy_from_slice(&(self.graph.node_count() as Ix).to_le_bytes());
        ebuf[IX_SIZE..(2 * IX_SIZE)].copy_from_slice(&(self.graph.edge_count() as Ix).to_le_bytes());
        buf.write_all(&ebuf[..16])?;
        let mut table = vec![None; self.graph.node_count()];
        let mut c: Ix = 0;
        for node in self.graph.node_indices() {
            table[node.index()] = Some(c);
            let atom = &self.graph[node];
            buf.write_all(&[
                atom.protons,
                unsafe { *(&atom.charge as *const i8 as *const u8) },
                atom.chirality as u8,
                atom.isotope.map_or(0, |i| (i & 255) as u8),
                atom.isotope.map_or(0, |i| (i >> 8) as u8),
            ])?;
            c += 1;
        }
        for edge in self.graph.edge_references() {
            ebuf[0..IX_SIZE].copy_from_slice(&table[edge.source().index()].unwrap().to_le_bytes());
            ebuf[IX_SIZE..(2 * IX_SIZE)].copy_from_slice(&table[edge.target().index()].unwrap().to_le_bytes());
            ebuf[2 * IX_SIZE] = unsafe { *(edge.weight() as *const _ as *const u8) };
            buf.write_all(&ebuf)?;
        }
        Ok(())
    }
}

/// A `Container` corresponds to a reaction vessel. It's at this layer that actual reactions are
/// handled.
#[derive(Debug, Clone, Copy)]
pub struct Container {}

/// A `Molecule` is a graph, and can have graph algorithms used on it. It's immutable, with all
/// mutationss making (efficient) copies.
#[derive(Debug, Clone, Copy)]
pub struct Molecule<R> {
    arena: R,
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

    /// ## Safety
    /// While in most cases, this is safe, it's up to the caller that multiple accesses never take
    /// place
    unsafe fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        Self: 'a;
}

/// This trait allows access to a backing arena.
pub trait ArenaAccessible {
    type Access<'a>: ArenaAccessor + 'a
    where
        Self: 'a;
    type AccessMut<'a>: ArenaAccessorMut + 'a
    where
        Self: 'a;

    fn get_accessor(&self) -> Self::Access<'_>;
    fn get_accessor_mut(&mut self) -> Self::AccessMut<'_>;
}

/// This trait also allows access to a graph through an internally mutable type.
pub trait ArenaAccessibleMut: ArenaAccessible {
    type InternalAccess<'a>: ArenaAccessorMut + 'a
    where
        Self: 'a;

    fn get_accessor_imut(&self) -> Self::InternalAccess<'_>;
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

// `PhantomData` to impl `!Sync`. Use a mutex if you want to use this across threads.
#[doc(hidden)]
#[derive(Debug)]
pub struct RefMutAcc<'a>(&'a mut Arena, PhantomData<RefCell<()>>);

impl<'a> ArenaAccessor for RefMutAcc<'a> {
    type Ref<'b> = &'b Arena where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0
    }
}
impl<'a> ArenaAccessorMut for RefMutAcc<'a> {
    type RefMut<'b> = &'b mut Arena where 'a: 'b;

    unsafe fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        std::ptr::read(&self.0)
    }
}

impl ArenaAccessible for Arena {
    type Access<'a> = RefAcc<'a>;
    type AccessMut<'a> = RefMutAcc<'a>;

    fn get_accessor(&self) -> RefAcc {
        RefAcc(self)
    }
    fn get_accessor_mut(&mut self) -> RefMutAcc {
        RefMutAcc(self, PhantomData)
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

    unsafe fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.borrow_mut()
    }
}

impl ArenaAccessible for RefCell<Arena> {
    type Access<'a> = RefCellAcc<'a>;
    type AccessMut<'a> = RefMutAcc<'a>;

    fn get_accessor(&self) -> RefCellAcc {
        RefCellAcc(self)
    }
    fn get_accessor_mut(&mut self) -> RefMutAcc {
        RefMutAcc(self.get_mut(), PhantomData)
    }
}
impl ArenaAccessibleMut for RefCell<Arena> {
    type InternalAccess<'a> = RefCellAcc<'a>;

    fn get_accessor_imut(&self) -> RefCellAcc {
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

    unsafe fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.write()
    }
}

impl ArenaAccessible for RwLock<Arena> {
    type Access<'a> = RwLockAcc<'a>;
    type AccessMut<'a> = RefMutAcc<'a>;

    fn get_accessor(&self) -> RwLockAcc {
        RwLockAcc(self)
    }
    fn get_accessor_mut(&mut self) -> RefMutAcc {
        RefMutAcc(self.get_mut(), PhantomData)
    }
}
impl ArenaAccessibleMut for RwLock<Arena> {
    type InternalAccess<'a> = RwLockAcc<'a>;

    fn get_accessor_imut(&self) -> RwLockAcc {
        RwLockAcc(self)
    }
}

impl<T: ArenaAccessibleMut> ArenaAccessible for Rc<T> {
    type Access<'a> = T::Access<'a> where T: 'a;
    type AccessMut<'a> = T::InternalAccess<'a> where T: 'a;

    fn get_accessor(&self) -> Self::Access<'_> {
        T::get_accessor(self)
    }
    fn get_accessor_mut(&mut self) -> Self::AccessMut<'_> {
        T::get_accessor_imut(&**self)
    }
}
impl<T: ArenaAccessibleMut> ArenaAccessibleMut for Rc<T> {
    type InternalAccess<'a> = T::InternalAccess<'a> where T: 'a;

    fn get_accessor_imut(&self) -> T::InternalAccess<'_> {
        T::get_accessor_imut(&**self)
    }
}

impl<T: ArenaAccessibleMut> ArenaAccessible for Arc<T> {
    type Access<'a> = T::Access<'a> where T: 'a;
    type AccessMut<'a> = T::InternalAccess<'a> where T: 'a;

    fn get_accessor(&self) -> Self::Access<'_> {
        T::get_accessor(self)
    }
    fn get_accessor_mut(&mut self) -> Self::AccessMut<'_> {
        T::get_accessor_imut(&**self)
    }
}
impl<T: ArenaAccessibleMut> ArenaAccessibleMut for Arc<T> {
    type InternalAccess<'a> = T::InternalAccess<'a> where T: 'a;

    fn get_accessor_imut(&self) -> T::InternalAccess<'_> {
        T::get_accessor_imut(&**self)
    }
}
