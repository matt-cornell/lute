//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use crate::molecule::{Atom, Bond};
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::stable_graph::StableGraph;
use std::cell::{Ref, RefCell, RefMut};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Clone)]
pub struct Arena {
    graph: StableGraph<Atom, Bond>,
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
