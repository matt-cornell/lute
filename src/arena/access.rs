use super::Arena;
use parking_lot::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use petgraph::graph::IndexType;
use std::cell::{Ref, RefCell, RefMut};
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use std::sync::Arc;

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
