use super::Arena;
use lock_api::{
    MappedRwLockReadGuard, MappedRwLockWriteGuard, RawRwLock, RwLock, RwLockReadGuard,
    RwLockWriteGuard,
};
use petgraph::graph::IndexType;
use std::cell::{Ref, RefCell, RefMut};
use std::ops::{Deref, DerefMut};
use std::ptr::NonNull;
use std::rc::Rc;
use std::sync::Arc;

/// This trait handles the access to the backing arena. Rather than just passing around references,
/// this allows for lock guards to be used while not forcing them to live for as long as the
/// accessor.
pub trait ArenaAccessor: Copy {
    type Ix: IndexType;
    type Data;
    type Ref<'a>: Deref<Target = Arena<Self::Ix, Self::Data>>
    where
        Self: 'a,
        Self::Data: 'a;
    type MappedRef<'a, T: 'a>: Deref<Target = T>
    where
        Self: 'a;
    fn get_arena<'a>(&'a self) -> Self::Ref<'a>
    where
        Self: 'a,
        Self::Data: 'a;

    fn map_ref<'a, T: 'a, F: FnOnce(&Arena<Self::Ix, Self::Data>) -> &T>(
        r: Self::Ref<'a>,
        f: F,
    ) -> Self::MappedRef<'a, T>;
    fn map_mapped_ref<'a, T: 'a, U: 'a, F: Fn(&T) -> &U>(
        r: Self::MappedRef<'a, T>,
        f: F,
    ) -> Self::MappedRef<'a, U>;
}

/// This trait provides mutable access to the underlying arena.
pub trait ArenaAccessorMut: ArenaAccessor {
    type RefMut<'a>: DerefMut<Target = Arena<Self::Ix, Self::Data>>
    where
        Self: 'a,
        Self::Data: 'a;
    type MappedRefMut<'a, T: 'a>: DerefMut<Target = T>
    where
        Self: 'a;

    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        Self: 'a,
        Self::Data: 'a;

    fn map_ref_mut<'a, T: 'a, F: FnOnce(&mut Arena<Self::Ix, Self::Data>) -> &mut T>(
        r: Self::RefMut<'a>,
        f: F,
    ) -> Self::MappedRefMut<'a, T>;
    fn map_mapped_ref_mut<'a, T: 'a, U: 'a, F: FnOnce(&mut T) -> &mut U>(
        r: Self::MappedRefMut<'a, T>,
        f: F,
    ) -> Self::MappedRefMut<'a, U>;
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
#[derive(Debug)]
#[repr(transparent)]
pub struct PtrAcc<Ix: IndexType, D>(NonNull<Arena<Ix, D>>);
impl<Ix: IndexType, D> Clone for PtrAcc<Ix, D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<Ix: IndexType, D> Copy for PtrAcc<Ix, D> {}
impl<Ix: IndexType, D> ArenaAccessor for PtrAcc<Ix, D> {
    type Ix = Ix;
    type Data = D;
    type Ref<'a> = &'a Arena<Ix, D> where D: 'a;
    type MappedRef<'a, T: 'a> = &'a T where D: 'a;

    fn get_arena<'a>(&'a self) -> Self::Ref<'a>
    where
        D: 'a,
    {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe { &mut *self.0.as_ptr() }
    }

    fn map_ref<'a, T: 'a, F: FnOnce(&Arena<Self::Ix, Self::Data>) -> &T>(
        r: Self::Ref<'a>,
        f: F,
    ) -> Self::MappedRef<'a, T> {
        f(r)
    }
    fn map_mapped_ref<'a, T: 'a, U: 'a, F: FnOnce(&T) -> &U>(
        r: Self::MappedRef<'a, T>,
        f: F,
    ) -> Self::MappedRef<'a, U>
    where
        D: 'a,
    {
        f(r)
    }
}
impl<Ix: IndexType, D> ArenaAccessorMut for PtrAcc<Ix, D> {
    type RefMut<'a> = &'a mut Arena<Ix, D> where D: 'a;
    type MappedRefMut<'a, T: 'a> = &'a mut T where D: 'a;

    fn get_arena_mut<'a>(&'a self) -> Self::RefMut<'a>
    where
        D: 'a,
    {
        // Safety: this was created with an `unsafe`, so the caller must know about the invariants
        unsafe { &mut *self.0.as_ptr() }
    }

    fn map_ref_mut<'a, T: 'a, F: FnOnce(&mut Arena<Self::Ix, Self::Data>) -> &mut T>(
        r: Self::RefMut<'a>,
        f: F,
    ) -> Self::MappedRefMut<'a, T> {
        f(r)
    }
    fn map_mapped_ref_mut<'a, T: 'a, U: 'a, F: FnOnce(&mut T) -> &mut U>(
        r: Self::MappedRefMut<'a, T>,
        f: F,
    ) -> Self::MappedRefMut<'a, U>
    where
        D: 'a,
    {
        f(r)
    }
}

/// Wrapper type around a `*mut Arena`. It has an `unsafe` constructor because the
/// `ArenaAccessible` implementation can't be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct ArenaPtr<Ix: IndexType, D>(NonNull<Arena<Ix, D>>);
impl<Ix: IndexType, D> ArenaPtr<Ix, D> {
    /// # Safety
    /// This type basically wraps a pointer, but moves the `unsafe` to its construction. All
    /// pointer invariants must hold, as they won't be checked elsewhere.
    pub const unsafe fn new(ptr: NonNull<Arena<Ix, D>>) -> Self {
        Self(ptr)
    }
    pub fn get_ptr(self) -> *mut Arena<Ix, D> {
        self.0.as_ptr()
    }
}
impl<Ix: IndexType, D> ArenaAccessible for ArenaPtr<Ix, D> {
    type Ix = Ix;
    type Access<'a> = PtrAcc<Ix, D> where D: 'a;

    fn get_accessor(&self) -> PtrAcc<Ix, D> {
        PtrAcc(self.0)
    }
}
impl<Ix: IndexType, D> ArenaAccessibleMut for ArenaPtr<Ix, D> {
    type AccessMut<'a> = PtrAcc<Ix, D> where D: 'a;

    fn get_accessor_mut(&self) -> PtrAcc<Ix, D> {
        PtrAcc(self.0)
    }
}

#[doc(hidden)]
#[derive(Debug)]
pub struct RefAcc<'a, Ix: IndexType, D>(&'a Arena<Ix, D>);
impl<Ix: IndexType, D> Clone for RefAcc<'_, Ix, D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<Ix: IndexType, D> Copy for RefAcc<'_, Ix, D> {}
impl<'a, Ix: IndexType, D> ArenaAccessor for RefAcc<'a, Ix, D> {
    type Ix = Ix;
    type Data = D;
    type Ref<'b> = &'b Arena<Ix, D> where 'a: 'b;
    type MappedRef<'b, T: 'b> = &'b T where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0
    }

    fn map_ref<'b, T: 'b, F: FnOnce(&Arena<Self::Ix, Self::Data>) -> &T>(
        r: Self::Ref<'b>,
        f: F,
    ) -> Self::MappedRef<'b, T>
    where
        'a: 'b,
    {
        f(r)
    }
    fn map_mapped_ref<'b, T: 'b, U: 'b, F: FnOnce(&T) -> &U>(
        r: Self::MappedRef<'b, T>,
        f: F,
    ) -> Self::MappedRef<'b, U>
    where
        'a: 'b,
    {
        f(r)
    }
}

impl<Ix: IndexType, D> ArenaAccessible for Arena<Ix, D> {
    type Ix = Ix;
    type Access<'a> = RefAcc<'a, Ix, D> where D: 'a;

    fn get_accessor(&self) -> RefAcc<Ix, D> {
        RefAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug)]
pub struct RefCellAcc<'a, Ix: IndexType, D>(&'a RefCell<Arena<Ix, D>>);
impl<Ix: IndexType, D> Clone for RefCellAcc<'_, Ix, D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<Ix: IndexType, D> Copy for RefCellAcc<'_, Ix, D> {}
impl<'a, Ix: IndexType, D> ArenaAccessor for RefCellAcc<'a, Ix, D> {
    type Ix = Ix;
    type Data = D;
    type Ref<'b> = Ref<'b, Arena<Ix, D>> where 'a: 'b;
    type MappedRef<'b, T: 'b> = Ref<'b, T> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.borrow()
    }

    fn map_ref<'b, T: 'b, F: FnOnce(&Arena<Self::Ix, D>) -> &T>(
        r: Self::Ref<'b>,
        f: F,
    ) -> Self::MappedRef<'b, T>
    where
        'a: 'b,
    {
        Ref::map(r, f)
    }
    fn map_mapped_ref<'b, T: 'b, U: 'b, F: FnOnce(&T) -> &U>(
        r: Self::MappedRef<'b, T>,
        f: F,
    ) -> Self::MappedRef<'b, U>
    where
        'a: 'b,
    {
        Ref::map(r, f)
    }
}
impl<'a, Ix: IndexType, D> ArenaAccessorMut for RefCellAcc<'a, Ix, D> {
    type RefMut<'b> = RefMut<'b, Arena<Ix, D>> where 'a: 'b;
    type MappedRefMut<'b, T: 'b> = RefMut<'b, T> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.borrow_mut()
    }

    fn map_ref_mut<'b, T: 'b, F: FnOnce(&mut Arena<Self::Ix, D>) -> &mut T>(
        r: Self::RefMut<'b>,
        f: F,
    ) -> Self::MappedRefMut<'b, T>
    where
        'a: 'b,
    {
        RefMut::map(r, f)
    }
    fn map_mapped_ref_mut<'b, T: 'b, U: 'b, F: FnOnce(&mut T) -> &mut U>(
        r: Self::MappedRefMut<'b, T>,
        f: F,
    ) -> Self::MappedRefMut<'b, U>
    where
        'a: 'b,
    {
        RefMut::map(r, f)
    }
}

impl<Ix: IndexType, D> ArenaAccessible for RefCell<Arena<Ix, D>> {
    type Ix = Ix;
    type Access<'a> = RefCellAcc<'a, Ix, D> where D: 'a;

    fn get_accessor(&self) -> RefCellAcc<Ix, D> {
        RefCellAcc(self)
    }
}
impl<Ix: IndexType, D> ArenaAccessibleMut for RefCell<Arena<Ix, D>> {
    type AccessMut<'a> = RefCellAcc<'a, Ix, D> where D: 'a;

    fn get_accessor_mut(&self) -> RefCellAcc<Ix, D> {
        RefCellAcc(self)
    }
}

#[doc(hidden)]
#[derive(Debug)]
pub struct RwLockAcc<'a, Ix: IndexType, D, R: RawRwLock + 'a>(&'a RwLock<R, Arena<Ix, D>>);
// Manual impls of `Clone` and `Copy` because the derives add bounds to `R`
impl<'a, Ix: IndexType, D, R: RawRwLock + 'a> Clone for RwLockAcc<'a, Ix, D, R> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a, Ix: IndexType, D, R: RawRwLock + 'a> Copy for RwLockAcc<'a, Ix, D, R> {}
impl<'a, Ix: IndexType, D, R: RawRwLock + 'a> ArenaAccessor for RwLockAcc<'a, Ix, D, R> {
    type Ix = Ix;
    type Data = D;
    type Ref<'b> = RwLockReadGuard<'b, R, Arena<Ix, D>> where 'a: 'b, D: 'b;
    type MappedRef<'b, T: 'b> = MappedRwLockReadGuard<'b, R, T> where 'a: 'b;

    fn get_arena<'b>(&'b self) -> Self::Ref<'b>
    where
        'a: 'b,
    {
        self.0.read()
    }

    fn map_ref<'b, T: 'b, F: FnOnce(&Arena<Self::Ix, D>) -> &T>(
        r: Self::Ref<'b>,
        f: F,
    ) -> Self::MappedRef<'b, T>
    where
        'a: 'b,
    {
        RwLockReadGuard::map(r, f)
    }
    fn map_mapped_ref<'b, T: 'b, U: 'b, F: FnOnce(&T) -> &U>(
        r: Self::MappedRef<'b, T>,
        f: F,
    ) -> Self::MappedRef<'b, U>
    where
        'a: 'b,
    {
        MappedRwLockReadGuard::map(r, f)
    }
}
impl<'a, Ix: IndexType, D, R: RawRwLock + 'a> ArenaAccessorMut for RwLockAcc<'a, Ix, D, R> {
    type RefMut<'b> = RwLockWriteGuard<'b, R, Arena<Ix, D>> where 'a: 'b;
    type MappedRefMut<'b, T: 'b> = MappedRwLockWriteGuard<'b, R, T> where 'a: 'b;

    fn get_arena_mut<'b>(&'b self) -> Self::RefMut<'b>
    where
        'a: 'b,
    {
        self.0.write()
    }

    fn map_ref_mut<'b, T: 'b, F: FnOnce(&mut Arena<Self::Ix, D>) -> &mut T>(
        r: Self::RefMut<'b>,
        f: F,
    ) -> Self::MappedRefMut<'b, T>
    where
        'a: 'b,
    {
        RwLockWriteGuard::map(r, f)
    }
    fn map_mapped_ref_mut<'b, T: 'b, U: 'b, F: FnOnce(&mut T) -> &mut U>(
        r: Self::MappedRefMut<'b, T>,
        f: F,
    ) -> Self::MappedRefMut<'b, U>
    where
        'a: 'b,
    {
        MappedRwLockWriteGuard::map(r, f)
    }
}

impl<Ix: IndexType, D, R: RawRwLock> ArenaAccessible for RwLock<R, Arena<Ix, D>> {
    type Ix = Ix;
    type Access<'a> = RwLockAcc<'a, Ix, D, R> where R: 'a, D: 'a;

    fn get_accessor(&self) -> RwLockAcc<Ix, D, R> {
        RwLockAcc(self)
    }
}
impl<Ix: IndexType, D, R: RawRwLock> ArenaAccessibleMut for RwLock<R, Arena<Ix, D>> {
    type AccessMut<'a> = RwLockAcc<'a, Ix, D, R> where R: 'a, D: 'a;

    fn get_accessor_mut(&self) -> RwLockAcc<Ix, D, R> {
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
