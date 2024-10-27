use super::*;
use std::cell::{RefCell, Ref, RefMut};
use std::sync::*;

impl<T, I, C: PropertyRef<T, I>> Property<T, Inside<I>> for RefCell<C> {
    type Ref<'a> = Ref<'a, T> where Self: 'a, T: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        Ref::map(self.borrow(), |r| r.get_prop_ref())
    }
}
impl<T, I, C: PropertyMut<T, I> + PropertyRef<T, I>> PropertyMut<T, Inside<I>> for RefCell<C> {
    type RefMut<'a> = C::RefMut<'a> where Self: 'a, T: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.get_mut().get_prop_mut()
    }
}
impl<T, I, C: PropertyMutRef<T, I> + PropertyRef<T, I>> PropertyIMut<T, Inside<I>> for RefCell<C> {
    type RefIMut<'a> = RefMut<'a, T> where Self: 'a, T: 'a;

    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        RefMut::map(self.borrow_mut(), |r| r.get_prop_mut_ref())
    }
}

impl<T, I, C: PropertyRef<T, I>> Property<T, Inside<I>> for RwLock<C> {
    type Ref<'a> = MappedRwLockReadGuard<'a, T> where Self: 'a, T: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        RwLockReadGuard::map(self.read().unwrap(), |r| r.get_prop_ref())
    }
}
impl<T, I, C: PropertyMut<T, I> + PropertyRef<T, I>> PropertyMut<T, Inside<I>> for RwLock<C> {
    type RefMut<'a> = C::RefMut<'a> where Self: 'a, T: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.get_mut().unwrap().get_prop_mut()
    }
}
impl<T, I, C: PropertyMutRef<T, I> + PropertyRef<T, I>> PropertyIMut<T, Inside<I>> for RwLock<C> {
    type RefIMut<'a> = MappedRwLockWriteGuard<'a, T> where Self: 'a, T: 'a;

    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        RwLockWriteGuard::map(self.write().unwrap(), |r| r.get_prop_mut_ref())
    }
}
