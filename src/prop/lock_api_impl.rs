use super::*;
use lock_api::*;

impl<T, I, C: PropertyRef<T, I>, R: RawRwLock> Property<T, Inside<I>> for RwLock<R, C> {
    type Ref<'a> = MappedRwLockReadGuard<'a, R, T> where Self: 'a, T: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        RwLockReadGuard::map(self.read(), |r| r.get_prop_ref())
    }
}
impl<T, I, C: PropertyMut<T, I> + PropertyRef<T, I>, R: RawRwLock> PropertyMut<T, Inside<I>> for RwLock<R, C> {
    type RefMut<'a> = C::RefMut<'a> where Self: 'a, T: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.get_mut().get_prop_mut()
    }
}
impl<T, I, C: PropertyMutRef<T, I> + PropertyRef<T, I>, R: RawRwLock> PropertyIMut<T, Inside<I>> for RwLock<R, C> {
    type RefIMut<'a> = MappedRwLockWriteGuard<'a, R, T> where Self: 'a, T: 'a;

    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        RwLockWriteGuard::map(self.write(), |r| r.get_prop_mut_ref())
    }
}
