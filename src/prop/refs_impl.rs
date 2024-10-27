use super::*;
use std::rc::Rc;
use std::sync::Arc;

impl<T, I, C: Property<T, I>> Property<T, Inside<I>> for &C {
    type Ref<'a>
        = C::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyMut<T, Inside<I>> for &C {
    type RefMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_imut()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyIMut<T, Inside<I>> for &C {
    type RefIMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        (**self).get_prop_imut()
    }
}

impl<T, I, C: Property<T, I>> Property<T, Inside<I>> for &mut C {
    type Ref<'a>
        = C::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, C: PropertyMut<T, I>> PropertyMut<T, Inside<I>> for &mut C {
    type RefMut<'a>
        = C::RefMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_mut()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyIMut<T, Inside<I>> for &mut C {
    type RefIMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        (**self).get_prop_imut()
    }
}

impl<T, I, C: Property<T, I>> Property<T, Inside<I>> for Box<C> {
    type Ref<'a>
        = C::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, C: PropertyMut<T, I>> PropertyMut<T, Inside<I>> for Box<C> {
    type RefMut<'a>
        = C::RefMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_mut()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyIMut<T, Inside<I>> for Box<C> {
    type RefIMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        (**self).get_prop_imut()
    }
}

impl<T, I, C: Property<T, I>> Property<T, Inside<I>> for Rc<C> {
    type Ref<'a>
        = C::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyMut<T, Inside<I>> for Rc<C> {
    type RefMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_imut()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyIMut<T, Inside<I>> for Rc<C> {
    type RefIMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        (**self).get_prop_imut()
    }
}

impl<T, I, C: Property<T, I>> Property<T, Inside<I>> for Arc<C> {
    type Ref<'a>
        = C::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyMut<T, Inside<I>> for Arc<C> {
    type RefMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_imut()
    }
}
impl<T, I, C: PropertyIMut<T, I>> PropertyIMut<T, Inside<I>> for Arc<C> {
    type RefIMut<'a>
        = C::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        (**self).get_prop_imut()
    }
}
