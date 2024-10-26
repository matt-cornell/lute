// use crate::arena::{ArenaAccessor, ArenaAccessorMut, Molecule};
use frunk::hlist;
use frunk::HCons;
use indices::*;
// use petgraph::graph::IndexType;
use std::ops::{Deref, DerefMut};

/// Marker indices for recursive property access
pub mod indices {
    use std::marker::PhantomData;
    /// This current value is the property we want
    pub struct Current;
    /// The next value in the `HList` has the property we want
    pub struct Next<T> {
        _marker: PhantomData<T>,
    }
    /// The current value has the property we want
    pub struct Inside<T> {
        _marker: PhantomData<T>,
    }
    /// The derefed value has the property we want
    pub struct Derefed<T> {
        _marker: PhantomData<T>,
    }
}
pub trait Property<T, I> {
    /// Quack quack! I need this for the generic `Molecule` to be able to have properties
    type Ref<'a>: Deref<Target = T> + 'a
    where
        Self: 'a;

    fn get_prop(&self) -> Self::Ref<'_>;
}
pub trait PropertyMut<T, I>: Property<T, I> {
    /// Quack quack! I need this for the generic `Molecule` to be able to have properties
    type RefMut<'a>: DerefMut<Target = T> + 'a
    where
        Self: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_>;
}

impl<T> Property<T, Current> for T {
    type Ref<'a>
        = &'a T
    where
        Self: 'a;

    fn get_prop(&self) -> &T {
        self
    }
}
impl<T> PropertyMut<T, Current> for T {
    type RefMut<'a>
        = &'a mut T
    where
        Self: 'a;

    fn get_prop_mut(&mut self) -> &mut T {
        self
    }
}
impl<T, I, Head: Property<T, I>, Tail> Property<T, Inside<I>> for HCons<Head, Tail> {
    type Ref<'a>
        = Head::Ref<'a>
    where
        Self: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        self.head.get_prop()
    }
}
impl<T, I, Head: PropertyMut<T, I>, Tail> PropertyMut<T, Inside<I>> for HCons<Head, Tail> {
    type RefMut<'a>
        = Head::RefMut<'a>
    where
        Self: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.head.get_prop_mut()
    }
}
impl<T, I, Head, Tail: Property<T, I>> Property<T, Next<I>> for HCons<Head, Tail> {
    type Ref<'a>
        = Tail::Ref<'a>
    where
        Self: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        self.tail.get_prop()
    }
}
impl<T, I, Head, Tail: PropertyMut<T, I>> PropertyMut<T, Next<I>> for HCons<Head, Tail> {
    type RefMut<'a>
        = Tail::RefMut<'a>
    where
        Self: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.tail.get_prop_mut()
    }
}
impl<T, I, R: Deref> Property<T, Derefed<I>> for R
where
    R::Target: Property<T, I>,
{
    type Ref<'a>
        = <<R as Deref>::Target as Property<T, I>>::Ref<'a>
    where
        Self: 'a;

    fn get_prop(&self) -> Self::Ref<'_> {
        (**self).get_prop()
    }
}
impl<T, I, R: DerefMut> PropertyMut<T, Derefed<I>> for R
where
    R::Target: PropertyMut<T, I>,
{
    type RefMut<'a>
        = <<R as Deref>::Target as PropertyMut<T, I>>::RefMut<'a>
    where
        Self: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        (**self).get_prop_mut()
    }
}

/// Convenience trait to put the type in the function parameter rather than the trait
pub trait PropertyExt {
    type RefOf<'a, T, I>: Deref<Target = T> + 'a
    where
        Self: Property<T, I> + 'a;
    type RefMutOf<'a, T, I>: DerefMut<Target = T> + 'a
    where
        Self: PropertyMut<T, I> + 'a;

    fn get_prop_of<T, I>(&self) -> Self::Ref<'_>
    where
        Self: Property<T, I>,
    {
        self.get_prop()
    }
    fn get_prop_mut_of<T, I>(&mut self) -> Self::RefMut<'_>
    where
        Self: PropertyMut<T, I>,
    {
        self.get_prop_mut()
    }
}
impl<U> PropertyExt for U {
    type RefOf<'a, T, I>
        = <Self as Property<T, I>>::Ref<'a>
    where
        Self: Property<T, I> + 'a;
    type RefMutOf<'a, T, I>
        = <Self as PropertyMut<T, I>>::RefMut<'a>
    where
        Self: PropertyMut<T, I> + 'a;
}

/// Just check that everything works the way we want
fn _type_check() {
    struct Prop1;
    struct Prop2;
    struct Prop3;
    let _prop: &Prop1 = Prop1.get_prop();
    let prop_store = hlist![(), Prop1, 0i32, hlist![Prop2, "a string", &Prop3]];
    let _prop = prop_store.get_prop_of::<Prop1, _>();
    let _prop = prop_store.get_prop_of::<Prop2, _>();
    let _prop = prop_store.get_prop_of::<Prop3, _>();
}
