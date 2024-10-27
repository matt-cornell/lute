//! Type-based property access.
//!
//! Uses the same basic idea as [`frunk`] does in its `HList` methods, but adds in recursive searches.
//! Looking at the implementors of these traits is a mess, so to summarize:
//! - Every type can access itself
//! - `HList`s can access any of their elements
//! - Tuples up to 12 elements can access any of their elements
//! - `&mut T` and `Box` forward access to their contained type
//! - `&T`, `Rc`, and `Arc` all provide immutable and interior-mutable access
//! - `RefCell`, `std::sync::Mutex`, and `lock_api::Mutex` all give interior-mutable access
//!
//! ## Examples
//! ```
//! use frunk::hlist;
//! use lute::prelude::*;
//! use petgraph::prelude::UnGraph;
//!
//! struct Prop1;
//! struct Prop2;
//! struct Prop3;
//! let _prop = Prop1.get_prop_of::<Prop1, _>(); // reflexive access
//! let prop_store = hlist![(), Prop1, 0i32, hlist![Prop2, "a string", &Prop3]];
//! let _prop = prop_store.get_prop_of::<Prop1, _>(); // access to a value in a `HList`
//! let _prop = prop_store.get_prop_of::<Prop2, _>();
//! let _prop = prop_store.get_prop_of::<Prop3, _>(); // works recursively too
//!
//! let graph = smiles!("CCC");
//! let mut prop = String::new();
//! needs_graph::<&UnGraph<_, _>, _, _, _>(hlist![&graph, &mut prop]); // pass in an `HList` with the properties
//! needs_graph::<&UnGraph<_, _>, _, _, _>((&graph, &mut prop)); // tuples work too!
//!
//! let mut arena = Arena::<u32, std::cell::RefCell<String>>::new();
//! let mol = arena.insert_mol(&graph).in_arena(&arena);
//! needs_graph::<Molecule<_, _>, _, _, _>(mol); // and a `Molecule` can pass any properties that associated data has
//!
//! fn needs_graph<
//!     G: petgraph::visit::IntoNodeIdentifiers<NodeId: std::fmt::Debug>,
//!     P: prop::Property<G, IG> + prop::PropertyMut<String, IS>,
//!     IG, // this is automatically determined, as long as it's not ambiguous
//!     IS, // same as above
//! >(
//!     mut args: P,
//! ) {
//!     use itertools::Itertools;
//!     let nodes = args
//!         .get_prop_of::<G, IG>()
//!         .node_identifiers()
//!         .map(|i| format!("{i:?}"))
//!         .join(",");
//!     *args.get_prop_mut_of::<String, IS>() = nodes;
//! }
//! ```

use crate::arena::{ArenaAccessorRef, Molecule};
use frunk::hlist;
use frunk::HCons;
use indices::*;
use petgraph::graph::IndexType;
use std::ops::{Deref, DerefMut};

mod lock_api_impl;
mod refs_impl;
mod std_imut_impl;
mod tuples_impl;

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
    /// The value at the nth position has the property we want
    pub struct Nth<T, const N: usize> {
        _marker: PhantomData<[T; N]>,
    }
}

/// Immutable property access
///
/// Takes `&self` and returns a type similar to `&T`.
pub trait Property<T, I> {
    /// Reference type. We can't just return a reference because of interior mutability.
    type Ref<'a>: Deref<Target = T> + 'a
    where
        Self: 'a,
        T: 'a;

    fn get_prop(&self) -> Self::Ref<'_>;
}

/// Mutable property access
///
/// Takes `&mut self` and returns a type similar to `&mut T`.
pub trait PropertyMut<T, I>: Property<T, I> {
    /// Mutable reference type. Similarly to [`Property::Ref`], we can't return a reference because of the possibility of there being guards.
    /// Even though a `RefCell` or an `RwLock` can return a mutable reference, it could be behind a `Rc` or `Arc`.
    type RefMut<'a>: DerefMut<Target = T> + 'a
    where
        Self: 'a,
        T: 'a;

    fn get_prop_mut(&mut self) -> Self::RefMut<'_>;
}

/// Interior-mutable property access
///
/// Takes `&self` and returns a type similar to `&mut T`.
pub trait PropertyIMut<T, I>: PropertyMut<T, I> {
    /// This type right here, officer! This is what's been driving me insane!
    /// This could in theory be the same as `RefMut`, which would mean one less associated type, but then we couldn't take advantage of `get_mut` methods.
    type RefIMut<'a>: DerefMut<Target = T> + 'a
    where
        Self: 'a,
        T: 'a;

    fn get_prop_imut(&self) -> Self::RefIMut<'_>;
}

/// Trait for properties that just return a reference with no additional drop guards.
///
/// I'd love to have implemented this in terms of bounds on the GAT but Rust isn't there yet.
pub trait PropertyRef<T, I>: Property<T, I> {
    /// This should be a no-op.
    fn extract_ref(r: Self::Ref<'_>) -> &T;
    #[inline]
    fn get_prop_ref(&self) -> &T {
        Self::extract_ref(self.get_prop())
    }
}
/// Same as [`PropertyRef`], but for the mutable reference type.
///
/// There's no interior mutable counterpart because there's no safe way to convert `&T` to `&mut T` without a drop guard.
pub trait PropertyMutRef<T, I>: PropertyMut<T, I> {
    /// This should be a no-op.
    fn extract_mut_ref(r: Self::RefMut<'_>) -> &mut T;
    #[inline]
    fn get_prop_mut_ref(&mut self) -> &mut T {
        Self::extract_mut_ref(self.get_prop_mut())
    }
}

impl<T> Property<T, Current> for T {
    type Ref<'a>
        = &'a T
    where
        Self: 'a;

    #[inline(always)]
    fn get_prop(&self) -> &T {
        self
    }
}
impl<T> PropertyMut<T, Current> for T {
    type RefMut<'a>
        = &'a mut T
    where
        Self: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> &mut T {
        self
    }
}
impl<T> PropertyRef<T, Current> for T {
    #[inline(always)]
    fn extract_ref(r: Self::Ref<'_>) -> &T {
        r
    }
}
impl<T> PropertyMutRef<T, Current> for T {
    #[inline(always)]
    fn extract_mut_ref(r: Self::RefMut<'_>) -> &mut T {
        r
    }
}
impl<T, I, Head: Property<T, I>, Tail> Property<T, Inside<I>> for HCons<Head, Tail> {
    type Ref<'a>
        = Head::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        self.head.get_prop()
    }
}
impl<T, I, Head: PropertyMut<T, I>, Tail> PropertyMut<T, Inside<I>> for HCons<Head, Tail> {
    type RefMut<'a>
        = Head::RefMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.head.get_prop_mut()
    }
}
impl<T, I, Head: PropertyIMut<T, I>, Tail> PropertyIMut<T, Inside<I>> for HCons<Head, Tail> {
    type RefIMut<'a>
        = Head::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        self.head.get_prop_imut()
    }
}
impl<T, I, Head: PropertyRef<T, I>, Tail> PropertyRef<T, Inside<I>> for HCons<Head, Tail> {
    #[inline(always)]
    fn extract_ref<'a>(r: Self::Ref<'a>) -> &'a T
    where
        Self: 'a,
    {
        Head::extract_ref(r)
    }
}
impl<T, I, Head: PropertyMutRef<T, I>, Tail> PropertyMutRef<T, Inside<I>> for HCons<Head, Tail> {
    #[inline(always)]
    fn extract_mut_ref<'a>(r: Self::RefMut<'a>) -> &'a mut T
    where
        Self: 'a,
    {
        Head::extract_mut_ref(r)
    }
}
impl<T, I, Head, Tail: Property<T, I>> Property<T, Next<I>> for HCons<Head, Tail> {
    type Ref<'a>
        = Tail::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        self.tail.get_prop()
    }
}
impl<T, I, Head, Tail: PropertyMut<T, I>> PropertyMut<T, Next<I>> for HCons<Head, Tail> {
    type RefMut<'a>
        = Tail::RefMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        self.tail.get_prop_mut()
    }
}
impl<T, I, Head, Tail: PropertyIMut<T, I>> PropertyIMut<T, Next<I>> for HCons<Head, Tail> {
    type RefIMut<'a>
        = Tail::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        self.tail.get_prop_imut()
    }
}
impl<T, I, Head, Tail: PropertyRef<T, I>> PropertyRef<T, Next<I>> for HCons<Head, Tail> {
    #[inline(always)]
    fn extract_ref<'a>(r: Self::Ref<'a>) -> &'a T
    where
        Self: 'a,
    {
        Tail::extract_ref(r)
    }
}
impl<T, I, Head, Tail: PropertyMutRef<T, I>> PropertyMutRef<T, Next<I>> for HCons<Head, Tail> {
    #[inline(always)]
    fn extract_mut_ref<'a>(r: Self::RefMut<'a>) -> &'a mut T
    where
        Self: 'a,
    {
        Tail::extract_mut_ref(r)
    }
}

impl<T, I, Ix: IndexType, R: ArenaAccessorRef<Ix = Ix, Data: Property<T, I>>> Property<T, Inside<I>>
    for Molecule<Ix, R>
{
    type Ref<'a>
        = <R::Data as Property<T, I>>::Ref<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop(&self) -> Self::Ref<'_> {
        R::extract_mapped_ref(self.data()).get_prop()
    }
}
impl<T, I, Ix: IndexType, R: ArenaAccessorRef<Ix = Ix, Data: PropertyIMut<T, I>>>
    PropertyMut<T, Inside<I>> for Molecule<Ix, R>
{
    type RefMut<'a>
        = <R::Data as PropertyIMut<T, I>>::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
        R::extract_mapped_ref(self.data()).get_prop_imut()
    }
}
impl<T, I, Ix: IndexType, R: ArenaAccessorRef<Ix = Ix, Data: PropertyIMut<T, I>>>
    PropertyIMut<T, Inside<I>> for Molecule<Ix, R>
{
    type RefIMut<'a>
        = <R::Data as PropertyIMut<T, I>>::RefIMut<'a>
    where
        Self: 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_imut(&self) -> Self::RefIMut<'_> {
        R::extract_mapped_ref(self.data()).get_prop_imut()
    }
}

/// Convenience trait to put the type in the function parameter rather than the trait
pub trait PropertyExt {
    type RefOf<'a, T, I>: Deref<Target = T> + 'a
    where
        Self: Property<T, I> + 'a,
        T: 'a;
    type RefMutOf<'a, T, I>: DerefMut<Target = T> + 'a
    where
        Self: PropertyMut<T, I> + 'a,
        T: 'a;
    type RefIMutOf<'a, T, I>: DerefMut<Target = T> + 'a
    where
        Self: PropertyIMut<T, I> + 'a,
        T: 'a;

    #[inline(always)]
    fn get_prop_of<T, I>(&self) -> Self::Ref<'_>
    where
        Self: Property<T, I>,
    {
        self.get_prop()
    }
    #[inline(always)]
    fn get_prop_mut_of<T, I>(&mut self) -> Self::RefMut<'_>
    where
        Self: PropertyMut<T, I>,
    {
        self.get_prop_mut()
    }
    #[inline(always)]
    fn get_prop_imut_of<T, I>(&self) -> Self::RefIMut<'_>
    where
        Self: PropertyIMut<T, I>,
    {
        self.get_prop_imut()
    }
    #[inline(always)]
    fn get_prop_ref_of<T, I>(&self) -> &T
    where
        Self: PropertyRef<T, I>,
    {
        self.get_prop_ref()
    }
    #[inline(always)]
    fn get_prop_mut_ref_of<T, I>(&mut self) -> &mut T
    where
        Self: PropertyMutRef<T, I>,
    {
        self.get_prop_mut_ref()
    }
}
impl<U> PropertyExt for U {
    type RefOf<'a, T, I>
        = <Self as Property<T, I>>::Ref<'a>
    where
        Self: Property<T, I> + 'a,
        T: 'a;
    type RefMutOf<'a, T, I>
        = <Self as PropertyMut<T, I>>::RefMut<'a>
    where
        Self: PropertyMut<T, I> + 'a,
        T: 'a;
    type RefIMutOf<'a, T, I>
        = <Self as PropertyIMut<T, I>>::RefIMut<'a>
    where
        Self: PropertyIMut<T, I> + 'a,
        T: 'a;
}

/// Just check that everything works the way we want
#[allow(dead_code)]
fn type_check() {
    struct Prop1;
    struct Prop2;
    struct Prop3;
    let _prop: &Prop1 = Prop1.get_prop();
    let prop_store = hlist![(), Prop1, 0i32, hlist![Prop2, "a string", &Prop3]];
    let _prop = prop_store.get_prop_of::<Prop1, _>();
    let _prop = prop_store.get_prop_of::<Prop2, _>();
    let _prop = prop_store.get_prop_of::<Prop3, _>();

    let graph = crate::smiles!("CCC");
    let mut prop = String::new();
    needs_graph::<&petgraph::prelude::UnGraph<_, _>, _, _, _>(hlist![&graph, &mut prop]);

    let mut arena = crate::arena::Arena::<u32, std::cell::RefCell<String>>::new();
    let mol = arena.insert_mol(&graph).in_arena(&arena);
    needs_graph::<Molecule<_, _>, _, _, _>(mol);
}

#[allow(dead_code)]
fn needs_graph<
    G: petgraph::visit::IntoNodeIdentifiers<NodeId: std::fmt::Debug>,
    P: Property<G, IG> + PropertyMut<String, IS>,
    IG,
    IS,
>(
    mut args: P,
) {
    use itertools::Itertools;
    let nodes = args
        .get_prop_of::<G, IG>()
        .node_identifiers()
        .map(|i| format!("{i:?}"))
        .join(",");
    *args.get_prop_mut_of::<String, IS>() = nodes;
}
