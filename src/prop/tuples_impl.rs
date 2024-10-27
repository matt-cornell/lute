use super::*;

const fn array_len<const N: usize>(_: [(); N]) -> usize {
    N
}

macro_rules! replace_idents {
    ($_id:ident $($sub:tt)*) => {
        $($sub)*
    };
}

macro_rules! impl_tuple {
    ($single:ident) => {
        impl_tuple!(@inner $single;);
    };
    ($first:ident, $($multi:ident),*) => {
        impl_tuple!(@inner $first $(,$multi)*;);
        impl_tuple!($($multi),*);
    };
    (@inner ; $($head:ident,)*) => {};
    (@inner $curr:ident $(,$tail:ident)*; $($head:ident,)*) => {
        impl<T, I, $($head,)* $($tail,)* $curr: Property<T, I>> Property<T, Nth<I, { array_len([$(replace_idents!($head ()),)*]) }>> for ($($head,)* $curr, $($tail,)*) {
            type Ref<'a> = $curr::Ref<'a> where Self: 'a, T: 'a;

            fn get_prop(&self) -> Self::Ref<'_> {
                let ($(replace_idents!($head _),)* curr, $(replace_idents!($tail _),)*) = self;
                curr.get_prop()
            }
        }
        impl<T, I, $($head,)* $($tail,)* $curr: PropertyMut<T, I>> PropertyMut<T, Nth<I, { array_len([$(replace_idents!($head ()),)*]) }>> for ($($head,)* $curr, $($tail,)*) {
            type RefMut<'a> = $curr::RefMut<'a> where Self: 'a, T: 'a;

            fn get_prop_mut(&mut self) -> Self::RefMut<'_> {
                let ($(replace_idents!($head _),)* curr, $(replace_idents!($tail _),)*) = self;
                curr.get_prop_mut()
            }
        }
        impl<T, I, $($head,)* $($tail,)* $curr: PropertyIMut<T, I>> PropertyIMut<T, Nth<I, { array_len([$(replace_idents!($head ()),)*]) }>> for ($($head,)* $curr, $($tail,)*) {
            type RefIMut<'a> = $curr::RefIMut<'a> where Self: 'a, T: 'a;

            fn get_prop_imut(&self) -> Self::RefIMut<'_> {
                let ($(replace_idents!($head _),)* curr, $(replace_idents!($tail _),)*) = self;
                curr.get_prop_imut()
            }
        }
        impl<T, I, $($head,)* $($tail,)* $curr: PropertyRef<T, I>> PropertyRef<T, Nth<I, { array_len([$(replace_idents!($head ()),)*]) }>> for ($($head,)* $curr, $($tail,)*) {
            fn extract_ref<'a>(r: Self::Ref<'a>) -> &'a T where Self: 'a {
                $curr::extract_ref(r)
            }
        }
        impl<T, I, $($head,)* $($tail,)* $curr: PropertyMutRef<T, I>> PropertyMutRef<T, Nth<I, { array_len([$(replace_idents!($head ()),)*]) }>> for ($($head,)* $curr, $($tail,)*) {
            fn extract_mut_ref<'a>(r: Self::RefMut<'a>) -> &'a mut T where Self: 'a {
                $curr::extract_mut_ref(r)
            }
        }
        impl_tuple!(@inner $($tail),*; $curr, $($head,)*);
    };
}
impl_tuple!(T0, T1, T2, T4, T5, T6, T7, T8, T9, T10, T11, T12);

#[allow(dead_code)]
fn type_check() {
    struct A;
    struct B;
    struct C;
    struct D;
    struct E;
    struct F;
    let prop_store = (A, B, C, D, (E, F));
    let _v = prop_store.get_prop_of::<A, _>();
    let _v = prop_store.get_prop_of::<B, _>();
    let _v = prop_store.get_prop_of::<C, _>();
    let _v = prop_store.get_prop_of::<D, _>();
    let _v = prop_store.get_prop_of::<E, _>();
    let _v = prop_store.get_prop_of::<F, _>();
}
