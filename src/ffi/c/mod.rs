#![allow(clippy::missing_safety_doc)]
use super::*;
use std::ffi::c_char;
use std::mem::MaybeUninit;
use std::ptr::NonNull;

type Molecule<Ix = u32> = lib::Molecule<Ix, crate::arena::access::PtrAcc<Ix>>;

macro_rules! gen_defs {
    (alloc $ty:ty as $name:ident) => {
        concat_idents::concat_idents!(fn_name = LUTE_alloc_, $name {
            #[doc = concat!("Allocate a ", stringify!($ty), " without initialization.")]
            #[no_mangle]
            pub extern "C" fn fn_name() -> *mut MaybeUninit<$ty> {
                // On nightly, we can use `Box::try_new` for fallible allocation. On stable, we have to catch the panic.
                #[cfg(feature = "nightly")]
                let res = Box::try_new(MaybeUninit::uninit()).map(|b| Box::leak(b) as _);

                #[cfg(not(feature = "nightly"))]
                let res = std::panic::catch_unwind(|| Box::leak(Box::new(MaybeUninit::uninit())) as _);

                res.unwrap_or(std::ptr::null_mut())
            }
        });
    };
    (free $ty:ty as $name:ident) => {
        concat_idents::concat_idents!(fn_name = LUTE_free_, $name {
            #[doc = concat!("Free a ", stringify!($ty), ". Assumes it has already had `deinit` called.")]
            #[no_mangle]
            pub unsafe extern "C" fn fn_name(ptr: NonNull<MaybeUninit<$ty>>) {
                std::mem::drop(Box::from_raw(ptr.as_ptr()));
            }
        });
    };
    (init $ty:ty as $name:ident) => {
        concat_idents::concat_idents!(fn_name = LUTE_init_, $name {
            #[doc = concat!("Initialize an already initialized ", stringify!($ty), ". May leak memory if called on an already initialized one.")]
            #[no_mangle]
            pub unsafe extern "C" fn fn_name(mut ptr: NonNull<MaybeUninit<$ty>>) {
                *ptr.as_mut() = MaybeUninit::new($ty::default());
            }
        });
    };
    (deinit $ty:ty as $name:ident) => {
        concat_idents::concat_idents!(fn_name = LUTE_deinit_, $name {
            #[doc = concat!("Clean up, but don't free a ", stringify!($ty), ".")]
            pub unsafe extern "C" fn fn_name(mut ptr: NonNull<MaybeUninit<$ty>>) {
                std::mem::replace(ptr.as_mut(), MaybeUninit::uninit()).assume_init_drop();
            }
        });
    };
    (noinit $ty:ty as $name:ident) => {
        gen_defs!((alloc deinit free) $ty as $name);
    };
    (all $ty:ty as $name:ident) => {
        gen_defs!((alloc init deinit free) $ty as $name);
    };
    (($inst:ident) $ty:ty as $name:ident) => {
        gen_defs!($inst $ty as $name);
    };
    (($inst:ident $($insts:ident)*) $ty:ty as $name:ident) => {
        gen_defs!($inst $ty as $name);
        gen_defs!(($($insts)*) $ty as $name);
    };
}

pub mod arena;
pub mod graph;
pub mod molecule;
