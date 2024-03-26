use super::*;
use crate::arena::access::ArenaPtr;

gen_defs!(noinit Molecule as molecule);

#[no_mangle]
pub unsafe extern "C" fn CHEMSIM_new_molecule(arena: NonNull<Arena>, index: u32) -> Molecule {
    Molecule::from_mut_arena(&ArenaPtr::new(arena), index)
}

#[no_mangle]
pub unsafe extern "C" fn CHEMSIM_molecule_svg(
    mol: Molecule,
    alloc: unsafe extern "C" fn(len: usize, context: *mut ()) -> *mut c_char,
    context: *mut (),
) -> *mut c_char {
    let s = lib::fmt_as_svg(&crate::graph::GraphCompactor::<Molecule>::new(mol)).to_string();
    let len = s.len();
    let ptr = alloc(len, context);
    if !ptr.is_null() {
        std::slice::from_raw_parts_mut(ptr, len)
            .copy_from_slice(std::slice::from_raw_parts(s.as_ptr() as *const c_char, len));
    }
    ptr
}
