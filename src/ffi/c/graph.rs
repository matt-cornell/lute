use super::*;

gen_defs!(all Graph as graph);

#[no_mangle]
pub unsafe extern "C" fn LUTE_graph_svg(
    mol: NonNull<Graph>,
    alloc: unsafe extern "C" fn(len: usize, context: *mut ()) -> *mut c_char,
    context: *mut (),
) -> *mut c_char {
    let s = lib::fmt_as_svg(mol.as_ref()).to_string();
    let len = s.len();
    let ptr = alloc(len, context);
    if !ptr.is_null() {
        std::slice::from_raw_parts_mut(ptr, len)
            .copy_from_slice(std::slice::from_raw_parts(s.as_ptr() as *const c_char, len));
    }
    ptr
}
