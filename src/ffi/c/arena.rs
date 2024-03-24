use super::*;

gen_defs!(all Arena as arena);

#[no_mangle]
pub unsafe extern "C" fn CHEMSIM_add_moleclue(
    mut arena: NonNull<Arena>,
    graph: NonNull<Graph>,
) -> u32 {
    arena.as_mut().insert_mol(graph.as_ref())
}
