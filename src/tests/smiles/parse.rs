use petgraph::visit::{IntoNodeReferences, NodeRef};

use crate::prelude::*;

#[test]
fn simple() {
    // methane
    smiles!("C");
    // methylene
    smiles!("[CH2]");
    // acetylene
    smiles!("C#C");
    // benzene
    smiles!("c1ccccc1");
    // carbon monoxoide
    smiles!("[C-]#[O+]");
}

#[test]
fn glucose() {
    use crate::graph::algo::isomorphism::is_isomorphic_matching;
    let canon = smiles!("C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O");
    let isomer =
        smiles!("C([C@@H]1[C@@H](C([C@](O1)(CO)O[C@@H]2C([C@H]([C@@H](C(O2)CO)O)O)O)O)O)O");
    assert!(is_isomorphic_matching(
        &canon,
        &isomer,
        &mut Atom::matches,
        &mut PartialEq::eq,
        false,
    ));
    let chirals = isomer
        .node_references()
        .filter_map(|r| {
            r.weight()
                .data
                .chirality()
                .is_chiral()
                .then_some(r.id().index())
        })
        .collect::<Vec<_>>();
    assert_eq!(chirals, &[1, 2, 4, 9, 11, 12]);
}
