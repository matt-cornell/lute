use crate::prelude::*;
use petgraph::visit::*;

#[test]
fn simple() {
    let mut arena = Arena::<u32>::new();
    let methane = smiles!("C");
    let idx = arena.insert_mol(&methane);
    let mol = arena.molecule(idx);
    assert_eq!(methane.mass(), mol.mass());
}

#[test]
fn alcohols() {
    let mut arena = arena!(
        of u32:
        smiles!("RO"),     // alcohol
        smiles!("CCO"),    // ethanol
        smiles!("CC(C)O"), // isopropanol
    );
    let orig = arena.parts.clone();
    eprintln!("arena: {:#?}", arena.expose_parts());

    let alcohol = arena.insert_mol(&smiles!("OR"));
    assert_eq!(alcohol, 0);
    eprintln!("arena (should be the same): {:#?}", arena.expose_parts());
    assert!(orig == arena.parts);

    let isp_idx = arena.insert_mol(&smiles!("CC(C)O"));
    eprintln!("arena (should be the same): {:#?}", arena.expose_parts());
    let isp = arena.molecule(isp_idx);

    // if these aren't equal, the output is just ugly
    assert!(orig == arena.parts);

    assert!(isp.contains(alcohol));
    eprintln!(
        "atoms: {:?}",
        isp.node_references()
            .map(|n| n.atom.protons)
            .collect::<Vec<_>>()
    );

    // this demonstrates how the ordering changes
    assert_eq!(isp.get_atom(0.into()).map(|a| a.protons), Some(8));
    assert_eq!(isp.get_atom(1.into()).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(2.into()).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(3.into()).map(|a| a.protons), Some(6));

    for (l, r) in itertools::iproduct!(0..4, 0..4) {
        if l < r {
            eprintln!("{l}--{r}: {:?}", isp.get_bond((l, r).into()));
        }
    }

    for i in 0..4 {
        eprintln!(
            "neighbors of {i}: {:?}",
            isp.neighbors(i.into()).collect::<Vec<_>>()
        );
    }

    // assert_eq!(isp.get_bond((0, 3).into()), Some(Bond::Single));
    // assert_eq!(isp.get_bond((1, 3).into()), Some(Bond::Single));
    // assert_eq!(isp.get_bond((2, 3).into()), Some(Bond::Single));
}
