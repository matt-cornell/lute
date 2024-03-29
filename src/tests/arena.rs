use crate::graph::algo::isomorphism::*;
use crate::prelude::*;
use petgraph::prelude::*;
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
fn double_insert() {
    let mut arena = Arena::<u32>::new();
    let methane = smiles!("C");
    let m1 = arena.insert_mol(&methane);
    let m2 = arena.insert_mol(&methane);
    assert_eq!(m1, m2);
}

#[test]
fn atomic_lookup() {
    let mut arena = Arena::<u32>::new();
    let ethanol = smiles!("CCO");
    let eth_idx = arena.insert_mol(&ethanol);
    let eth = arena.molecule(eth_idx);
    assert_eq!(eth.get_atom(0).map(|a| a.protons), Some(6));
    assert_eq!(eth.get_atom(1).map(|a| a.protons), Some(6));
    assert_eq!(eth.get_atom(2).map(|a| a.protons), Some(8));
    assert_eq!(eth.get_bond((0, 1)), Some(Bond::Single));
    assert_eq!(eth.get_bond((1, 2)), Some(Bond::Single));

    for i in 0..3 {
        eprintln!("atom {i}:");
        eprintln!(
            "  incoming: {:?}",
            eth.neighbors_directed(i.into(), Direction::Incoming)
                .map(|i| i.0)
                .collect::<Vec<_>>()
        );
        eprintln!(
            "  outgoing: {:?}",
            eth.neighbors_directed(i.into(), Direction::Outgoing)
                .map(|i| i.0)
                .collect::<Vec<_>>()
        );
        eprintln!(
            "  total: {:?}",
            eth.neighbors(i.into()).map(|i| i.0).collect::<Vec<_>>()
        );
    }

    assert!(is_isomorphic_matching(
        &ethanol,
        &ethanol,
        PartialEq::eq,
        PartialEq::eq
    ));
    assert!(is_isomorphic_matching(
        eth,
        eth,
        PartialEq::eq,
        PartialEq::eq
    ));
    assert!(is_isomorphic_matching(
        eth,
        &ethanol,
        PartialEq::eq,
        PartialEq::eq,
    ));
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

    let alcohol = arena.insert_mol(&smiles!("OR"));
    assert_eq!(alcohol, 0);
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");

    let isp_idx = arena.insert_mol(&smiles!("CC(C)O"));
    let isp = arena.molecule(isp_idx);

    // if these aren't equal, the output is just ugly
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");
    assert!(isp.contains(alcohol));
    eprintln!(
        "atoms: {:?}",
        isp.node_references()
            .map(|n| n.atom.protons)
            .collect::<Vec<_>>()
    );

    // this demonstrates how the ordering changes
    assert_eq!(isp.get_atom(0).map(|a| a.protons), Some(8));
    assert_eq!(isp.get_atom(1).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(2).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(3).map(|a| a.protons), Some(6));

    for (l, r) in itertools::iproduct!(0..4, 0..4) {
        if l < r {
            eprintln!("{l}--{r}: {:?}", isp.get_bond((l, r)));
        }
    }

    for i in 0..4 {
        eprintln!(
            "neighbors of {i}: {:?}",
            isp.neighbors(i.into()).collect::<Vec<_>>()
        );
    }

    // assert_eq!(isp.get_bond((0, 3)), Some(Bond::Single));
    // assert_eq!(isp.get_bond((1, 3)), Some(Bond::Single));
    // assert_eq!(isp.get_bond((2, 3)), Some(Bond::Single));
}
