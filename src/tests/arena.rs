use crate::graph::algo::isomorphism::*;
use crate::prelude::*;
use petgraph::prelude::*;
use petgraph::visit::*;

macro_rules! trace_capture {
    () => {
        use tracing_subscriber::filter::{LevelFilter, Targets};
        use tracing_subscriber::fmt::*;
        use tracing_subscriber::prelude::*;

        let targets = Targets::new()
            .with_target("chem_sim::parse", LevelFilter::INFO)
            .with_target("chem_sim::arena::molecule", LevelFilter::DEBUG)
            .with_target("chem_sim::arena::arena", LevelFilter::TRACE);

        let formatter = layer().compact().with_writer(TestWriter::new());

        let _guard = tracing_subscriber::registry()
            .with(targets)
            .with(formatter)
            .set_default();
    };
}

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
    trace_capture!();
    let mut arena = Arena::<u32>::new();
    let ethanol = smiles!("CCO");
    let eth_idx = arena.insert_mol(&ethanol);
    let eth = arena.molecule(eth_idx);
    assert_eq!(eth.get_atom(0).map(|a| a.protons), Some(6));
    assert_eq!(eth.get_atom(1).map(|a| a.protons), Some(6));
    assert_eq!(eth.get_atom(2).map(|a| a.protons), Some(8));
    assert_eq!(eth.get_bond((0, 1)), Some(Bond::Single));
    assert_eq!(eth.get_bond((1, 2)), Some(Bond::Single));

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
    trace_capture!();
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

    // this demonstrates how the ordering changes
    assert_eq!(isp.get_atom(0).map(|a| a.protons), Some(8));
    assert_eq!(isp.get_atom(1).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(2).map(|a| a.protons), Some(6));
    assert_eq!(isp.get_atom(3).map(|a| a.protons), Some(6));


    // assert_eq!(isp.get_bond((0, 3)), Some(Bond::Single));
    // assert_eq!(isp.get_bond((1, 3)), Some(Bond::Single));
    // assert_eq!(isp.get_bond((2, 3)), Some(Bond::Single));
}
