use super::*;

mod order;

#[test]
fn double_insert() {
    trace_capture!();
    let mut arena = Arena::<u32>::new();
    let benzene = smiles!("c1ccccc1");
    let m1 = arena.insert_mol(&benzene);
    let m2 = arena.insert_mol(&benzene);
    assert_eq!(m1, m2);
}

#[test]
fn isomorphism() {
    use crate::graph::algo::isomorphism::*;
    use petgraph::visit::*;
    use petgraph::{Incoming, Outgoing};
    let mut arena = Arena::<u32>::new();
    let acid = smiles!("OS(=O)(=O)O");
    let ix = arena.insert_mol(&acid);
    let ins_acid = arena.molecule(ix);
    println!("before");
    for node in acid.node_identifiers() {
        println!("{node:?}");
        for edge in acid.edges_directed(node, Incoming) {
            println!("=> {} -> {}", edge.source().index(), edge.target().index());
        }
        for edge in acid.edges_directed(node, Outgoing) {
            println!("<= {} -> {}", edge.source().index(), edge.target().index());
        }
    }
    println!("after");
    for node in ins_acid.node_identifiers() {
        println!("{node:?}");
        for edge in ins_acid.edges_directed(node, Incoming) {
            println!("=> {} -> {}", edge.source().0, edge.target().0);
        }
        for edge in ins_acid.edges_directed(node, Outgoing) {
            println!("<= {} -> {}", edge.source().0, edge.target().0);
        }
    }
    assert!(is_isomorphic_matching(
        &acid,
        ins_acid,
        PartialEq::eq,
        PartialEq::eq,
    ));
}

#[test]
fn alcohols() {
    trace_capture!();
    let mut arena = arena!(
        of u32:
        smiles!("O&"),     // alcohol
        smiles!("CCO"),    // ethanol
        smiles!("CC(C)O"), // isopropanol
    );
    tracing::info!("arena created");
    let orig = arena.parts.clone();

    let alcohol = arena.insert_mol(&smiles!("O&"));
    tracing::info!("alcohol inserted");
    assert_eq!(alcohol, 0);
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");

    let isp_idx = arena.insert_mol(&smiles!("CC(C)O"));
    tracing::info!("isopropanol inserted");
    let isp = arena.molecule(isp_idx);

    println!("{:#?}", arena.parts);

    // if these aren't equal, the output is just ugly
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");

    assert!(isp.contains(alcohol));
}
