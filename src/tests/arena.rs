use crate::prelude::*;

macro_rules! trace_capture {
    () => {
        use tracing_subscriber::filter::{LevelFilter, Targets};
        use tracing_subscriber::prelude::*;

        let targets = Targets::new()
            .with_target("lute::parse", LevelFilter::INFO)
            .with_target("lute::arena::molecule", LevelFilter::DEBUG)
            .with_target("lute::arena::arena", LevelFilter::TRACE);

        let formatter = tracing_subscriber::fmt::layer().with_test_writer();

        let _guard = tracing_subscriber::registry()
            .with(targets)
            .with(formatter)
            .set_default();
    };
}

#[test]
fn simple() {
    let mut arena = Arena::<u32>::new();
    let methane = smiles!("CC(=O)O");
    let idx = arena.insert_mol(&methane);
    let mol = arena.molecule(idx);
    assert_eq!(methane.mass(), mol.mass());
}

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
    trace_capture!();
    let mut arena = Arena::<u32>::new();
    let acid = smiles!("OS(=O)(=O)O");
    let ix = arena.insert_mol(&acid);
    let ins_acid = arena.molecule(ix);
    assert!(is_isomorphic_matching(
        &acid,
        ins_acid,
        PartialEq::eq,
        PartialEq::eq
    ));
}

#[test]
fn atomic_lookup() {
    trace_capture!();
    let mut arena = Arena::<u32>::new();
    for (mol, atoms) in [
        ("CCO", &[6, 6, 8][..]),
        ("CN", &[6, 7]),
        ("O=C=O", &[8, 6, 8]),
    ] {
        let mol_idx = arena.insert_mol(&SmilesParser::new(mol).parse().unwrap());
        let ins = arena.molecule(mol_idx);
        for (n, &protons) in atoms.iter().enumerate() {
            assert_eq!(
                ins.get_atom(n as u32).map(|a| a.protons),
                Some(protons),
                "assertion failed, mol={mol:?}, n={n}, p={protons}"
            );
        }
    }
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
    let orig = arena.parts.clone();

    let alcohol = arena.insert_mol(&smiles!("O&"));
    assert_eq!(alcohol, 0);
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");

    let isp_idx = arena.insert_mol(&smiles!("CC(C)O"));
    let isp = arena.molecule(isp_idx);

    // if these aren't equal, the output is just ugly
    assert_eq!(orig.len(), arena.parts.len(), "arena length changed");
    assert!(orig == arena.parts, "arena changed");
    
    assert!(isp.contains(alcohol));
}
