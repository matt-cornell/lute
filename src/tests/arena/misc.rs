use super::*;

#[test]
fn simple() {
    let mut arena = Arena::<u32>::new();
    let formic_acid = smiles!("CC(=O)O");
    let idx = arena.insert_mol(&formic_acid);
    let mol = arena.molecule(idx);
    assert_eq!(formic_acid.mass(), mol.mass());
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
