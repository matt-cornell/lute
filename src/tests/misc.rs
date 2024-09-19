//! Unit tests to affirm my sanity

use crate::core::Atom;

#[test]
fn atom_check() {
    // h2o
    let mut atom1 = Atom::new(8);
    atom1.add_hydrogens(2).unwrap();
    // roh
    let mut atom2 = Atom::new(8);
    atom2.add_hydrogens(1).unwrap();
    atom2.add_rs(1).unwrap();
    // -oh (in a molecule)
    let mut atom3 = Atom::new(8);
    atom3.add_hydrogens(1).unwrap();
    atom3.set_single_bonds(1).unwrap();

    // ror
    let mut atom4 = Atom::new(8);
    atom4.add_rs(2).unwrap();

    assert_ne!(atom1, atom2);
    assert!(!atom1.matches(&atom2));
    assert!(atom2.matches(&atom1));
    assert!(atom2.matches(&atom3));
    assert!(atom4.matches(&atom2));
}

/// I promise I'm not insane
#[test]
fn atom_eq_check() {
    let mut atom = Atom::new(6);
    atom.add_hydrogens(3).unwrap();
    atom.add_rs(1).unwrap();
    assert_eq!(atom, atom);
}
