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

    assert_ne!(atom1, atom2);
    assert!(atom2.matches(&atom1));
    assert!(!atom1.matches(&atom2));
}
