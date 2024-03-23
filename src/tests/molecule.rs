use crate::prelude::*;

#[test]
fn mass() {
    let methane = smiles!("C");
    assert_eq!(methane.mass().round(), 16.0);
}
