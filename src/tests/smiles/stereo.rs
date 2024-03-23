use crate::prelude::*;

#[test]
fn ez() {
    use petgraph::graph::EdgeIndex;
    let but2ene = smiles!("CC=CC");
    assert_eq!(but2ene[EdgeIndex::new(1)], Bond::Double);
    let cis1 = smiles!("C/C=C/C");
    assert_eq!(cis1[EdgeIndex::new(1)], Bond::DoubleZ);
    let cis2 = smiles!("C\\C=C\\C");
    assert_eq!(cis2[EdgeIndex::new(1)], Bond::DoubleZ);
    let trans1 = smiles!("C/C=C\\C");
    assert_eq!(trans1[EdgeIndex::new(1)], Bond::DoubleE);
    let trans2 = smiles!("C\\C=C/C");
    assert_eq!(trans2[EdgeIndex::new(1)], Bond::DoubleE);
}

#[test]
fn rs() {
    use petgraph::graph::NodeIndex;
    let cfclbri = smiles!("FC(Cl)(Br)I");
    assert_eq!(cfclbri[NodeIndex::new(1)].data.chirality(), Chirality::None);
    let r = smiles!("F[C@](Cl)(Br)I");
    assert_eq!(r[NodeIndex::new(1)].data.chirality(), Chirality::R);
    let s = smiles!("F[C@@](Cl)(Br)I");
    assert_eq!(s[NodeIndex::new(1)].data.chirality(), Chirality::S);
}
