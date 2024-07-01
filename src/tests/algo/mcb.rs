use crate::graph::algo::mcb::*;
use crate::smiles;

/// Rings could be rotated or reversed, this checks all of them, panicking if none match.
#[track_caller]
fn assert_ring_eq<T: PartialEq + std::fmt::Debug>(a: &[T], b: &mut [T]) {
    if a.len() == b.len() {
        for _ in 0..a.len() {
            b.rotate_right(1);
            if a == b {
                return;
            }
        }
        b.reverse();
        for _ in 1..a.len() {
            b.rotate_right(1);
            if a == b {
                return;
            }
        }
    }
    assert_eq!(a, b);
}

#[test]
fn acylcic() {
    let methylacetate = smiles!("COC(=O)C");
    assert_eq!(num_cycles(&methylacetate), 0);
    let mut it = CycleBasis::new_struct(&methylacetate);
    assert_eq!(it.next(), None);
}

#[test]
fn monocycles() {
    {
        let benzene = smiles!("c1ccccc1");
        assert_eq!(num_cycles(&benzene), 1);
        let mut benzene_it = CycleBasis::new_struct(&benzene);
        assert_ring_eq(
            &benzene_it.next().expect("benzene is cyclic").1,
            &mut [0, 1, 2, 3, 4, 5],
        );
    }
    {
        let toluene = smiles!("c1ccccc1C");
        assert_eq!(num_cycles(&toluene), 1);
        let mut toluene_it = CycleBasis::new_struct(&toluene);
        assert_ring_eq(
            &toluene_it.next().expect("toluene is cyclic").1,
            &mut [0, 1, 2, 3, 4, 5],
        );
    }
    {
        let bisphenol_a = smiles!("Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C");
        assert_eq!(num_cycles(&bisphenol_a), 2);
        let bpa_it = CycleBasis::new_struct(&bisphenol_a);
        let mut cycles: Vec<_> = bpa_it.collect();
        cycles.sort();
        assert_eq!(cycles.len(), 2);
        assert_ring_eq(&cycles[0].1, &mut [1, 2, 3, 4, 5, 6]);
        assert_ring_eq(&cycles[1].1, &mut [8, 9, 10, 11, 13, 14]);
    }
}

#[test]
fn fused_cycles() {
    {
        let napthalene = smiles!("c1c2ccccc2ccc1");
        assert_eq!(num_cycles(&napthalene), 2);
        let naptha_it = CycleBasis::new_struct(&napthalene);
        let mut cycles: Vec<_> = naptha_it.collect();
        cycles.sort();
        assert_eq!(cycles.len(), 2);
        // not gonna try to figure this out right now
    }
    {
        let propellane = smiles!("C1(C2)(C3)C23C1");
        assert_eq!(num_cycles(&propellane), 3);
        let prop_it = CycleBasis::new_struct(&propellane);
        assert_eq!(prop_it.count(), 3); // probably never gonna try to get this working
    }
    {
        let cubane = smiles!("C12C3C4C1C5C2C3C45");
        assert_eq!(num_cycles(&cubane), 5);
        let cube_it = CycleBasis::new_struct(&cubane);
        assert_eq!(cube_it.count(), 5); // probably never gonna try to get this working
    }
}
