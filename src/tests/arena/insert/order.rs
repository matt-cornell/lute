//! Tests for the arena insert order

use super::*;

#[test]
fn small_large() {
    let mut arena = Arena::<u8>::new();
    let ether = arena.insert_mol(&smiles!("O&&"));
    let dimethyl_ether = arena.insert_mol(&smiles!("COC"));
    assert!(arena.contains_group(dimethyl_ether, ether));
}

#[test]
fn large_small() {
    trace_capture!();
    let mut arena = Arena::<u8>::new();
    let dimethyl_ether = arena.insert_mol(&smiles!("COC"));
    let ether = arena.insert_mol(&smiles!("O&&"));
    println!("{:#?}", arena.expose_frags());
    assert!(arena.contains_group(dimethyl_ether, ether));
}

#[test]
fn small_med_large() {
    // trace_capture!();
    let mut arena = Arena::<u8>::new();
    let ether = arena.insert_mol(&smiles!("O&&"));
    let ketone = arena.insert_mol(&smiles!("C&&=O"));
    let ester = arena.insert_mol(&smiles!("C&(=O)O&"));
    let methyl_acetate = arena.insert_mol(&smiles!("CC(=O)OC"));

    assert!(arena.contains_group(ester, ether));
    assert!(arena.contains_group(ester, ketone));
    assert!(arena.contains_group(methyl_acetate, ether));
    assert!(arena.contains_group(methyl_acetate, ketone));
    assert!(arena.contains_group(methyl_acetate, ester));
}

#[test]
fn small_large_med() {
    // trace_capture!();
    let mut arena = Arena::<u8>::new();
    let ether = arena.insert_mol(&smiles!("O&&"));
    let methyl_acetate = arena.insert_mol(&smiles!("CC(=O)OC"));
    let ketone = arena.insert_mol(&smiles!("C&&=O"));
    let ester = arena.insert_mol(&smiles!("C&(=O)O&"));

    assert!(arena.contains_group(ester, ether));
    assert!(arena.contains_group(ester, ketone));
    assert!(arena.contains_group(methyl_acetate, ether));
    assert!(arena.contains_group(methyl_acetate, ketone));
    assert!(arena.contains_group(methyl_acetate, ester));
}

#[test]
fn large_med_small() {
    // trace_capture!();
    let mut arena = Arena::<u8>::new();
    let methyl_acetate = arena.insert_mol(&smiles!("CC(=O)OC"));
    let ester = arena.insert_mol(&smiles!("C&(=O)O&"));
    let ether = arena.insert_mol(&smiles!("O&&"));
    let ketone = arena.insert_mol(&smiles!("C&&=O"));

    assert!(arena.contains_group(ester, ether));
    assert!(arena.contains_group(ester, ketone));
    assert!(arena.contains_group(methyl_acetate, ether));
    assert!(arena.contains_group(methyl_acetate, ketone));
    assert!(arena.contains_group(methyl_acetate, ester));
}

#[test]
fn large_small_med() {
    // trace_capture!();
    let mut arena = Arena::<u8>::new();
    let methyl_acetate = arena.insert_mol(&smiles!("CC(=O)OC"));
    let ether = arena.insert_mol(&smiles!("O&&"));
    let ketone = arena.insert_mol(&smiles!("C&&=O"));
    let ester = arena.insert_mol(&smiles!("C&(=O)O&"));

    assert!(arena.contains_group(ester, ether));
    assert!(arena.contains_group(ester, ketone));
    assert!(arena.contains_group(methyl_acetate, ether));
    assert!(arena.contains_group(methyl_acetate, ketone));
    assert!(arena.contains_group(methyl_acetate, ester));
}

#[test]
fn large_small_med_small() {
    // trace_capture!();
    let mut arena = Arena::<u8>::new();
    let methyl_acetate = arena.insert_mol(&smiles!("CC(=O)OC"));
    let ether = arena.insert_mol(&smiles!("O&&"));
    let ester = arena.insert_mol(&smiles!("C&(=O)O&"));
    let ketone = arena.insert_mol(&smiles!("C&&=O"));

    assert!(arena.contains_group(ester, ether));
    assert!(arena.contains_group(ester, ketone));
    assert!(arena.contains_group(methyl_acetate, ether));
    assert!(arena.contains_group(methyl_acetate, ketone));
    assert!(arena.contains_group(methyl_acetate, ester));
}
