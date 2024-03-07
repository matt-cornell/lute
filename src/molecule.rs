use crate::core::*;
use crate::empirical::*;
use modular_bitfield::prelude::*;
use petgraph::data::*;
use petgraph::visit::*;

const ELECTRON_COUNTS: &[u8] = &[0, 2, 10, 18, 30, 36, 54, 86, 118];

fn unshared_electrons(protons: u8, charge: i8, bonded: u8) -> u8 {
    let elecs = (protons as i8 - charge) as u8;
    let shell = ELECTRON_COUNTS.binary_search(&elecs).unwrap_or_else(|e| e);
    if shell == 0 {
        return 0;
    }
    let lower = ELECTRON_COUNTS[shell];
    let upper = ELECTRON_COUNTS[shell + 1];
    // noble gas
    if elecs == lower {
        return 8 - bonded;
    }
    // p-block
    if let Some(res) = elecs.checked_sub(upper - 6) {
        return res + 2 - bonded;
    }
    // s-block
    if elecs <= lower + 2 {
        return elecs - lower - bonded;
    }
    // some metal
    2_u8.saturating_sub(bonded)
}

#[bitfield]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomInfoFields {
    pub bonded: B7,
    pub neighbors_pi: bool,
}

/// `AtomInfo` holds some more data about the atoms this atom is bonded to. It can calculate things
/// like unshared electrons and hybridization because of this.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomInfo {
    pub atom: Atom,
    pub fields: AtomInfoFields,
}

impl AtomInfo {
    /// The number of unshared (non-bonding) valence electrons.
    pub fn unshared(&self) -> u8 {
        unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded())
    }
    /// The hybridization of this atom. 1 is s, 2 is sp, 3 is sp2, etc.
    pub fn hybridization(&self) -> u8 {
        let unshared =
            unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded());
        let groups = self.atom.data.total_bonds() + (unshared + 1) / 2;
        match groups {
            4 => {
                if self.fields.neighbors_pi() {
                    3
                } else {
                    4
                }
            }
            0..=6 => groups,
            g => panic!("atom has {g} bonding groups, which shouldn't be possible"),
        }
    }
    /// The number of delocalized pi electrons available for an aromatic system.
    pub fn pi_donor(&self) -> u8 {
        let unshared =
            unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded());
        let mut groups = self.atom.data.total_bonds() + unshared.div_ceil(2);
        if self.fields.neighbors_pi() && groups == 4 {
            groups -= 1;
        }
        unshared.saturating_sub(groups * 2)
    }
}

/// It's a lot easier to return an atom by value, since we have to handle borrows and locking, etc.
/// We rely on this instead of `DataMap` so that `Molecule` can implement the full set of functions.
pub trait ValueMolecule: Data<NodeWeight = Atom, EdgeWeight = Bond> {
    fn get_atom(&self, id: Self::NodeId) -> Option<Atom>;
    fn get_bond(&self, id: Self::EdgeId) -> Option<Bond>;
}
impl<T: DataMap<NodeWeight = Atom, EdgeWeight = Bond>> ValueMolecule for T {
    fn get_atom(&self, id: Self::NodeId) -> Option<Atom> {
        self.node_weight(id).copied()
    }
    fn get_bond(&self, id: Self::EdgeId) -> Option<Bond> {
        self.edge_weight(id).copied()
    }
}

/// Most of the molecule operations. This should be implemented for whatever is necessary, and
/// impossible to implement for additional types.
pub trait Molecule:
    GraphProp<EdgeType = petgraph::Undirected> + Data<NodeWeight = Atom, EdgeWeight = Bond>
{
    /// Find the mass of a molecule.
    fn mass(&self) -> f32
    where
        Self: IntoNodeReferences,
    {
        self.node_references().map(|a| a.weight().mass()).sum()
    }

    fn empirical(&self) -> EmpiricalFormula
    where
        Self: IntoNodeReferences,
    {
        self.node_references().map(|a| *a.weight()).collect()
    }

    /// Get the `AtomData` for a molecule in the graph.
    fn atom_data(&self, id: Self::NodeId) -> AtomInfo
    where
        Self: ValueMolecule + IntoEdges,
    {
        let atom = self.get_atom(id).unwrap();
        let mut neighbors_pi = false;
        let mut bond_electrons = 0.0;
        for edge in self.edges(id) {
            bond_electrons += edge.weight().bond_count();
            let neighbor = if edge.source() == id {
                if edge.target() == id {
                    panic!("Molecule should not contain self-loops!")
                } else {
                    edge.target()
                }
            } else {
                edge.source()
            };
            if !neighbors_pi {
                neighbors_pi = self.edges(neighbor).any(|e| e.weight().bond_count() > 1.0);
            }
        }
        AtomInfo {
            atom,
            fields: AtomInfoFields::new()
                .with_bonded(bond_electrons.floor() as _)
                .with_neighbors_pi(neighbors_pi),
        }
    }

    /// Resolve tautomers. This also ensures that aromaticicty is tracked, and rings have the
    /// "correct" conjugation.
    fn resolve_tautomer(&mut self)
    where
        Self: IntoEdges + DataMapMut,
    {
        todo!()
    }
}
impl<
        T: GraphProp<EdgeType = petgraph::Undirected> + Data<NodeWeight = Atom, EdgeWeight = Bond>,
    > Molecule for T
{
}
