use crate::core::*;
use crate::empirical::*;
use petgraph::data::*;
use petgraph::visit::*;

const ELECTRON_COUNTS: &[u8] = &[0, 2, 10, 18, 30, 36, 54, 86, 118];

fn unshared_electrons(protons: u8, charge: i8, bonded: u8) -> u8 {
    let elecs = (protons as i8 - charge) as u8;
    let shell = ELECTRON_COUNTS
        .binary_search(&elecs)
        .map_or_else(|e| e, |t| t);
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

    /// Find the hybridization of an atom. Panics if the atom doesn't exist, or if it has more than
    /// 6 bonds/lone pairs.
    fn hybridization(&self, id: Self::NodeId) -> u8
    where
        Self: DataMap + IntoEdges,
    {
        let atom = self.node_weight(id).unwrap();
        let bonds = atom.data.total_bonds();
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
        let bonded = bond_electrons.floor() as u8;
        let unshared = unshared_electrons(atom.protons, atom.charge, bonded);
        let groups = bonds + (unshared + 1) / 2;
        match groups {
            4 => {
                if neighbors_pi {
                    3
                } else {
                    4
                }
            }
            0..=6 => groups,
            g => panic!("atom has {g} bonding groups, which shouldn't be possible"),
        }
    }

    /// Number of pi electrons donated to an aromatic system.
    fn pi_donor(&self, id: Self::NodeId) -> u8
    where
        Self: DataMap + IntoEdges,
    {
        let atom = self.node_weight(id).unwrap();
        let bonds = atom.data.total_bonds();
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
        let bonded = bond_electrons.floor() as u8;
        let unshared = unshared_electrons(atom.protons, atom.charge, bonded);
        let mut groups = bonds + unshared.div_ceil(2);
        if neighbors_pi && groups == 4 {
            groups -= 1;
        }
        unshared.saturating_sub(groups * 2)
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
