use crate::core::*;
use petgraph::data::*;
use petgraph::visit::*;

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
