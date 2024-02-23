use petgraph::visit::*;
use super::*;

impl<Ix: Copy + PartialEq, R> GraphBase for Molecule<Ix, R> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
}

impl<Ix: Copy + PartialEq, R: Copy> GraphRef for Molecule<Ix, R> {}

impl<Ix: Copy + PartialEq, R> GraphProp for Molecule<Ix, R> {
    type EdgeType = petgraph::Undirected;
}
