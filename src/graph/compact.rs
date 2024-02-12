//! Create a compact indexable form of a graph.

use petgraph::data::*;
use petgraph::visit::*;
use std::collections::HashMap;
use std::hash::Hash;

/// Override `NodeIndexable` so that they're always compact
#[derive(Debug, Clone)]
pub struct GraphCompactor<G: GraphBase> {
    graph: G,
    node_map: Vec<G::NodeId>,
    inv_map: HashMap<G::NodeId, usize>,
}
impl<N: Hash + Eq + Copy, G: GraphBase<NodeId = N> + GraphRef + IntoNodeIdentifiers>
    GraphCompactor<G>
{
    pub fn new(graph: G) -> Self {
        let node_map: Vec<G::NodeId> = graph.node_identifiers().collect();
        let inv_map = node_map.iter().enumerate().map(|(v, &k)| (k, v)).collect();
        Self {
            graph,
            node_map,
            inv_map,
        }
    }
}

impl<G: GraphBase> GraphBase for GraphCompactor<G> {
    type NodeId = G::NodeId;
    type EdgeId = G::EdgeId;
}

impl<G: GraphProp> GraphProp for GraphCompactor<G> {
    type EdgeType = G::EdgeType;
    fn is_directed(&self) -> bool {
        self.graph.is_directed()
    }
}

impl<G: GraphBase> NodeCount for GraphCompactor<G> {
    fn node_count(&self) -> usize {
        self.node_map.len()
    }
}

impl<N: Hash + Eq + Copy, G: GraphBase<NodeId = N> + NodeIndexable> NodeIndexable
    for GraphCompactor<G>
{
    fn from_index(&self, i: usize) -> Self::NodeId {
        self.node_map[i]
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        self.inv_map[&a]
    }
    fn node_bound(&self) -> usize {
        self.node_map.len()
    }
}

impl<N: Hash + Eq + Copy, G: GraphBase<NodeId = N> + NodeIndexable> NodeCompactIndexable
    for GraphCompactor<G>
{
}

impl<G: EdgeCount> EdgeCount for GraphCompactor<G> {
    fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
}

impl<G: GetAdjacencyMatrix> GetAdjacencyMatrix for GraphCompactor<G> {
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        self.graph.is_adjacent(matrix, a, b)
    }
}

impl<G: Data> Data for GraphCompactor<G> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}

impl<G: DataMap> DataMap for GraphCompactor<G> {
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        self.graph.node_weight(id)
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}
