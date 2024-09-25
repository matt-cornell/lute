//! Create a compact indexable form of a graph.

use crate::prelude::DataValueMap;
use ahash::*;
use petgraph::data::*;
use petgraph::prelude::*;
use petgraph::visit::*;
use std::hash::Hash;

/// Override `NodeIndexable` so that they're always compact
#[derive(Debug, Clone)]
pub struct GraphCompactor<G: GraphBase> {
    pub graph: G,
    pub node_map: Vec<G::NodeId>,
    pub inv_map: HashMap<G::NodeId, usize>,
}
impl<N: Hash + Eq + Copy, G: GraphBase<NodeId = N>> GraphCompactor<G>
where
    for<'a> &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
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

impl<N: Hash + Eq + Copy, G: GraphBase<NodeId = N>> NodeIndexable for GraphCompactor<G> {
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

impl<'a, G: Data> IntoNodeIdentifiers for &'a GraphCompactor<G>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
{
    type NodeIdentifiers = <&'a G as IntoNodeIdentifiers>::NodeIdentifiers;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        self.graph.node_identifiers()
    }
}

impl<'a, G: Data> IntoNodeReferences for &'a GraphCompactor<G>
where
    &'a G: IntoNodeReferences<NodeId = G::NodeId, NodeWeight = G::NodeWeight>,
{
    type NodeRef = <&'a G as IntoNodeReferences>::NodeRef;
    type NodeReferences = <&'a G as IntoNodeReferences>::NodeReferences;

    fn node_references(self) -> Self::NodeReferences {
        self.graph.node_references()
    }
}

impl<'a, G: Data> IntoEdges for &'a GraphCompactor<G>
where
    &'a G: IntoEdges<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type Edges = <&'a G as IntoEdges>::Edges;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        self.graph.edges(a)
    }
}

impl<'a, G: Data> IntoEdgeReferences for &'a GraphCompactor<G>
where
    &'a G: IntoEdgeReferences<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgeRef = <&'a G as IntoEdgeReferences>::EdgeRef;
    type EdgeReferences = <&'a G as IntoEdgeReferences>::EdgeReferences;

    fn edge_references(self) -> Self::EdgeReferences {
        self.graph.edge_references()
    }
}

impl<'a, G: Data> IntoNeighbors for &'a GraphCompactor<G>
where
    &'a G: IntoNeighbors<EdgeId = G::EdgeId, NodeId = G::NodeId>,
{
    type Neighbors = <&'a G as IntoNeighbors>::Neighbors;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        self.graph.neighbors(a)
    }
}

impl<'a, G: Data> IntoNeighborsDirected for &'a GraphCompactor<G>
where
    &'a G: IntoNeighborsDirected<EdgeId = G::EdgeId, NodeId = G::NodeId>,
{
    type NeighborsDirected = <&'a G as IntoNeighborsDirected>::NeighborsDirected;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        self.graph.neighbors_directed(a, dir)
    }
}

impl<'a, G: Data> IntoEdgesDirected for &'a GraphCompactor<G>
where
    &'a G: IntoEdgesDirected<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgesDirected = <&'a G as IntoEdgesDirected>::EdgesDirected;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        self.graph.edges_directed(a, dir)
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

impl<G: DataValueMap> DataValueMap for GraphCompactor<G> {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        self.graph.node_weight(id)
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}

impl<G: Visitable> Visitable for GraphCompactor<G> {
    type Map = G::Map;

    fn reset_map(&self, map: &mut Self::Map) {
        self.graph.reset_map(map)
    }
    fn visit_map(&self) -> Self::Map {
        self.graph.visit_map()
    }
}
