use petgraph::data::*;
use petgraph::visit::*;
use petgraph::Direction;
use std::fmt::{self, Debug, Formatter};

/// Like `petgraph`'s, but with `GetAdjacencyMatrix`.
#[derive(Clone)]
pub struct NodeFilter<G, F> {
    pub graph: G,
    pub filter: F,
}
impl<G, F> NodeFilter<G, F> {
    pub const fn new(graph: G, filter: F) -> Self {
        Self { graph, filter }
    }
}

impl<G: Debug, F> Debug for NodeFilter<G, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("NodeFilter")
            .field("graph", &self.graph)
            .finish_non_exhaustive()
    }
}

impl<G: GraphBase, F> GraphBase for NodeFilter<G, F> {
    type EdgeId = G::EdgeId;
    type NodeId = G::NodeId;
}

impl<G: GraphProp, F> GraphProp for NodeFilter<G, F> {
    type EdgeType = G::EdgeType;
    fn is_directed(&self) -> bool {
        self.graph.is_directed()
    }
}

impl<G: Data, F> Data for NodeFilter<G, F> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}

impl<G: DataMap, F: Fn(G::NodeId) -> bool> DataMap for NodeFilter<G, F> {
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        (self.filter)(id).then(|| self.graph.node_weight(id))?
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}

impl<G: NodeIndexable, F> NodeIndexable for NodeFilter<G, F> {
    fn from_index(&self, i: usize) -> Self::NodeId {
        self.graph.from_index(i)
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        self.graph.to_index(a)
    }
    fn node_bound(&self) -> usize {
        self.graph.node_bound()
    }
}

impl<G: GetAdjacencyMatrix, F> GetAdjacencyMatrix for NodeFilter<G, F> {
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        self.graph.is_adjacent(matrix, a, b)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoNodeIdentifiers for &'a NodeFilter<G, F>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
{
    type NodeIdentifiers =
        iter::NodeIdFilter<'a, <&'a G as IntoNodeIdentifiers>::NodeIdentifiers, &'a G, F>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::NodeIdFilter(self.graph.node_identifiers(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoNodeReferences for &'a NodeFilter<G, F>
where
    &'a G: IntoNodeReferences<NodeWeight = G::NodeWeight, NodeId = G::NodeId>,
{
    type NodeRef = <&'a G as IntoNodeReferences>::NodeRef;
    type NodeReferences =
        iter::NodeRefFilter<'a, <&'a G as IntoNodeReferences>::NodeReferences, &'a G, F>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeRefFilter(self.graph.node_references(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoEdgeReferences for &'a NodeFilter<G, F>
where
    &'a G: IntoEdgeReferences<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgeRef = <&'a G as IntoEdgeReferences>::EdgeRef;
    type EdgeReferences =
        iter::EdgeRefFilter<'a, <&'a G as IntoEdgeReferences>::EdgeReferences, &'a G, F>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeRefFilter(self.graph.edge_references(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoEdges for &'a NodeFilter<G, F>
where
    &'a G: IntoEdges<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type Edges = iter::EdgeRefFilter<'a, <&'a G as IntoEdges>::Edges, &'a G, F>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::EdgeRefFilter(self.graph.edges(a), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoNeighbors for &'a NodeFilter<G, F>
where
    &'a G: IntoNeighbors<EdgeId = G::EdgeId, NodeId = G::NodeId>,
{
    type Neighbors = iter::NodeIdFilter<'a, <&'a G as IntoNeighbors>::Neighbors, &'a G, F>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::NodeIdFilter(self.graph.neighbors(a), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoNeighborsDirected for &'a NodeFilter<G, F>
where
    &'a G: IntoNeighborsDirected<EdgeId = G::EdgeId, NodeId = G::NodeId>,
{
    type NeighborsDirected =
        iter::NodeIdFilter<'a, <&'a G as IntoNeighborsDirected>::NeighborsDirected, &'a G, F>;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        iter::NodeIdFilter(
            self.graph.neighbors_directed(a, dir),
            &self.graph,
            &self.filter,
        )
    }
}

impl<'a, G: Data, F: Fn(G::NodeId) -> bool> IntoEdgesDirected for &'a NodeFilter<G, F>
where
    &'a G: IntoEdgesDirected<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgesDirected =
        iter::EdgeRefFilter<'a, <&'a G as IntoEdgesDirected>::EdgesDirected, &'a G, F>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgeRefFilter(self.graph.edges_directed(a, dir), &self.graph, &self.filter)
    }
}

#[doc(hidden)]
pub mod iter {
    use super::*;

    pub struct NodeIdFilter<'a, I, G, F>(pub I, pub G, pub &'a F);
    impl<G: Data, I: Iterator<Item = G::NodeId>, F: Fn(G::NodeId) -> bool> Iterator
        for NodeIdFilter<'_, I, G, F>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0.by_ref().find(|&i| (self.2)(i))
        }
    }

    pub struct NodeRefFilter<'a, I, G, F>(pub I, pub G, pub &'a F);
    impl<G: IntoNodeReferences, I: Iterator<Item = G::NodeRef>, F: Fn(G::NodeId) -> bool> Iterator
        for NodeRefFilter<'_, I, G, F>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0.by_ref().find(|i| (self.2)(i.id()))
        }
    }

    pub struct EdgeRefFilter<'a, I, G, F>(pub I, pub G, pub &'a F);
    impl<G: IntoEdgeReferences, I: Iterator<Item = G::EdgeRef>, F: Fn(G::NodeId) -> bool> Iterator
        for EdgeRefFilter<'_, I, G, F>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0
                .by_ref()
                .find(|i| (self.2)(i.source()) && (self.2)(i.target()))
        }
    }
}
