use super::misc::DataValueMap;
use petgraph::data::*;
use petgraph::visit::*;
use petgraph::Direction;

#[derive(Debug, Clone, Copy)]
pub struct RangeFiltered<G> {
    pub graph: G,
    pub start: usize,
    pub end: usize,
}
impl<G> RangeFiltered<G> {
    pub const fn new(graph: G, start: usize, end: usize) -> Self {
        Self { graph, start, end }
    }
}

impl<G: GraphBase> GraphBase for RangeFiltered<G> {
    type EdgeId = G::EdgeId;
    type NodeId = G::NodeId;
}

impl<G: GraphProp> GraphProp for RangeFiltered<G> {
    type EdgeType = G::EdgeType;
    fn is_directed(&self) -> bool {
        self.graph.is_directed()
    }
}

impl<G: Data> Data for RangeFiltered<G> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}

impl<G: DataMap + NodeIndexable> DataMap for RangeFiltered<G> {
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        let i = self.graph.to_index(id);
        (self.start..self.end)
            .contains(&i)
            .then(|| self.graph.node_weight(id))?
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}
impl<G: GraphBase + NodeIndexable> NodeCount for RangeFiltered<G> {
    fn node_count(&self) -> usize {
        self.end - self.start
    }
}
impl<G: DataValueMap + NodeIndexable> DataValueMap for RangeFiltered<G> {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        let i = self.graph.to_index(id);
        (self.start..self.end)
            .contains(&i)
            .then(|| self.graph.node_weight(id))?
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}

impl<G: NodeIndexable> NodeIndexable for RangeFiltered<G> {
    fn from_index(&self, i: usize) -> Self::NodeId {
        self.graph.from_index(i + self.start)
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        self.graph.to_index(a) - self.start
    }
    fn node_bound(&self) -> usize {
        self.end - self.start
    }
}
impl<G: NodeIndexable> NodeCompactIndexable for RangeFiltered<G> {}

impl<G: GetAdjacencyMatrix> GetAdjacencyMatrix for RangeFiltered<G> {
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        self.graph.is_adjacent(matrix, a, b)
    }
}

impl<'a, G: Data> IntoNodeIdentifiers for &'a RangeFiltered<G>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId> + NodeIndexable,
{
    type NodeIdentifiers =
        iter::NodeIdFilter<<&'a G as IntoNodeIdentifiers>::NodeIdentifiers, &'a G>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::NodeIdFilter(
            self.graph.node_identifiers(),
            &self.graph,
            self.start,
            self.end,
        )
    }
}

impl<'a, G: Data> IntoNodeReferences for &'a RangeFiltered<G>
where
    &'a G: IntoNodeReferences<NodeWeight = G::NodeWeight, NodeId = G::NodeId> + NodeIndexable,
{
    type NodeRef = <&'a G as IntoNodeReferences>::NodeRef;
    type NodeReferences = iter::NodeRefFilter<<&'a G as IntoNodeReferences>::NodeReferences, &'a G>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeRefFilter(
            self.graph.node_references(),
            &self.graph,
            self.start,
            self.end,
        )
    }
}

impl<'a, G: Data> IntoEdgeReferences for &'a RangeFiltered<G>
where
    &'a G: IntoEdgeReferences<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type EdgeRef = <&'a G as IntoEdgeReferences>::EdgeRef;
    type EdgeReferences = iter::EdgeRefFilter<<&'a G as IntoEdgeReferences>::EdgeReferences, &'a G>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeRefFilter(
            self.graph.edge_references(),
            &self.graph,
            self.start,
            self.end,
        )
    }
}

impl<'a, G: Data> IntoEdges for &'a RangeFiltered<G>
where
    &'a G: IntoEdges<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type Edges = iter::EdgeRefFilter<<&'a G as IntoEdges>::Edges, &'a G>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::EdgeRefFilter(self.graph.edges(a), &self.graph, self.start, self.end)
    }
}

impl<'a, G: Data> IntoNeighbors for &'a RangeFiltered<G>
where
    &'a G: IntoNeighbors<EdgeId = G::EdgeId, NodeId = G::NodeId> + NodeIndexable,
{
    type Neighbors = iter::NodeIdFilter<<&'a G as IntoNeighbors>::Neighbors, &'a G>;
    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::NodeIdFilter(self.graph.neighbors(a), &self.graph, self.start, self.end)
    }
}

impl<'a, G: Data> IntoNeighborsDirected for &'a RangeFiltered<G>
where
    &'a G: IntoNeighborsDirected<EdgeId = G::EdgeId, NodeId = G::NodeId> + NodeIndexable,
{
    type NeighborsDirected =
        iter::NodeIdFilter<<&'a G as IntoNeighborsDirected>::NeighborsDirected, &'a G>;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        iter::NodeIdFilter(
            self.graph.neighbors_directed(a, dir),
            &self.graph,
            self.start,
            self.end,
        )
    }
}

impl<'a, G: Data> IntoEdgesDirected for &'a RangeFiltered<G>
where
    &'a G: IntoEdgesDirected<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type EdgesDirected = iter::EdgeRefFilter<<&'a G as IntoEdgesDirected>::EdgesDirected, &'a G>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgeRefFilter(
            self.graph.edges_directed(a, dir),
            &self.graph,
            self.start,
            self.end,
        )
    }
}

#[doc(hidden)]
pub mod iter {
    use super::*;

    pub struct NodeIdFilter<I, G>(pub I, pub G, pub usize, pub usize);
    impl<G: NodeIndexable, I: Iterator<Item = G::NodeId>> Iterator for NodeIdFilter<I, G> {
        type Item = G::NodeId;
        fn next(&mut self) -> Option<G::NodeId> {
            self.0
                .by_ref()
                .find(|&i| (self.2..self.3).contains(&self.1.to_index(i)))
        }
    }

    pub struct NodeRefFilter<I, G>(pub I, pub G, pub usize, pub usize);
    impl<G: IntoNodeReferences + NodeIndexable, I: Iterator<Item = G::NodeRef>> Iterator
        for NodeRefFilter<I, G>
    {
        type Item = G::NodeRef;
        fn next(&mut self) -> Option<G::NodeRef> {
            self.0
                .by_ref()
                .find(|i| (self.2..self.3).contains(&self.1.to_index(i.id())))
        }
    }

    pub struct EdgeRefFilter<I, G>(pub I, pub G, pub usize, pub usize);
    impl<G: IntoEdgeReferences + NodeIndexable, I: Iterator<Item = G::EdgeRef>> Iterator
        for EdgeRefFilter<I, G>
    {
        type Item = G::EdgeRef;
        fn next(&mut self) -> Option<G::EdgeRef> {
            self.0.by_ref().find(|i| {
                let range = self.2..self.3;
                range.contains(&self.1.to_index(i.source()))
                    && range.contains(&self.1.to_index(i.target()))
            })
        }
    }
}
