use super::misc::DataValueMap;
use crate::utils::bitset::BitSet;
use num_traits::PrimInt;
use petgraph::data::*;
use petgraph::visit::*;
use petgraph::Direction;
use std::fmt::{self, Binary, Debug, Formatter};

#[derive(Clone)]
pub struct BitFiltered<G, T, const N: usize, const C: bool> {
    pub graph: G,
    pub filter: BitSet<T, N>,
}
impl<G, T, const N: usize, const C: bool> BitFiltered<G, T, N, C> {
    pub const fn new(graph: G, filter: BitSet<T, N>) -> Self {
        Self { graph, filter }
    }
}

impl<G: Debug, T: Binary, const N: usize, const C: bool> Debug for BitFiltered<G, T, N, C> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("BitFiltered")
            .field("graph", &self.graph)
            .field("filter", &self.filter)
            .field("compact", &C)
            .finish()
    }
}

impl<G: GraphBase, T, const N: usize, const C: bool> GraphBase for BitFiltered<G, T, N, C> {
    type EdgeId = G::EdgeId;
    type NodeId = G::NodeId;
}

impl<G: GraphProp, T, const N: usize, const C: bool> GraphProp for BitFiltered<G, T, N, C> {
    type EdgeType = G::EdgeType;
    fn is_directed(&self) -> bool {
        self.graph.is_directed()
    }
}

impl<G: Data, T, const N: usize, const C: bool> Data for BitFiltered<G, T, N, C> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}

impl<G: DataMap + NodeIndexable, T: PrimInt, const N: usize, const C: bool> DataMap
    for BitFiltered<G, T, N, C>
{
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        self.filter
            .get(self.graph.to_index(id))
            .then(|| self.graph.node_weight(id))?
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}
impl<G: GraphBase, T: PrimInt, const N: usize, const C: bool> NodeCount
    for BitFiltered<G, T, N, C>
{
    fn node_count(&self) -> usize {
        self.filter.count_ones()
    }
}
impl<G: DataValueMap + NodeIndexable, T: PrimInt, const N: usize, const C: bool> DataValueMap
    for BitFiltered<G, T, N, C>
{
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        self.filter
            .get(self.graph.to_index(id))
            .then(|| self.graph.node_weight(id))?
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}

impl<G: NodeIndexable, T, const N: usize> NodeIndexable for BitFiltered<G, T, N, false> {
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
impl<G: NodeIndexable, T: PrimInt, const N: usize> NodeIndexable for BitFiltered<G, T, N, true> {
    fn from_index(&self, i: usize) -> Self::NodeId {
        self.graph.from_index(self.filter.nth(i).unwrap())
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        self.filter.index(self.graph.to_index(a)).unwrap()
    }
    fn node_bound(&self) -> usize {
        self.filter.count_ones()
    }
}
impl<G: NodeIndexable, T: PrimInt, const N: usize> NodeCompactIndexable
    for BitFiltered<G, T, N, true>
{
}

impl<G: GetAdjacencyMatrix, T, const N: usize, const C: bool> GetAdjacencyMatrix
    for BitFiltered<G, T, N, C>
{
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        self.graph.is_adjacent(matrix, a, b)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoNodeIdentifiers
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId> + NodeIndexable,
{
    type NodeIdentifiers =
        iter::NodeIdFilter<'a, <&'a G as IntoNodeIdentifiers>::NodeIdentifiers, &'a G, T, { N }>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::NodeIdFilter(self.graph.node_identifiers(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoNodeReferences
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoNodeReferences<NodeWeight = G::NodeWeight, NodeId = G::NodeId> + NodeIndexable,
{
    type NodeRef = <&'a G as IntoNodeReferences>::NodeRef;
    type NodeReferences =
        iter::NodeRefFilter<'a, <&'a G as IntoNodeReferences>::NodeReferences, &'a G, T, { N }>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeRefFilter(self.graph.node_references(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoEdgeReferences
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoEdgeReferences<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type EdgeRef = <&'a G as IntoEdgeReferences>::EdgeRef;
    type EdgeReferences =
        iter::EdgeRefFilter<'a, <&'a G as IntoEdgeReferences>::EdgeReferences, &'a G, T, { N }>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeRefFilter(self.graph.edge_references(), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoEdges
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoEdges<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type Edges = iter::EdgeRefFilter<'a, <&'a G as IntoEdges>::Edges, &'a G, T, { N }>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::EdgeRefFilter(self.graph.edges(a), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoNeighbors
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoNeighbors<EdgeId = G::EdgeId, NodeId = G::NodeId> + NodeIndexable,
{
    type Neighbors = iter::NodeIdFilter<'a, <&'a G as IntoNeighbors>::Neighbors, &'a G, T, { N }>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::NodeIdFilter(self.graph.neighbors(a), &self.graph, &self.filter)
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoNeighborsDirected
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoNeighborsDirected<EdgeId = G::EdgeId, NodeId = G::NodeId> + NodeIndexable,
{
    type NeighborsDirected = iter::NodeIdFilter<
        'a,
        <&'a G as IntoNeighborsDirected>::NeighborsDirected,
        &'a G,
        T,
        { N },
    >;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        iter::NodeIdFilter(
            self.graph.neighbors_directed(a, dir),
            &self.graph,
            &self.filter,
        )
    }
}

impl<'a, G: Data, T: PrimInt, const N: usize, const C: bool> IntoEdgesDirected
    for &'a BitFiltered<G, T, N, C>
where
    &'a G: IntoEdgesDirected<EdgeId = G::EdgeId, NodeId = G::NodeId, EdgeWeight = G::EdgeWeight>
        + NodeIndexable,
{
    type EdgesDirected =
        iter::EdgeRefFilter<'a, <&'a G as IntoEdgesDirected>::EdgesDirected, &'a G, T, { N }>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgeRefFilter(self.graph.edges_directed(a, dir), &self.graph, &self.filter)
    }
}

#[doc(hidden)]
pub mod iter {
    use super::*;

    pub struct NodeIdFilter<'a, I, G, T, const N: usize>(pub I, pub G, pub &'a BitSet<T, N>);
    impl<G: NodeIndexable, I: Iterator<Item = G::NodeId>, T: PrimInt, const N: usize> Iterator
        for NodeIdFilter<'_, I, G, T, N>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0.by_ref().find(|&i| self.2.get(self.1.to_index(i)))
        }
    }

    pub struct NodeRefFilter<'a, I, G, T, const N: usize>(pub I, pub G, pub &'a BitSet<T, N>);
    impl<
            G: NodeIndexable + IntoNodeReferences,
            I: Iterator<Item = G::NodeRef>,
            T: PrimInt,
            const N: usize,
        > Iterator for NodeRefFilter<'_, I, G, T, N>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0
                .by_ref()
                .find(|i| self.2.get(self.1.to_index(i.id())))
        }
    }

    pub struct EdgeRefFilter<'a, I, G, T, const N: usize>(pub I, pub G, pub &'a BitSet<T, N>);
    impl<
            G: NodeIndexable + IntoEdgeReferences,
            I: Iterator<Item = G::EdgeRef>,
            T: PrimInt,
            const N: usize,
        > Iterator for EdgeRefFilter<'_, I, G, T, N>
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<I::Item> {
            self.0.by_ref().find(|i| {
                self.2.get(self.1.to_index(i.source())) && self.2.get(self.1.to_index(i.target()))
            })
        }
    }
}
