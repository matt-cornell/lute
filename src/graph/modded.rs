use super::misc::DataValueMap;
use petgraph::data::*;
use petgraph::prelude::Direction;
use petgraph::visit::*;
use small_map::ASmallMap;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct ModdedGraph<G: GraphBase + Data, const N: usize> {
    /// Base graph
    pub graph: G,
    /// Modified nodes
    pub mods: ASmallMap<N, G::NodeId, G::NodeWeight>,
}

impl<G: GraphBase + Data, const N: usize> ModdedGraph<G, N>
where
    G::NodeId: Hash + Eq,
{
    pub fn new(graph: G) -> Self {
        Self {
            graph,
            mods: ASmallMap::new(),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub enum NodeReference<'a, R: NodeRef>
where
    R::NodeId: Copy,
{
    Original(R),
    Modded(R::NodeId, &'a R::Weight),
}
impl<R: NodeRef> Clone for NodeReference<'_, R>
where
    R::NodeId: Copy,
{
    fn clone(&self) -> Self {
        *self
    }
}
impl<R: NodeRef> Copy for NodeReference<'_, R> where R::NodeId: Copy {}
impl<R: NodeRef> NodeRef for NodeReference<'_, R>
where
    R::NodeId: Copy,
{
    type NodeId = R::NodeId;
    type Weight = R::Weight;

    fn id(&self) -> Self::NodeId {
        match *self {
            Self::Original(id) => id.id(),
            Self::Modded(id, _) => id,
        }
    }
    fn weight(&self) -> &Self::Weight {
        match *self {
            Self::Original(ref id) => id.weight(),
            Self::Modded(_, w) => w,
        }
    }
}

impl<G: GraphBase + Data, const N: usize> GraphBase for ModdedGraph<G, N> {
    type NodeId = G::NodeId;
    type EdgeId = G::EdgeId;
}
impl<G: GraphBase + Data + GraphProp, const N: usize> GraphProp for ModdedGraph<G, N> {
    type EdgeType = G::EdgeType;

    fn is_directed(&self) -> bool {
        self.graph.is_directed()
    }
}
impl<G: GraphBase + Data, const N: usize> Data for ModdedGraph<G, N> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}
impl<G: DataMap, const N: usize> DataMap for ModdedGraph<G, N>
where
    G::NodeId: Hash + Eq,
{
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        self.mods.get(&id).or_else(|| self.graph.node_weight(id))
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}
impl<G: DataValueMap, const N: usize> DataValueMap for ModdedGraph<G, N>
where
    G::NodeId: Hash + Eq,
    G::NodeWeight: Copy,
    G::EdgeWeight: Copy,
{
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        self.mods
            .get(&id)
            .copied()
            .or_else(|| self.graph.node_weight(id))
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        self.graph.edge_weight(id)
    }
}
impl<G: NodeCount + Data, const N: usize> NodeCount for ModdedGraph<G, N> {
    fn node_count(&self) -> usize {
        self.graph.node_count()
    }
}
impl<G: EdgeCount + Data, const N: usize> EdgeCount for ModdedGraph<G, N> {
    fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
}
impl<G: NodeIndexable + Data, const N: usize> NodeIndexable for ModdedGraph<G, N> {
    fn node_bound(&self) -> usize {
        self.graph.node_bound()
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        self.graph.to_index(a)
    }
    fn from_index(&self, i: usize) -> Self::NodeId {
        self.graph.from_index(i)
    }
}
impl<G: EdgeIndexable + Data, const N: usize> EdgeIndexable for ModdedGraph<G, N> {
    fn edge_bound(&self) -> usize {
        self.graph.edge_bound()
    }
    fn to_index(&self, a: Self::EdgeId) -> usize {
        self.graph.to_index(a)
    }
    fn from_index(&self, i: usize) -> Self::EdgeId {
        self.graph.from_index(i)
    }
}
impl<G: NodeCompactIndexable + Data, const N: usize> NodeCompactIndexable for ModdedGraph<G, N> {}

impl<G: GetAdjacencyMatrix + Data, const N: usize> GetAdjacencyMatrix for ModdedGraph<G, N> {
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        self.graph.is_adjacent(matrix, a, b)
    }
}

impl<'a, G: Data, const N: usize> IntoNodeIdentifiers for &'a ModdedGraph<G, N>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
{
    type NodeIdentifiers = <&'a G as IntoNodeIdentifiers>::NodeIdentifiers;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        self.graph.node_identifiers()
    }
}
impl<'a, G: Data, const N: usize> IntoNodeReferences for &'a ModdedGraph<G, N>
where
    &'a G: IntoNodeReferences<NodeId = G::NodeId, NodeWeight = G::NodeWeight>,
    G::NodeId: Hash + Eq,
{
    type NodeRef = NodeReference<'a, <&'a G as IntoNodeReferences>::NodeRef>;
    type NodeReferences = iter::NodeReferences<'a, G, N>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeReferences {
            iter: self.graph.node_references(),
            graph: self,
        }
    }
}
impl<'a, G: Data, const N: usize> IntoEdgeReferences for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdgeReferences<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgeRef = <&'a G as IntoEdgeReferences>::EdgeRef;
    type EdgeReferences = <&'a G as IntoEdgeReferences>::EdgeReferences;

    fn edge_references(self) -> Self::EdgeReferences {
        self.graph.edge_references()
    }
}

impl<'a, G: Data, const N: usize> IntoEdges for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdges<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type Edges = <&'a G as IntoEdges>::Edges;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        self.graph.edges(a)
    }
}
impl<'a, G: Data, const N: usize> IntoEdgesDirected for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdgesDirected<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgesDirected = <&'a G as IntoEdgesDirected>::EdgesDirected;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        self.graph.edges_directed(a, dir)
    }
}
impl<'a, G: Data, const N: usize> IntoNeighbors for &'a ModdedGraph<G, N>
where
    &'a G: IntoNeighbors<NodeId = G::NodeId>,
{
    type Neighbors = <&'a G as IntoNeighbors>::Neighbors;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        self.graph.neighbors(a)
    }
}

impl<'a, G: Data, const N: usize> IntoNeighborsDirected for &'a ModdedGraph<G, N>
where
    &'a G: IntoNeighborsDirected<NodeId = G::NodeId>,
{
    type NeighborsDirected = <&'a G as IntoNeighborsDirected>::NeighborsDirected;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        self.graph.neighbors_directed(a, dir)
    }
}

pub mod iter {
    use super::*;

    pub struct NodeReferences<'a, G: GraphBase + Data, const N: usize>
    where
        &'a G: IntoNodeReferences,
    {
        pub(super) iter: <&'a G as IntoNodeReferences>::NodeReferences,
        pub(super) graph: &'a ModdedGraph<G, N>,
    }
    impl<'a, G: GraphBase + Data, const N: usize> Iterator for NodeReferences<'a, G, N>
    where
        &'a G: IntoNodeReferences<NodeId = G::NodeId, NodeWeight = G::NodeWeight>,
        G::NodeId: Hash + Eq,
    {
        type Item = NodeReference<'a, <&'a G as IntoNodeReferences>::NodeRef>;

        fn next(&mut self) -> Option<Self::Item> {
            let r = self.iter.next()?;
            let i = r.id();
            if let Some(w) = self.graph.mods.get(&r.id()) {
                Some(NodeReference::Modded(i, w))
            } else {
                Some(NodeReference::Original(r))
            }
        }
    }
}
