use super::misc::DataValueMap;
use hybridmap::HybridMap;
use itertools::Itertools;
use petgraph::data::*;
use petgraph::prelude::Direction;
use petgraph::visit::*;
use petgraph::Undirected;
use smallvec::SmallVec;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct ModdedGraph<G: GraphBase + Data, const N: usize> {
    /// Base graph
    pub graph: G,
    /// Modified nodes
    pub mods: HybridMap<G::NodeId, G::NodeWeight, N>,
    /// Additional nodes- this should be sorted, and contains the node it's attached to, the new
    /// weight, and the edge weight connecting it
    pub additional: SmallVec<(G::NodeId, G::NodeWeight, G::EdgeWeight), N>,
}

impl<G: GraphBase + Data, const N: usize> ModdedGraph<G, N>
where
    G::NodeId: Hash + Eq,
{
    pub fn new(graph: G) -> Self {
        Self {
            graph,
            mods: HybridMap::new(),
            additional: SmallVec::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NodeId<N> {
    Original(N),
    Added(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EdgeId<E> {
    Original(E),
    Added(usize),
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub enum NodeReference<'a, R: NodeRef>
where
    R::NodeId: Copy,
{
    Original(R),
    Modded(R::NodeId, &'a R::Weight),
    Added(usize, &'a R::Weight),
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
    type NodeId = NodeId<R::NodeId>;
    type Weight = R::Weight;

    fn id(&self) -> Self::NodeId {
        match *self {
            Self::Original(id) => NodeId::Original(id.id()),
            Self::Modded(id, _) => NodeId::Original(id),
            Self::Added(id, _) => NodeId::Added(id),
        }
    }
    fn weight(&self) -> &Self::Weight {
        match *self {
            Self::Original(ref id) => id.weight(),
            Self::Modded(_, w) | Self::Added(_, w) => w,
        }
    }
}
#[derive(Debug, PartialEq, Eq, Hash)]
pub enum EdgeReference<'a, R: EdgeRef> {
    Original(R),
    Added {
        source: usize,
        target: R::NodeId,
        weight: &'a R::Weight,
    },
}
impl<R: EdgeRef> Clone for EdgeReference<'_, R>
where
    R::NodeId: Copy,
{
    fn clone(&self) -> Self {
        *self
    }
}
impl<R: EdgeRef> Copy for EdgeReference<'_, R> where R::NodeId: Copy {}
impl<R: EdgeRef> EdgeRef for EdgeReference<'_, R>
where
    R::NodeId: Copy,
{
    type NodeId = NodeId<R::NodeId>;
    type EdgeId = EdgeId<R::EdgeId>;
    type Weight = R::Weight;

    fn id(&self) -> Self::EdgeId {
        match *self {
            Self::Original(id) => EdgeId::Original(id.id()),
            Self::Added { source, .. } => EdgeId::Added(source),
        }
    }
    fn source(&self) -> Self::NodeId {
        match *self {
            Self::Original(id) => NodeId::Original(id.source()),
            Self::Added { source, .. } => NodeId::Added(source),
        }
    }
    fn target(&self) -> Self::NodeId {
        match *self {
            Self::Original(id) => NodeId::Original(id.target()),
            Self::Added { target, .. } => NodeId::Original(target),
        }
    }
    fn weight(&self) -> &Self::Weight {
        match *self {
            Self::Original(ref id) => id.weight(),
            Self::Added { weight, .. } => weight,
        }
    }
}

impl<G: GraphBase + Data, const N: usize> GraphBase for ModdedGraph<G, N> {
    type NodeId = NodeId<G::NodeId>;
    type EdgeId = EdgeId<G::EdgeId>;
}
impl<G: GraphBase + Data + GraphProp<EdgeType = Undirected>, const N: usize> GraphProp
    for ModdedGraph<G, N>
{
    type EdgeType = Undirected;

    fn is_directed(&self) -> bool {
        false
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
        match id {
            NodeId::Original(n) => self.mods.get(&n).or_else(|| self.graph.node_weight(n)),
            NodeId::Added(i) => self.additional.get(i).map(|n| &n.1),
        }
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        match id {
            EdgeId::Original(n) => self.graph.edge_weight(n),
            EdgeId::Added(i) => self.additional.get(i).map(|n| &n.2),
        }
    }
}
impl<G: DataValueMap, const N: usize> DataValueMap for ModdedGraph<G, N>
where
    G::NodeId: Hash + Eq,
    G::NodeWeight: Copy,
    G::EdgeWeight: Copy,
{
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        match id {
            NodeId::Original(n) => self
                .mods
                .get(&n)
                .copied()
                .or_else(|| self.graph.node_weight(n)),
            NodeId::Added(i) => self.additional.get(i).map(|n| n.1),
        }
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        match id {
            EdgeId::Original(n) => self.graph.edge_weight(n),
            EdgeId::Added(i) => self.additional.get(i).map(|n| n.2),
        }
    }
}
impl<G: NodeCount + Data, const N: usize> NodeCount for ModdedGraph<G, N> {
    fn node_count(&self) -> usize {
        self.graph.node_count() + self.additional.len()
    }
}
impl<G: EdgeCount + Data, const N: usize> EdgeCount for ModdedGraph<G, N> {
    fn edge_count(&self) -> usize {
        self.graph.edge_count() + self.additional.len()
    }
}
impl<G: NodeIndexable + Data, const N: usize> NodeIndexable for ModdedGraph<G, N> {
    fn node_bound(&self) -> usize {
        self.graph.node_bound() + self.additional.len()
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        match a {
            NodeId::Original(n) => self.graph.to_index(n),
            NodeId::Added(i) => i + self.graph.node_bound(),
        }
    }
    fn from_index(&self, i: usize) -> Self::NodeId {
        i.checked_sub(self.graph.node_bound())
            .map_or_else(|| NodeId::Original(self.graph.from_index(i)), NodeId::Added)
    }
}
impl<G: EdgeIndexable + Data, const N: usize> EdgeIndexable for ModdedGraph<G, N> {
    fn edge_bound(&self) -> usize {
        self.graph.edge_bound() + self.additional.len()
    }
    fn to_index(&self, a: Self::EdgeId) -> usize {
        match a {
            EdgeId::Original(n) => self.graph.to_index(n),
            EdgeId::Added(i) => i + self.graph.edge_bound(),
        }
    }
    fn from_index(&self, i: usize) -> Self::EdgeId {
        i.checked_sub(self.graph.edge_bound())
            .map_or_else(|| EdgeId::Original(self.graph.from_index(i)), EdgeId::Added)
    }
}
impl<G: NodeCompactIndexable + Data, const N: usize> NodeCompactIndexable for ModdedGraph<G, N> {}

impl<G: GetAdjacencyMatrix + Data, const N: usize> GetAdjacencyMatrix for ModdedGraph<G, N> {
    type AdjMatrix = G::AdjMatrix;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        self.graph.adjacency_matrix()
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        match (a, b) {
            (NodeId::Original(a), NodeId::Original(b)) => self.graph.is_adjacent(matrix, a, b),
            (NodeId::Original(o), NodeId::Added(i)) | (NodeId::Added(i), NodeId::Original(o)) => {
                self.additional.get(i).map_or(false, |a| a.0 == o)
            }
            _ => false,
        }
    }
}

impl<'a, G: Data, const N: usize> IntoNodeIdentifiers for &'a ModdedGraph<G, N>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
{
    type NodeIdentifiers = iter::NodeIdentifiers<<&'a G as IntoNodeIdentifiers>::NodeIdentifiers>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::NodeIdentifiers(
            iter::IterState::Original(self.graph.node_identifiers()),
            self.additional.len(),
        )
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
            iter: iter::IterState::Original(self.graph.node_references()),
            graph: self,
        }
    }
}
impl<'a, G: Data, const N: usize> IntoEdgeReferences for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdgeReferences<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgeRef = EdgeReference<'a, <&'a G as IntoEdgeReferences>::EdgeRef>;
    type EdgeReferences = iter::EdgeReferences<'a, G, N>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeReferences {
            iter: iter::IterState::Original(self.graph.edge_references()),
            graph: self,
        }
    }
}

impl<'a, G: Data, const N: usize> IntoEdges for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdges<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type Edges = iter::Edges<'a, G, <&'a G as IntoEdges>::Edges, N>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        match a {
            NodeId::Original(a) => {
                let start = self
                    .additional
                    .iter()
                    .position(|i| i.0 == a)
                    .unwrap_or(usize::MAX);
                iter::Edges::Inner(iter::EdgesInner {
                    iter: self.graph.edges(a),
                    range: start,
                    graph: self,
                })
            }
            NodeId::Added(i) => {
                if let Some(&(n, _, ref w)) = self.additional.get(i) {
                    iter::Edges::Single(EdgeReference::Added {
                        source: i,
                        target: n,
                        weight: w,
                    })
                } else {
                    iter::Edges::Exhausted
                }
            }
        }
    }
}
impl<'a, G: Data, const N: usize> IntoEdgesDirected for &'a ModdedGraph<G, N>
where
    &'a G: IntoEdgesDirected<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgesDirected = iter::Edges<'a, G, <&'a G as IntoEdgesDirected>::EdgesDirected, N>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        match a {
            NodeId::Original(a) => {
                let start = if dir == Direction::Outgoing {
                    self.additional
                        .iter()
                        .position(|i| i.0 == a)
                        .unwrap_or(usize::MAX)
                } else {
                    usize::MAX
                };
                iter::Edges::Inner(iter::EdgesInner {
                    iter: self.graph.edges_directed(a, dir),
                    range: start,
                    graph: self,
                })
            }
            NodeId::Added(i) => {
                if dir == Direction::Incoming {
                    if let Some(&(n, _, ref w)) = self.additional.get(i) {
                        iter::Edges::Single(EdgeReference::Added {
                            source: i,
                            target: n,
                            weight: w,
                        })
                    } else {
                        iter::Edges::Exhausted
                    }
                } else {
                    iter::Edges::Exhausted
                }
            }
        }
    }
}
impl<'a, G: Data, const N: usize> IntoNeighbors for &'a ModdedGraph<G, N>
where
    &'a G: IntoNeighbors<NodeId = G::NodeId>,
{
    type Neighbors = iter::Neighbors<<&'a G as IntoNeighbors>::Neighbors>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        match a {
            NodeId::Original(a) => {
                let mut it = self.additional.iter().positions(|i| i.0 == a);
                let (range, end) = if let Some(start) = it.next() {
                    (start, it.next().unwrap_or(self.additional.len()))
                } else {
                    (usize::MAX, usize::MAX)
                };
                iter::Neighbors::Inner(iter::NeighborsInner {
                    iter: self.graph.neighbors(a),
                    range,
                    end,
                })
            }
            NodeId::Added(i) => {
                if let Some(r) = self.additional.get(i) {
                    iter::Neighbors::Single(r.0)
                } else {
                    iter::Neighbors::Exhausted
                }
            }
        }
    }
}

impl<'a, G: Data, const N: usize> IntoNeighborsDirected for &'a ModdedGraph<G, N>
where
    &'a G: IntoNeighborsDirected<NodeId = G::NodeId>,
{
    type NeighborsDirected = iter::Neighbors<<&'a G as IntoNeighborsDirected>::NeighborsDirected>;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        match a {
            NodeId::Original(a) => {
                let (range, end) = if dir == Direction::Outgoing {
                    let mut it = self.additional.iter().positions(|i| i.0 == a);
                    if let Some(start) = it.next() {
                        (start, it.next().unwrap_or(self.additional.len()))
                    } else {
                        (usize::MAX, usize::MAX)
                    }
                } else {
                    (usize::MAX, usize::MAX)
                };
                iter::Neighbors::Inner(iter::NeighborsInner {
                    iter: self.graph.neighbors_directed(a, dir),
                    range,
                    end,
                })
            }
            NodeId::Added(i) => {
                if dir == Direction::Incoming {
                    if let Some(r) = self.additional.get(i) {
                        iter::Neighbors::Single(r.0)
                    } else {
                        iter::Neighbors::Exhausted
                    }
                } else {
                    iter::Neighbors::Exhausted
                }
            }
        }
    }
}

pub mod iter {
    use super::*;

    #[derive(Debug, Clone)]
    pub(super) enum IterState<I> {
        Original(I),
        Range(usize),
    }
    #[derive(Debug, Clone)]
    pub struct NodeIdentifiers<I>(pub(super) IterState<I>, pub(super) usize);
    impl<I: Iterator> Iterator for NodeIdentifiers<I> {
        type Item = NodeId<I::Item>;
        fn next(&mut self) -> Option<Self::Item> {
            match &mut self.0 {
                IterState::Original(it) => {
                    if let Some(ret) = it.next() {
                        Some(NodeId::Original(ret))
                    } else if self.1 == 0 {
                        self.0 = IterState::Range(0);
                        None
                    } else {
                        self.0 = IterState::Range(1);
                        Some(NodeId::Added(0))
                    }
                }
                IterState::Range(i) => (*i < self.1).then(|| {
                    *i += 1;
                    NodeId::Added(*i - 1)
                }),
            }
        }
    }
    pub struct NodeReferences<'a, G: GraphBase + Data, const N: usize>
    where
        &'a G: IntoNodeReferences,
    {
        pub(super) iter: IterState<<&'a G as IntoNodeReferences>::NodeReferences>,
        pub(super) graph: &'a ModdedGraph<G, N>,
    }
    impl<'a, G: GraphBase + Data, const N: usize> Iterator for NodeReferences<'a, G, N>
    where
        &'a G: IntoNodeReferences<NodeId = G::NodeId, NodeWeight = G::NodeWeight>,
        G::NodeId: Hash + Eq,
    {
        type Item = NodeReference<'a, <&'a G as IntoNodeReferences>::NodeRef>;

        fn next(&mut self) -> Option<Self::Item> {
            match &mut self.iter {
                IterState::Original(it) => {
                    if let Some(ret) = it.next() {
                        let id = ret.id();
                        Some(
                            self.graph
                                .mods
                                .get(&id)
                                .map_or(NodeReference::Original(ret), |r| {
                                    NodeReference::Modded(id, r)
                                }),
                        )
                    } else if let Some(r) = self.graph.additional.first() {
                        self.iter = IterState::Range(1);
                        Some(NodeReference::Added(0, &r.1))
                    } else {
                        self.iter = IterState::Range(0);
                        None
                    }
                }
                IterState::Range(i) => self.graph.additional.get(*i).map(|r| {
                    *i += 1;
                    NodeReference::Added(*i - 1, &r.1)
                }),
            }
        }
    }
    pub struct EdgeReferences<'a, G: GraphBase + Data, const N: usize>
    where
        &'a G: IntoEdgeReferences,
    {
        pub(super) iter: IterState<<&'a G as IntoEdgeReferences>::EdgeReferences>,
        pub(super) graph: &'a ModdedGraph<G, N>,
    }
    impl<'a, G: GraphBase + Data, const N: usize> Iterator for EdgeReferences<'a, G, N>
    where
        &'a G:
            IntoEdgeReferences<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
    {
        type Item = EdgeReference<'a, <&'a G as IntoEdgeReferences>::EdgeRef>;

        fn next(&mut self) -> Option<Self::Item> {
            match &mut self.iter {
                IterState::Original(it) => {
                    if let Some(ret) = it.next() {
                        Some(EdgeReference::Original(ret))
                    } else if let Some(&(n, _, ref w)) = self.graph.additional.first() {
                        self.iter = IterState::Range(1);
                        Some(EdgeReference::Added {
                            source: 0,
                            target: n,
                            weight: w,
                        })
                    } else {
                        self.iter = IterState::Range(0);
                        None
                    }
                }
                IterState::Range(i) => self.graph.additional.get(*i).map(|&(n, _, ref w)| {
                    *i += 1;
                    EdgeReference::Added {
                        source: *i - 1,
                        target: n,
                        weight: w,
                    }
                }),
            }
        }
    }
    pub struct EdgesInner<'a, G: GraphBase + Data, I, const N: usize> {
        pub(super) iter: I,
        pub(super) range: usize,
        pub(super) graph: &'a ModdedGraph<G, N>,
    }
    impl<'a, G: GraphBase + Data, I: Iterator, const N: usize> Iterator for EdgesInner<'a, G, I, N>
    where
        I::Item: EdgeRef<NodeId = G::NodeId, EdgeId = G::EdgeId, Weight = G::EdgeWeight>,
    {
        type Item = EdgeReference<'a, I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(ret) = self.iter.next() {
                Some(EdgeReference::Original(ret))
            } else if self.range == usize::MAX {
                None
            } else {
                let (n, _, ref w) = self.graph.additional[self.range];
                self.range += 1;
                if self
                    .graph
                    .additional
                    .get(self.range)
                    .map_or(true, |r| r.0 != n)
                {
                    self.range = usize::MAX;
                }
                Some(EdgeReference::Added {
                    source: self.range,
                    target: n,
                    weight: w,
                })
            }
        }
    }
    pub enum Edges<'a, G: GraphBase + Data, I: Iterator, const N: usize>
    where
        I::Item: EdgeRef,
    {
        Inner(EdgesInner<'a, G, I, N>),
        Single(EdgeReference<'a, I::Item>),
        Exhausted,
    }
    impl<'a, G: GraphBase + Data, I: Iterator, const N: usize> Iterator for Edges<'a, G, I, N>
    where
        I::Item: EdgeRef<NodeId = G::NodeId, EdgeId = G::EdgeId, Weight = G::EdgeWeight>,
    {
        type Item = EdgeReference<'a, I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            if let Self::Inner(it) = self {
                it.next()
            } else if let Self::Single(r) = std::mem::replace(self, Self::Exhausted) {
                Some(r)
            } else {
                None
            }
        }
    }
    pub struct NeighborsInner<I> {
        pub(super) iter: I,
        pub(super) range: usize,
        pub(super) end: usize,
    }
    impl<I: Iterator> Iterator for NeighborsInner<I>
    where
        I::Item: Copy,
    {
        type Item = NodeId<I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(ret) = self.iter.next() {
                Some(NodeId::Original(ret))
            } else if self.range == self.end {
                None
            } else {
                self.range += 1;
                Some(NodeId::Added(self.range - 1))
            }
        }
    }
    pub enum Neighbors<I: Iterator>
    where
        I::Item: Copy,
    {
        Inner(NeighborsInner<I>),
        Single(I::Item),
        Exhausted,
    }
    impl<I: Iterator> Iterator for Neighbors<I>
    where
        I::Item: Copy,
    {
        type Item = NodeId<I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            if let Self::Inner(it) = self {
                it.next()
            } else if let Self::Single(r) = std::mem::replace(self, Self::Exhausted) {
                Some(NodeId::Original(r))
            } else {
                None
            }
        }
    }
}
