use crate::graph::misc::DataValueMap;
use itertools::Itertools;
use petgraph::data::*;
use petgraph::graph::*;
use petgraph::visit::*;
use petgraph::{Direction, EdgeType};
use std::ops::{Index, IndexMut};

/// A trait for types acting like `Option`. Mostly here for its special impl for `Atom`, which uses a scratch bit.
pub trait Optional {
    type Inner;

    fn is_some(&self) -> bool;
    /// Create a version of this type that's empty
    fn none() -> Self;
    /// Create a version of this type
    fn some(val: Self::Inner) -> Self;

    /// Get a reference to the inner value, panicking if this is in an "empty" state.
    fn unwrap_ref(&self) -> &Self::Inner;

    /// Get a mutable reference to the inner value, panicking if this is in an "empty" state.
    fn unwrap_mut(&mut self) -> &mut Self::Inner;
}
impl Optional for crate::core::Atom {
    type Inner = Self;

    fn is_some(&self) -> bool {
        self.data.scratch() & (1 << 5) == 0
    }
    fn none() -> Self {
        Self::new_scratch(0, 1 << 5)
    }
    fn some(mut val: Self::Inner) -> Self {
        val.map_scratch(|s| s & ((1 << 5) - 1));
        val
    }
    fn unwrap_ref(&self) -> &Self::Inner {
        debug_assert!(self.is_some(), "\"unwrapped\" an \"empty\" molecule!");
        self
    }
    fn unwrap_mut(&mut self) -> &mut Self::Inner {
        debug_assert!(self.is_some(), "\"unwrapped\" an \"empty\" molecule!");
        self
    }
}
impl<T> Optional for Option<T> {
    type Inner = T;

    fn is_some(&self) -> bool {
        Option::is_some(self)
    }
    fn none() -> Self {
        None
    }
    fn some(val: T) -> Self {
        Some(val)
    }
    fn unwrap_ref(&self) -> &Self::Inner {
        self.as_ref().unwrap()
    }
    fn unwrap_mut(&mut self) -> &mut Self::Inner {
        self.as_mut().unwrap()
    }
}

/// A "semi-sparse" graph. Node indices are stable across removals, but unlike `StableGraph`, you can also allocate a contiguous range of nodes.
#[derive(Debug, Clone)]
pub struct SemiSparseGraph<N, E, Ty: EdgeType, Ix: IndexType> {
    /// The underlying graph we use for storage
    inner: Graph<N, E, Ty, Ix>,
    /// Holes we've created in the graph
    holes: [[Ix; 2]; 8],
    node_count: usize,
}
impl<N, E, Ty: EdgeType, Ix: IndexType> SemiSparseGraph<N, E, Ty, Ix> {
    pub fn new() -> Self {
        Self::with_capacity(0, 0)
    }

    /// Create a new semi-sparse graph with a given capacity
    pub fn with_capacity(nodes: usize, edges: usize) -> Self {
        Self {
            inner: Graph::with_capacity(nodes, edges),
            holes: [[IndexType::max(); 2]; 8],
            node_count: 0,
        }
    }

    /// Create a new semi-sparse graph from an underlying graph, assuming that all elements are in the "filled" state.
    pub fn from_compact(graph: Graph<N, E, Ty, Ix>) -> Self {
        Self {
            node_count: graph.node_count(),
            inner: graph,
            holes: [[IndexType::max(); 2]; 8],
        }
    }

    /// Clear all of the nodes and edges in this graph.
    pub fn clear(&mut self) {
        self.inner.clear();
        self.holes = [[IndexType::max(); 2]; 8];
        self.node_count = 0;
    }

    /// Get a reference to the underlying graph.
    pub fn graph(&self) -> &Graph<N, E, Ty, Ix> {
        &self.inner
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> SemiSparseGraph<N, E, Ty, Ix> {
    /// Allocate a contiguous range of nodes, returns the index of the first one.
    pub fn allocate_range(&mut self, size: usize, mut fill: impl FnMut() -> N::Inner) -> usize {
        if size == 0 {
            return 0;
        }
        self.node_count += size;
        let best = self
            .holes
            .iter()
            .enumerate()
            .filter_map(|(n, &r @ [from, to])| {
                let sz = to.index() - from.index();
                (sz >= size).then_some((n, r, sz))
            })
            .min_by_key(|x| size.abs_diff(x.2 * 2));
        if let Some((i, [from, _to], sz)) = best {
            let start = from.index();
            for i in start..(start + size) {
                let elem = self.inner.node_weight_mut(Ix::new(i).into()).unwrap();
                debug_assert!(!elem.is_some(), "overwrite of element detected!");
                *elem = N::some(fill());
            }
            if sz > size {
                self.holes[i][0] = Ix::new(start + size);
            } else {
                self.holes[i..].rotate_left(1);
                self.update_holes(1);
            }
            for i in start..(start + size) {
                self.inner[NodeIndex::new(i)] = N::some(fill());
            }
            start
        } else {
            let start = self.inner.node_count();
            self.inner.reserve_nodes(size);
            for _ in 0..size {
                self.inner.add_node(N::some(fill()));
            }
            start
        }
    }

    /// Free a range of nodes. It doesn't have to be the same as the one passed in.
    pub fn free_range(&mut self, start: usize, end: usize) {
        for i in start..end {
            let old = std::mem::replace(&mut self.inner[NodeIndex::new(i)], N::none());
            if !old.is_some() {
                self.node_count -= 1;
            }
            self.detach_node(i);
        }
        let mut merge_to: Option<&mut Ix> = None;
        let mut rotate = None;
        for (i, [s, e]) in self.holes.iter_mut().enumerate() {
            if s.index() == end {
                if let Some(s2) = merge_to {
                    *s2 = *s;
                    rotate = Some(i);
                    break;
                } else {
                    *s = Ix::new(start);
                    merge_to = Some(e);
                }
            } else if e.index() == start {
                if let Some(e2) = merge_to {
                    *e2 = *e;
                    rotate = Some(i);
                    break;
                } else {
                    *e = Ix::new(end);
                    merge_to = Some(s);
                }
            }
        }
        if let Some(i) = rotate {
            self.holes[i..].rotate_left(1);
        }
        self.update_holes(1);
    }

    /// Remove all edges from a node.
    fn detach_node(&mut self, idx: usize) {
        let mut next = self.inner.raw_nodes()[idx].next_edge(Direction::Outgoing);
        while next != EdgeIndex::end() {
            next = self.inner.raw_nodes()[idx].next_edge(Direction::Outgoing);
            self.inner.remove_edge(next);
        }
        next = self.inner.raw_nodes()[idx].next_edge(Direction::Incoming);
        while next != EdgeIndex::end() {
            next = self.inner.raw_nodes()[idx].next_edge(Direction::Incoming);
            self.inner.remove_edge(next);
        }
    }

    /// Update our hole list, with last being the number of elements from the end that we need to check
    fn update_holes(&mut self, last: usize) {
        let hlen = self.holes.len();
        let start_ix = self.holes[hlen - last - 1][0].index();
        let empty = [IndexType::max(); 2];
        if start_ix >= self.inner.node_count() {
            self.holes[(hlen - last)..].fill(empty);
            return;
        }
        let chunks = self.inner.raw_nodes()[start_ix..]
            .iter()
            .enumerate()
            .chunk_by(|n| !n.1.weight.is_some());
        let mut iter = chunks.into_iter().filter_map(|(p, mut g)| {
            p.then(|| -> Option<_> {
                let start = g.next()?.0;
                let end = g.last()?.0;
                Some([start + start_ix, end + start_ix])
            })
            .flatten()
        });
        for hole in &mut self.holes[(hlen - last)..] {
            *hole = iter.next().map_or(empty, |arr| arr.map(Ix::new));
        }
    }

    /// Add a node to the graph, the same way it does in the base `Graph`.
    pub fn add_node(&mut self, weight: N::Inner) -> NodeIndex<Ix> {
        let mut val = Some(weight);
        NodeIndex::new(self.allocate_range(1, || val.take().unwrap()))
    }
    /// Add an edge to the graph, the same way it does in the base `Graph`.
    pub fn add_edge(&mut self, a: NodeIndex<Ix>, b: NodeIndex<Ix>, weight: E) -> EdgeIndex<Ix> {
        #[cfg(debug_assertions)]
        {
            debug_assert!(self.inner.raw_nodes()[a.index()].weight.is_some());
            debug_assert!(self.inner.raw_nodes()[b.index()].weight.is_some());
        }
        self.inner.add_edge(a, b, weight)
    }
    /// Find an edge between two nodes.
    pub fn find_edge(&self, a: NodeIndex<Ix>, b: NodeIndex<Ix>) -> Option<EdgeIndex<Ix>> {
        self.inner.raw_nodes()[a.index()]
            .weight
            .is_some()
            .then_some(())?;
        self.inner.raw_nodes()[b.index()]
            .weight
            .is_some()
            .then_some(())?;
        self.inner.find_edge(a, b)
    }
}
impl<N, E, Ty: EdgeType, Ix: IndexType> Default for SemiSparseGraph<N, E, Ty, Ix> {
    fn default() -> Self {
        Self::new()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> Index<NodeIndex<Ix>>
    for SemiSparseGraph<N, E, Ty, Ix>
{
    type Output = N::Inner;

    fn index(&self, index: NodeIndex<Ix>) -> &Self::Output {
        self.inner[index].unwrap_ref()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> IndexMut<NodeIndex<Ix>>
    for SemiSparseGraph<N, E, Ty, Ix>
{
    fn index_mut(&mut self, index: NodeIndex<Ix>) -> &mut Self::Output {
        self.inner[index].unwrap_mut()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> Index<EdgeIndex<Ix>>
    for SemiSparseGraph<N, E, Ty, Ix>
{
    type Output = E;

    fn index(&self, index: EdgeIndex<Ix>) -> &Self::Output {
        &self.inner[index]
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> IndexMut<EdgeIndex<Ix>>
    for SemiSparseGraph<N, E, Ty, Ix>
{
    fn index_mut(&mut self, index: EdgeIndex<Ix>) -> &mut Self::Output {
        &mut self.inner[index]
    }
}

impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> GraphBase for SemiSparseGraph<N, E, Ty, Ix> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> GraphProp for SemiSparseGraph<N, E, Ty, Ix> {
    type EdgeType = Ty;
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> Data for SemiSparseGraph<N, E, Ty, Ix> {
    type NodeWeight = N::Inner;
    type EdgeWeight = E;
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> NodeIndexable for SemiSparseGraph<N, E, Ty, Ix> {
    fn node_bound(&self) -> usize {
        self.inner.node_bound()
    }
    fn from_index(&self, i: usize) -> Self::NodeId {
        NodeIndex::new(i)
    }
    fn to_index(&self, a: Self::NodeId) -> usize {
        a.index()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> EdgeIndexable for SemiSparseGraph<N, E, Ty, Ix> {
    fn edge_bound(&self) -> usize {
        self.inner.edge_bound()
    }
    fn from_index(&self, i: usize) -> Self::EdgeId {
        EdgeIndex::new(i)
    }
    fn to_index(&self, a: Self::EdgeId) -> usize {
        a.index()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> NodeCount for SemiSparseGraph<N, E, Ty, Ix> {
    fn node_count(&self) -> usize {
        self.node_count
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> EdgeCount for SemiSparseGraph<N, E, Ty, Ix> {
    fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> DataMap for SemiSparseGraph<N, E, Ty, Ix> {
    fn node_weight(&self, id: Self::NodeId) -> Option<&Self::NodeWeight> {
        let w = self.inner.node_weight(id)?;
        w.is_some().then(|| w.unwrap_ref())
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<&Self::EdgeWeight> {
        self.inner.edge_weight(id)
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> DataMapMut for SemiSparseGraph<N, E, Ty, Ix> {
    fn node_weight_mut(&mut self, id: Self::NodeId) -> Option<&mut Self::NodeWeight> {
        let w = self.inner.node_weight_mut(id)?;
        w.is_some().then(|| w.unwrap_mut())
    }
    fn edge_weight_mut(&mut self, id: Self::EdgeId) -> Option<&mut Self::EdgeWeight> {
        self.inner.edge_weight_mut(id)
    }
}
impl<N: Optional, E: Copy, Ty: EdgeType, Ix: IndexType> DataValueMap
    for SemiSparseGraph<N, E, Ty, Ix>
where
    N::Inner: Copy,
{
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        let w = self.inner.node_weight(id)?;
        w.is_some().then(|| *w.unwrap_ref())
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        self.inner.edge_weight(id).copied()
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> Visitable for SemiSparseGraph<N, E, Ty, Ix> {
    type Map = <Graph<N, E, Ty, Ix> as Visitable>::Map;

    fn visit_map(&self) -> Self::Map {
        self.inner.visit_map()
    }
    fn reset_map(&self, map: &mut Self::Map) {
        self.inner.reset_map(map);
    }
}

impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoNodeIdentifiers
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type NodeIdentifiers = iter::NodeIdFilter<NodeIndices<Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::NodeIdFilter(self.inner.node_identifiers(), &self.inner)
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoNodeReferences
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type NodeRef = (NodeIndex<Ix>, &'a N::Inner);
    type NodeReferences = iter::NodeRefFilter<NodeReferences<'a, N, Ix>>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeRefFilter(self.inner.node_references())
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoEdgeReferences
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type EdgeRef = EdgeReference<'a, E, Ix>;
    type EdgeReferences = iter::EdgeRefFilter<EdgeReferences<'a, E, Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeRefFilter(self.inner.edge_references(), &self.inner)
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoNeighbors
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type Neighbors = iter::NodeIdFilter<Neighbors<'a, E, Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::NodeIdFilter(self.inner.neighbors(a), &self.inner)
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoNeighborsDirected
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type NeighborsDirected = iter::NodeIdFilter<Neighbors<'a, E, Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::Neighbors {
        iter::NodeIdFilter(self.inner.neighbors_directed(a, dir), &self.inner)
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoEdges
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type Edges = iter::EdgeRefFilter<Edges<'a, E, Ty, Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::EdgeRefFilter(self.inner.edges(a), &self.inner)
    }
}
impl<'a, N: Optional, E, Ty: EdgeType, Ix: IndexType> IntoEdgesDirected
    for &'a SemiSparseGraph<N, E, Ty, Ix>
{
    type EdgesDirected = iter::EdgeRefFilter<Edges<'a, E, Ty, Ix>, &'a Graph<N, E, Ty, Ix>>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgeRefFilter(self.inner.edges_directed(a, dir), &self.inner)
    }
}

#[doc(hidden)]
pub mod iter {
    use super::*;
    use petgraph::data::DataMap;

    #[derive(Debug, Clone, PartialEq)]
    pub struct NodeIdFilter<I, G>(pub I, pub G);
    impl<G: DataMap, I: Iterator<Item = G::NodeId>> Iterator for NodeIdFilter<I, G>
    where
        G::NodeWeight: Optional,
    {
        type Item = I::Item;
        fn next(&mut self) -> Option<Self::Item> {
            self.0
                .find(|&e| self.1.node_weight(e).map_or(false, G::NodeWeight::is_some))
        }
    }

    #[derive(Debug, Clone, PartialEq)]
    pub struct NodeRefFilter<I>(pub I);
    impl<'a, Ix, N: Optional + 'a, I: Iterator<Item = (NodeIndex<Ix>, &'a N)>> Iterator
        for NodeRefFilter<I>
    {
        type Item = (NodeIndex<Ix>, &'a N::Inner);
        fn next(&mut self) -> Option<Self::Item> {
            self.0
                .find_map(|(i, e)| e.is_some().then(|| (i, e.unwrap_ref())))
        }
    }

    #[derive(Debug, Clone, PartialEq)]
    pub struct EdgeRefFilter<I, G>(pub I, pub G);
    impl<
            'a,
            Ix: IndexType,
            E: 'a,
            G: DataMap<NodeId = NodeIndex<Ix>>,
            I: Iterator<Item = petgraph::graph::EdgeReference<'a, E, Ix>>,
        > Iterator for EdgeRefFilter<I, G>
    where
        G::NodeWeight: Optional,
    {
        type Item = EdgeReference<'a, E, Ix>;
        fn next(&mut self) -> Option<Self::Item> {
            self.0.find_map(|e| {
                self.1.node_weight(e.source())?.is_some().then_some(())?;
                self.1.node_weight(e.target())?.is_some().then_some(())?;
                Some(e)
            })
        }
    }
}
