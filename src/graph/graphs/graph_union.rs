use crate::prelude::DataValueMap;
use petgraph::data::*;
use petgraph::visit::*;
use petgraph::Direction;

/// A union of multiple graphs that acts as one, disconnected graph.
#[derive(Debug, Default, Clone)]
pub struct GraphUnion<G>(pub Vec<G>);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Tagged<I> {
    graph_id: usize,
    inner: I,
}
impl<I: NodeRef> NodeRef for Tagged<I> {
    type NodeId = Tagged<I::NodeId>;
    type Weight = I::Weight;

    fn id(&self) -> Self::NodeId {
        Tagged {
            graph_id: self.graph_id,
            inner: self.inner.id(),
        }
    }
    fn weight(&self) -> &Self::Weight {
        self.inner.weight()
    }
}
impl<I: EdgeRef> EdgeRef for Tagged<I> {
    type NodeId = Tagged<I::NodeId>;
    type EdgeId = Tagged<I::EdgeId>;
    type Weight = I::Weight;

    fn id(&self) -> Self::EdgeId {
        Tagged {
            graph_id: self.graph_id,
            inner: self.inner.id(),
        }
    }
    fn source(&self) -> Self::NodeId {
        Tagged {
            graph_id: self.graph_id,
            inner: self.inner.source(),
        }
    }
    fn target(&self) -> Self::NodeId {
        Tagged {
            graph_id: self.graph_id,
            inner: self.inner.target(),
        }
    }
    fn weight(&self) -> &Self::Weight {
        self.inner.weight()
    }
}

impl<G: GraphBase> GraphBase for GraphUnion<G> {
    type NodeId = Tagged<G::NodeId>;
    type EdgeId = Tagged<G::EdgeId>;
}
impl<G: GraphProp> GraphProp for GraphUnion<G> {
    type EdgeType = G::EdgeType;
}
impl<G: Data> Data for GraphUnion<G> {
    type NodeWeight = G::NodeWeight;
    type EdgeWeight = G::EdgeWeight;
}
impl<G: DataMap> DataMap for GraphUnion<G> {
    fn node_weight(&self, id: Tagged<G::NodeId>) -> Option<&G::NodeWeight> {
        self.0[id.graph_id].node_weight(id.inner)
    }
    fn edge_weight(&self, id: Tagged<G::EdgeId>) -> Option<&G::EdgeWeight> {
        self.0[id.graph_id].edge_weight(id.inner)
    }
}
impl<G: DataMapMut> DataMapMut for GraphUnion<G> {
    fn node_weight_mut(&mut self, id: Tagged<G::NodeId>) -> Option<&mut G::NodeWeight> {
        self.0[id.graph_id].node_weight_mut(id.inner)
    }
    fn edge_weight_mut(&mut self, id: Tagged<G::EdgeId>) -> Option<&mut G::EdgeWeight> {
        self.0[id.graph_id].edge_weight_mut(id.inner)
    }
}
impl<G: DataValueMap> DataValueMap for GraphUnion<G> {
    fn node_weight(&self, id: Tagged<G::NodeId>) -> Option<G::NodeWeight> {
        self.0[id.graph_id].node_weight(id.inner)
    }
    fn edge_weight(&self, id: Tagged<G::EdgeId>) -> Option<G::EdgeWeight> {
        self.0[id.graph_id].edge_weight(id.inner)
    }
}
impl<G: NodeCount> NodeCount for GraphUnion<G> {
    fn node_count(&self) -> usize {
        self.0.iter().map(G::node_count).sum()
    }
}
impl<G: EdgeCount> EdgeCount for GraphUnion<G> {
    fn edge_count(&self) -> usize {
        self.0.iter().map(G::edge_count).sum()
    }
}
impl<G: NodeIndexable> NodeIndexable for GraphUnion<G> {
    fn to_index(&self, a: Self::NodeId) -> usize {
        (0..a.graph_id)
            .map(|i| self.0[i].node_bound())
            .sum::<usize>()
            + self.0[a.graph_id].to_index(a.inner)
    }
    fn from_index(&self, mut i: usize) -> Self::NodeId {
        let mut graph_id = 0;
        for g in &self.0 {
            if let Some(i2) = i.checked_sub(g.node_bound()) {
                i = i2;
                graph_id += 1;
            } else {
                break;
            }
        }
        Self::NodeId {
            graph_id,
            inner: self.0[graph_id].from_index(i),
        }
    }
    fn node_bound(&self) -> usize {
        self.0.iter().map(G::node_bound).sum()
    }
}
impl<G: EdgeIndexable> EdgeIndexable for GraphUnion<G> {
    fn to_index(&self, a: Self::EdgeId) -> usize {
        (0..a.graph_id)
            .map(|i| self.0[i].edge_bound())
            .sum::<usize>()
            + self.0[a.graph_id].to_index(a.inner)
    }
    fn from_index(&self, mut i: usize) -> Self::EdgeId {
        let mut graph_id = 0;
        for g in &self.0 {
            if let Some(i2) = i.checked_sub(g.edge_bound()) {
                i = i2;
                graph_id += 1;
            } else {
                break;
            }
        }
        Self::EdgeId {
            graph_id,
            inner: self.0[graph_id].from_index(i),
        }
    }
    fn edge_bound(&self) -> usize {
        self.0.iter().map(G::edge_bound).sum()
    }
}
impl<G: NodeCompactIndexable> NodeCompactIndexable for GraphUnion<G> {}

impl<'a, G: GraphBase> IntoNodeIdentifiers for &'a GraphUnion<G>
where
    &'a G: IntoNodeIdentifiers<NodeId = G::NodeId>,
{
    type NodeIdentifiers = iter::TaggedIter<'a, G, <&'a G as IntoNodeIdentifiers>::NodeIdentifiers>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        iter::TaggedIter {
            graphs: &self.0,
            index: 0,
            iter: self.0.first().map(IntoNodeIdentifiers::node_identifiers),
            map: IntoNodeIdentifiers::node_identifiers,
        }
    }
}
impl<'a, G: Data> IntoNodeReferences for &'a GraphUnion<G>
where
    &'a G: IntoNodeReferences<NodeId = G::NodeId, NodeWeight = G::NodeWeight>,
{
    type NodeRef = Tagged<<&'a G as IntoNodeReferences>::NodeRef>;
    type NodeReferences = iter::TaggedIter<'a, G, <&'a G as IntoNodeReferences>::NodeReferences>;

    fn node_references(self) -> Self::NodeReferences {
        iter::TaggedIter {
            graphs: &self.0,
            index: 0,
            iter: self.0.first().map(IntoNodeReferences::node_references),
            map: IntoNodeReferences::node_references,
        }
    }
}
impl<'a, G: Data> IntoEdgeReferences for &'a GraphUnion<G>
where
    &'a G: IntoEdgeReferences<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgeRef = Tagged<<&'a G as IntoEdgeReferences>::EdgeRef>;
    type EdgeReferences = iter::TaggedIter<'a, G, <&'a G as IntoEdgeReferences>::EdgeReferences>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::TaggedIter {
            graphs: &self.0,
            index: 0,
            iter: self.0.first().map(IntoEdgeReferences::edge_references),
            map: IntoEdgeReferences::edge_references,
        }
    }
}
impl<'a, G: GraphBase> IntoNeighbors for &'a GraphUnion<G>
where
    &'a G: IntoNeighbors<NodeId = G::NodeId>,
{
    type Neighbors = iter::MappedIter<<&'a G as IntoNeighbors>::Neighbors>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::MappedIter {
            graph_id: a.graph_id,
            iter: self.0[a.graph_id].neighbors(a.inner),
        }
    }
}
impl<'a, G: GraphBase> IntoNeighborsDirected for &'a GraphUnion<G>
where
    &'a G: IntoNeighborsDirected<NodeId = G::NodeId>,
{
    type NeighborsDirected = iter::MappedIter<<&'a G as IntoNeighborsDirected>::NeighborsDirected>;

    fn neighbors_directed(self, a: Self::NodeId, dir: Direction) -> Self::NeighborsDirected {
        iter::MappedIter {
            graph_id: a.graph_id,
            iter: self.0[a.graph_id].neighbors_directed(a.inner, dir),
        }
    }
}
impl<'a, G: Data> IntoEdges for &'a GraphUnion<G>
where
    &'a G: IntoEdges<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type Edges = iter::MappedIter<<&'a G as IntoEdges>::Edges>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::MappedIter {
            graph_id: a.graph_id,
            iter: self.0[a.graph_id].edges(a.inner),
        }
    }
}
impl<'a, G: Data> IntoEdgesDirected for &'a GraphUnion<G>
where
    &'a G: IntoEdgesDirected<NodeId = G::NodeId, EdgeId = G::EdgeId, EdgeWeight = G::EdgeWeight>,
{
    type EdgesDirected = iter::MappedIter<<&'a G as IntoEdgesDirected>::EdgesDirected>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::MappedIter {
            graph_id: a.graph_id,
            iter: self.0[a.graph_id].edges_directed(a.inner, dir),
        }
    }
}
impl<G: Visitable> Visitable for GraphUnion<G> {
    type Map = UnionVisitMap<G::Map>;

    fn visit_map(&self) -> Self::Map {
        UnionVisitMap(self.0.iter().map(G::visit_map).collect())
    }
    fn reset_map(&self, map: &mut Self::Map) {
        self.0
            .iter()
            .zip(&mut map.0)
            .for_each(|(g, m)| g.reset_map(m));
    }
}

pub struct UnionVisitMap<M>(Vec<M>);
impl<N, M: VisitMap<N>> VisitMap<Tagged<N>> for UnionVisitMap<M> {
    fn visit(&mut self, a: Tagged<N>) -> bool {
        self.0[a.graph_id].visit(a.inner)
    }
    fn is_visited(&self, a: &Tagged<N>) -> bool {
        self.0[a.graph_id].is_visited(&a.inner)
    }
}

#[doc(hidden)]
pub mod iter {
    use super::*;
    pub struct TaggedIter<'a, G, I> {
        pub(super) graphs: &'a [G],
        pub(super) index: usize,
        pub(super) iter: Option<I>,
        pub(super) map: fn(&'a G) -> I,
    }
    impl<G, I: Iterator> Iterator for TaggedIter<'_, G, I> {
        type Item = Tagged<I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            loop {
                if let Some(elem) = self.iter.as_mut()?.next() {
                    return Some(Tagged {
                        graph_id: self.index,
                        inner: elem,
                    });
                }
                self.index += 1;
                self.iter = self.graphs.get(self.index).map(self.map);
            }
        }
    }
    pub struct MappedIter<I> {
        pub(super) graph_id: usize,
        pub(super) iter: I,
    }
    impl<I: Iterator> Iterator for MappedIter<I> {
        type Item = Tagged<I::Item>;

        fn next(&mut self) -> Option<Self::Item> {
            self.iter.next().map(|inner| Tagged {
                graph_id: self.graph_id,
                inner,
            })
        }
    }
}
