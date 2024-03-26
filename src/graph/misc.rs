use petgraph::data::DataMap;
use petgraph::graph::{Graph, IndexType};
use petgraph::stable_graph::StableGraph;
use petgraph::visit::Data;
use petgraph::EdgeType;

pub trait DataValueMap: Data {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight>;
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight>;
}

impl<G: DataValueMap> DataValueMap for &G {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        G::node_weight(self, id)
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        G::edge_weight(self, id)
    }
}
impl<N: Copy, E: Copy, Ty: EdgeType, Ix: IndexType> DataValueMap for Graph<N, E, Ty, Ix> {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        DataMap::node_weight(&self, id).copied()
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        DataMap::edge_weight(&self, id).copied()
    }
}
impl<N: Copy, E: Copy, Ty: EdgeType, Ix: IndexType> DataValueMap for StableGraph<N, E, Ty, Ix> {
    fn node_weight(&self, id: Self::NodeId) -> Option<Self::NodeWeight> {
        DataMap::node_weight(&self, id).copied()
    }
    fn edge_weight(&self, id: Self::EdgeId) -> Option<Self::EdgeWeight> {
        DataMap::edge_weight(&self, id).copied()
    }
}
