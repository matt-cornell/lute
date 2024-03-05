use super::*;
use petgraph::visit::*;

impl<Ix: Copy + PartialEq, R> GraphBase for Molecule<Ix, R> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
}

impl<Ix: Copy + PartialEq, R: Copy> GraphRef for Molecule<Ix, R> {}

impl<Ix: Copy + PartialEq, R> GraphProp for Molecule<Ix, R> {
    type EdgeType = petgraph::Undirected;
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> NodeCount for Molecule<Ix, R> {
    fn node_count(&self) -> usize {
        self.arena.get_arena().parts[self.index.index()].1.index()
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> NodeIndexable for Molecule<Ix, R> {
    fn node_bound(&self) -> usize {
        self.arena.get_arena().parts[self.index.index()].1.index()
    }
    fn to_index(&self, a: NodeIndex<Ix>) -> usize {
        a.0.index()
    }
    fn from_index(&self, i: usize) -> NodeIndex<Ix> {
        NodeIndex(Ix::new(i))
    }
}
