use super::*;
use petgraph::visit::*;
use std::iter::Map;
use std::ops::Range;

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

impl<Ix: IndexType, R> Data for Molecule<Ix, R> {
    type NodeWeight = Atom;
    type EdgeWeight = Bond;
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeIdentifiers for Molecule<Ix, R> {
    type NodeIdentifiers = Map<Range<usize>, fn(usize) -> NodeIndex<Ix>>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        (0..self.node_count()).map(|ix| NodeIndex(Ix::new(ix)))
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeReferences for Molecule<Ix, R> {
    type NodeRef = NodeReference<Ix>;
    type NodeReferences = NodeReferences<Ix, R>;

    fn node_references(self) -> Self::NodeReferences {
        NodeReferences {
            range: 0..self.node_bound(),
            mol_idx: self.index,
            arena: self.arena,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NodeReferences<Ix, R> {
    range: Range<usize>,
    mol_idx: Ix,
    arena: R,
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for NodeReferences<Ix, R> {
    type Item = NodeReference<Ix>;

    fn next(&mut self) -> Option<Self::Item> {
        self.range
            .next()
            .map(|i| NodeReference::new(self.mol_idx, NodeIndex(Ix::new(i)), self.arena))
    }
}
