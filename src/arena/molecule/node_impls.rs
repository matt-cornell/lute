use super::*;
use petgraph::visit::{EdgeRef, NodeRef};

/// Node indices are compact indices over the atoms in a molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeIndex<Ix>(pub Ix);
impl<Ix> From<Ix> for NodeIndex<Ix> {
    fn from(value: Ix) -> Self {
        Self(value)
    }
}

/// Edge indices are tuple of the starting and ending nodes.
/// Because the graph is undirected, these are always stored with the lower index first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeIndex<Ix>(Ix, Ix);
impl<Ix: PartialOrd> EdgeIndex<Ix> {
    /// Create a new `EdgeIndex`. This enforces that it's sorted.
    #[inline(always)]
    pub fn new(source: Ix, target: Ix) -> Self {
        if target < source {
            Self(target, source)
        } else {
            Self(source, target)
        }
    }
}
impl<Ix> EdgeIndex<Ix> {
    pub fn rev(self) -> Self {
        Self(self.1, self.0)
    }
    pub fn as_ref(&self) -> EdgeIndex<&Ix> {
        EdgeIndex(&self.0, &self.1)
    }
}
impl<Ix: Copy> EdgeIndex<Ix> {
    /// Get the lower index
    #[inline(always)]
    pub const fn source(self) -> Ix {
        self.0
    }
    /// Get the higher index
    #[inline(always)]
    pub const fn target(self) -> Ix {
        self.1
    }
    /// Explode as a tuple for pattern matching
    #[inline(always)]
    pub const fn explode(self) -> (Ix, Ix) {
        (self.0, self.1)
    }
}
impl<Ix: IndexType> From<(Ix, Ix)> for EdgeIndex<Ix> {
    fn from((a, b): (Ix, Ix)) -> Self {
        Self::new(a, b)
    }
}
impl<Ix: IndexType> From<[Ix; 2]> for EdgeIndex<Ix> {
    fn from([a, b]: [Ix; 2]) -> Self {
        Self::new(a, b)
    }
}
impl<Ix: IndexType> From<(NodeIndex<Ix>, NodeIndex<Ix>)> for EdgeIndex<Ix> {
    fn from((a, b): (NodeIndex<Ix>, NodeIndex<Ix>)) -> Self {
        Self::new(a.0, b.0)
    }
}
impl<Ix: IndexType> From<[NodeIndex<Ix>; 2]> for EdgeIndex<Ix> {
    fn from([a, b]: [NodeIndex<Ix>; 2]) -> Self {
        Self::new(a.0, b.0)
    }
}

/// Petgraph stuff needs a reference, but we can't get one because a `Molecule` doesn't actually hold a reference. Since
/// `Atom`s are copyable though, this works and isn't much more expensive.
#[derive(Debug, Clone, Copy)]
pub struct NodeReference<Ix> {
    pub node_idx: NodeIndex<Ix>,
    pub atom: Atom,
}
impl<Ix> NodeReference<Ix> {
    pub fn new<R: ArenaAccessor<Ix = Ix>>(
        mol_idx: MolIndex<Ix>,
        node_idx: NodeIndex<Ix>,
        arena: R,
    ) -> Self
    where
        Ix: IndexType,
    {
        trace!(
            mi = mol_idx.index(),
            ni = node_idx.0.index(),
            "creating node reference"
        );
        Self {
            node_idx,
            atom: arena
                .get_arena()
                .molecule(mol_idx)
                .get_atom(node_idx)
                .unwrap(),
        }
    }
    pub fn with_weight(node_idx: NodeIndex<Ix>, weight: Atom) -> Self {
        Self {
            node_idx,
            atom: weight,
        }
    }
    pub const fn atom(&self) -> &Atom {
        &self.atom
    }
}
impl<Ix: Copy> NodeRef for NodeReference<Ix> {
    type NodeId = NodeIndex<Ix>;
    type Weight = Atom;

    fn id(&self) -> NodeIndex<Ix> {
        self.node_idx
    }
    fn weight(&self) -> &Atom {
        &self.atom
    }
}

#[derive(Debug, Clone, Copy)]
pub struct EdgeReference<Ix> {
    pub edge_idx: EdgeIndex<Ix>,
    pub bond: Bond,
    pub target_first: bool,
}
impl<Ix: PartialOrd> EdgeReference<Ix> {
    pub fn new<R: ArenaAccessor<Ix = Ix>>(
        mol_idx: MolIndex<Ix>,
        source: Ix,
        target: Ix,
        arena: R,
    ) -> Self
    where
        Ix: IndexType,
    {
        let edge_idx = EdgeIndex::new(source, target);
        Self {
            bond: arena
                .get_arena()
                .molecule(mol_idx)
                .get_bond(edge_idx)
                .unwrap(),
            target_first: target < source,
            edge_idx,
        }
    }
    pub fn with_weight(source: Ix, target: Ix, weight: Bond) -> Self {
        Self {
            target_first: target < source,
            edge_idx: EdgeIndex::new(source, target),
            bond: weight,
        }
    }
    pub const fn bond(&self) -> &Bond {
        &self.bond
    }
    pub fn rev(self) -> Self {
        Self {
            target_first: !self.target_first,
            ..self
        }
    }
}
impl<Ix: Copy> EdgeRef for EdgeReference<Ix> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
    type Weight = Bond;

    fn id(&self) -> EdgeIndex<Ix> {
        self.edge_idx
    }
    fn weight(&self) -> &Bond {
        &self.bond
    }
    fn source(&self) -> NodeIndex<Ix> {
        NodeIndex(if self.target_first {
            self.edge_idx.target()
        } else {
            self.edge_idx.source()
        })
    }
    fn target(&self) -> NodeIndex<Ix> {
        NodeIndex(if self.target_first {
            self.edge_idx.source()
        } else {
            self.edge_idx.target()
        })
    }
}
