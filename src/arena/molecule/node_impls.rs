use super::*;
use petgraph::visit::{EdgeRef, NodeRef};
use std::hash::{Hash, Hasher};

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
#[derive(Debug, Clone, Copy)]
pub struct EdgeIndex<Ix>(pub Ix, pub Ix);
impl<Ix> EdgeIndex<Ix> {
    /// Create a new `EdgeIndex`. This enforces that it's sorted.
    #[inline(always)]
    pub fn new(source: Ix, target: Ix) -> Self {
        Self(source, target)
    }
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
impl<Ix: Ord> PartialEq for EdgeIndex<Ix> {
    fn eq(&self, other: &Self) -> bool {
        let EdgeIndex(mut a, mut b) = self.as_ref();
        if a.cmp(b) != other.0.cmp(&other.1) {
            std::mem::swap(&mut a, &mut b)
        }
        *a == other.0 && *b == other.1
    }
}
impl<Ix: Ord> Eq for EdgeIndex<Ix> {}
impl<Ix: Ord> PartialOrd for EdgeIndex<Ix> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<Ix: Ord> Ord for EdgeIndex<Ix> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let EdgeIndex(mut sa, mut sb) = self.as_ref();
        let EdgeIndex(mut oa, mut ob) = other.as_ref();
        if sa > sb {
            std::mem::swap(&mut sa, &mut sb)
        }
        if oa > ob {
            std::mem::swap(&mut oa, &mut ob)
        }
        sa.cmp(sb).then(oa.cmp(ob))
    }
}
impl<Ix: Ord + Hash> Hash for EdgeIndex<Ix> {
    fn hash<H: Hasher>(&self, h: &mut H) {
        let EdgeIndex(mut a, mut b) = self.as_ref();
        if a > b {
            std::mem::swap(&mut a, &mut b)
        }
        a.hash(h);
        b.hash(h);
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
    pub fn new<R: ArenaAccessor<Ix = Ix>>(mol_idx: Ix, node_idx: NodeIndex<Ix>, arena: R) -> Self
    where
        Ix: IndexType,
    {
        println!("mi: {}, ni: {}", mol_idx.index(), node_idx.0.index());
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
}
impl<Ix> EdgeReference<Ix> {
    pub fn new<R: ArenaAccessor<Ix = Ix>>(mol_idx: Ix, edge_idx: EdgeIndex<Ix>, arena: R) -> Self
    where
        Ix: IndexType,
    {
        Self {
            edge_idx,
            bond: arena
                .get_arena()
                .molecule(mol_idx)
                .get_bond(edge_idx)
                .unwrap(),
        }
    }
    pub fn with_weight(edge_idx: EdgeIndex<Ix>, weight: Bond) -> Self {
        Self {
            edge_idx,
            bond: weight,
        }
    }
    pub const fn bond(&self) -> &Bond {
        &self.bond
    }
    pub fn rev(self) -> Self {
        let EdgeReference { edge_idx, bond } = self;
        Self {
            edge_idx: edge_idx.rev(),
            bond,
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
        NodeIndex(self.edge_idx.source())
    }
    fn target(&self) -> NodeIndex<Ix> {
        NodeIndex(self.edge_idx.target())
    }
}
