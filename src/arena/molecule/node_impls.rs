use super::*;
use petgraph::visit::NodeRef;

/// Node indices are compact indices over the atoms in a molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeIndex<Ix>(pub Ix);

/// Edge indices are tuple of the starting and ending nodes.
/// Because the graph is undirected, these are always stored with the lower index first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeIndex<Ix>(Ix, Ix);
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
impl<Ix: IndexType> EdgeIndex<Ix> {
    /// Create a new `EdgeIndex`. This enforces that it's sorted.
    pub fn new(source: Ix, target: Ix) -> Self {
        if target.index() < source.index() {
            Self(target, source)
        } else {
            Self(source, target)
        }
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
        Self {
            node_idx,
            atom: arena
                .get_arena()
                .molecule(mol_idx)
                .get_atom(node_idx)
                .unwrap(),
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
