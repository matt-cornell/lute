use super::*;
use petgraph::visit::NodeRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeIndex<Ix>(pub Ix);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeIndex<Ix>(Ix, Ix);
impl<Ix: Copy> EdgeIndex<Ix> {
    #[inline(always)]
    pub const fn source(self) -> Ix {
        self.0
    }
    #[inline(always)]
    pub const fn target(self) -> Ix {
        self.1
    }
    #[inline(always)]
    pub const fn explode(self) -> (Ix, Ix) {
        (self.0, self.1)
    }
}
impl<Ix: IndexType> EdgeIndex<Ix> {
    pub fn new(source: Ix, target: Ix) -> Self {
        if target.index() < source.index() {
            Self(target, source)
        } else {
            Self(source, target)
        }
    }
}

/// Petgraph stuff needs a reference, but we may not have one because of `MolRepr::Modded`. Since
/// `Atom`s are copyable though, this works and isn't much more expensive.
#[derive(Debug, Clone, Copy)]
pub struct NodeReference<Ix> {
    node_idx: NodeIndex<Ix>,
    atom: Atom,
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
