use super::*;
use petgraph::visit::NodeRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NodeIndex<Ix> {
    ix: Ix,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeIndex<Ix> {
    pub from: NodeIndex<Ix>,
    pub to: NodeIndex<Ix>,
}

/// Petgraph stuff needs a reference, but we may not have one because of `MolRepr::Modded`. Since
/// `Atom`s are copyable though, this works and isn't much more expensive.
#[derive(Debug, Clone, Copy)]
pub struct NodeReference<Ix> {
    mol_idx: Ix,
    node_idx: NodeIndex<Ix>,
    atom: Atom,
}
impl<Ix> NodeReference<Ix> {
    pub fn new<R: ArenaAccessor<Ix = Ix>>(mol_idx: Ix, node_idx: NodeIndex<Ix>, arena: R) -> Self
    where
        Ix: IndexType,
    {
        Self {
            mol_idx,
            node_idx,
            atom: arena.get_arena().molecule(mol_idx).get_atom(node_idx),
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
