//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use super::*;
use crate::graph::*;
use arena::access::RefAcc;
use hybridmap::HybridMap;
use petgraph::graph::DefaultIx;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallvec::SmallVec;
use std::collections::VecDeque;
use std::fmt::Debug;
use std::hash::Hash;

#[macro_export]
macro_rules! arena {
    (of $ty:ty: $($mol:expr),* $(,)?) => {
        {
            let mut arena = Arena::<$ty>::new();
            $(arena.insert_mol(&$mol);)*
            arena
        }
    };
    ($($mol:expr),* $(,)?) => {
        {
            let mut arena = Arena::new();
            $(arena.insert_mol(&$mol);)*
            arena
        }
    };
}

const ATOM_BIT_STORAGE: usize = 2;

type Graph<Ix> = StableGraph<Atom, Bond, Undirected, Ix>;
type BSType = crate::utils::bitset::BitSet<usize, ATOM_BIT_STORAGE>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(C)]
pub(crate) struct InterFragBond<Ix> {
    pub an: Ix,
    pub bn: Ix,
    pub ai: Ix,
    pub bi: Ix,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BrokenMol<Ix> {
    pub frags: SmallVec<Ix, 4>,
    pub bonds: SmallVec<InterFragBond<Ix>, 4>,
}

#[derive(Debug, Clone)]
pub(crate) struct ModdedMol<Ix> {
    pub base: Ix,
    pub patch: HybridMap<Ix, Atom, 4>,
}
impl<Ix: IndexType> PartialEq for ModdedMol<Ix> {
    fn eq(&self, other: &Self) -> bool {
        self.base == other.base && {
            let mut lhs = self.patch.iter().collect::<SmallVec<_, 4>>();
            let mut rhs = other.patch.iter().collect::<SmallVec<_, 4>>();
            lhs.sort_by_key(|x| x.0);
            rhs.sort_by_key(|x| x.0);
            lhs == rhs
        }
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum MolRepr<Ix: IndexType> {
    Atomic(BSType),
    Broken(BrokenMol<Ix>),
    Modify(ModdedMol<Ix>),
    Redirect(Ix),
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena<Ix: IndexType = DefaultIx> {
    graph: Graph<Ix>,
    pub(crate) parts: SmallVec<(MolRepr<Ix>, Ix), 16>,
}
impl<Ix: IndexType> Arena<Ix> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
    }

    /// Get a reference to the `parts` field. For debugging purposes only.
    pub fn expose_parts(&self) -> impl Debug + Copy + '_ {
        &self.parts
    }

    #[inline(always)]
    fn push_frag(&mut self, frag: (MolRepr<Ix>, Ix)) -> Ix {
        let max = <Ix as IndexType>::max().index();
        let idx = self.parts.len();
        assert!(
            idx < max,
            "too many fragments in molecule: limit of {idx} reached!"
        );
        self.parts.push(frag);
        Ix::new(idx)
    }

    fn contains_group_impl(&self, mol: Ix, group: Ix, seen: &mut BSType) -> bool {
        if mol == group {
            return true;
        }
        if seen.get(mol.index()) {
            return false;
        }
        seen.set(mol.index(), true);
        match self.parts.get(mol.index()) {
            Some((MolRepr::Broken(b), _)) => b
                .frags
                .iter()
                .any(|f| self.contains_group_impl(*f as _, group, seen)),
            Some((MolRepr::Redirect(r), _)) => self.contains_group_impl(*r, group, seen),
            _ => false,
        }
    }

    /// Check if `mol` contains `group`
    pub fn contains_group(&self, mol: Ix, mut group: Ix) -> bool {
        while let Some((MolRepr::Redirect(r), _)) = self.parts.get(group.index()) {
            group = *r;
        }
        let mut seen = BSType::with_capacity(self.parts.len());
        self.contains_group_impl(mol, group, &mut seen)
    }

    /// Get a graph of the molecule at the given index. Note that `Molecule::from_arena` could give
    /// better results as it can borrow from `RefCell`s and `RwLock`s.
    pub fn molecule(&self, mol: Ix) -> Molecule<Ix, access::RefAcc<Ix>> {
        Molecule::from_arena(self, mol)
    }

    /// Insert a molecule into the arena, deduplicating common parts. May misbehave if multiple
    /// disjoint molecules are in the graph.
    pub fn insert_mol<G>(&mut self, mol: G) -> Ix
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + misc::DataValueMap
            + GraphProp<EdgeType = Undirected>
            + GraphRef
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
        G::NodeId: Hash + Eq,
    {
        let max = <Ix as IndexType>::max().index();
        assert!(
            mol.node_count() < max,
            "molecule has too many atoms: {}, max is {max}",
            mol.node_count()
        );
        let mut frags = {
            let (mut frags, mut src) = self
                .parts
                .iter()
                .enumerate()
                .filter_map(|(n, frag)| -> Option<(_, &[_])> {
                    match &frag.0 {
                        MolRepr::Atomic(_) => Some((n, &[])),
                        MolRepr::Broken(b) => Some((n, &b.frags)),
                        _ => None,
                    }
                })
                .partition::<VecDeque<_>, _>(|v| v.1.is_empty());
            let mut edge = frags.iter().map(|i| i.0).collect::<Vec<_>>();
            let mut scratch = Vec::new();
            frags.reserve(src.len());
            while !edge.is_empty() {
                for (i, subs) in &mut src {
                    if *i == usize::MAX {
                        continue;
                    }
                    {
                        if !subs.iter().any(|n| edge.contains(&n.index())) {
                            continue;
                        }
                    }
                    scratch.push(*i);
                    frags.push_back((*i, *subs));
                    *i = usize::MAX;
                }
                std::mem::swap(&mut edge, &mut scratch);
                if !edge.is_empty() {
                    scratch.clear();
                }
            }
            debug_assert!(
                src.iter().all(|i| i.0 == 0),
                "cycle detected in arena fragments?"
            );
            for frag in &mut frags {
                frag.0 -= 1;
            }
            frags
        };
        // let mut found = Vec::new();
        // let mut matched = Vec::new();
        while let Some((i, _)) = frags.pop_front() {
            let cmp = GraphCompactor::<Molecule<Ix, RefAcc<Ix>>>::new(self.molecule(Ix::new(i)));
            let mut found = false;
            for ism in crate::graph::algo::subgraph_isomorphisms_iter(
                &&cmp,
                &mol,
                &mut Atom::matches,
                &mut PartialEq::eq,
            ) {
                found = true;
            }
            if !found {
                frags.retain(|(_, subs)| subs.contains(&Ix::new(i)));
            }
        }
        todo!()
    }
}
