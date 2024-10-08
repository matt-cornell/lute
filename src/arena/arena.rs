//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use super::*;
use crate::graph::*;
use misc::DataValueMap;
use petgraph::graph::DefaultIx;
use petgraph::prelude::*;
use petgraph::visit::*;
use semisparse::SemiSparseGraph;
use slab::Slab;
use smallvec::SmallVec;
use std::cell::{Cell, UnsafeCell};
use std::collections::VecDeque;
use std::fmt::{self, from_fn, Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicUsize, Ordering};

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

/// Wrapper around a range of values with some convenience methods to mimic a bitset
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct Range<Ix>(pub Ix, pub Ix);
impl<Ix: IndexType> Range<Ix> {
    pub fn from_start_len(start: usize, len: usize) -> Self {
        Self(Ix::new(start), Ix::new(start + len))
    }
    pub fn start(&self) -> usize {
        self.0.index()
    }
    pub fn end(&self) -> usize {
        self.1.index()
    }
    pub fn nth(&self, val: usize) -> Option<usize> {
        let out = self.0.index() + val;
        (out < self.1.index()).then_some(out)
    }
    pub fn index(&self, val: usize) -> Option<usize> {
        (val < self.1.index()).then_some(())?;
        val.checked_sub(self.0.index())
    }
}
impl<Ix: Debug> Debug for Range<Ix> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}..{:?}", self.0, self.1)
    }
}

const ATOM_BIT_STORAGE: usize = 2;
type Graph<Ix> = SemiSparseGraph<Atom, Bond, Undirected, Ix>;
type BSType = crate::utils::bitset::BitSet<u16, ATOM_BIT_STORAGE>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(C)]
pub(crate) struct InterFragBond<Ix> {
    pub an: Ix,
    pub bn: Ix,
    pub ai: Ix,
    pub bi: Ix,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BrokenMol<Ix> {
    pub frags: SmallVec<Ix, 2>,
    pub bonds: SmallVec<InterFragBond<Ix>, 2>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ModdedMol<Ix> {
    pub base: Ix,
    pub patch: SmallVec<(Ix, Atom), 4>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum MolRepr<Ix: IndexType> {
    Atomic(Range<Ix>),
    Broken(BrokenMol<Ix>),
    Modify(ModdedMol<Ix>),
    TempEmpty,
    Empty,
}

/// Wrapper around the index of a molecule
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MolIndex<Ix>(pub Ix);
impl<Ix: IndexType> MolIndex<Ix> {
    pub fn index(self) -> usize {
        self.0.index()
    }
    pub fn in_arena<'a, 'b: 'a, R: ArenaAccessible<Ix = Ix>>(
        &self,
        arena: &'b R,
    ) -> Molecule<Ix, R::Access<'a>> {
        Molecule::from_arena(arena, *self)
    }
}
impl<Ix> From<Ix> for MolIndex<Ix> {
    fn from(value: Ix) -> Self {
        Self(value)
    }
}
impl<Ix: Display> Display for MolIndex<Ix> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Ix::fmt(&self.0, f)
    }
}

#[derive(Debug, Clone)]
pub struct Fragment<Ix: IndexType, D> {
    pub custom: D,
    pub seen: bool,
    pub(crate) repr: MolRepr<Ix>,
    pub(crate) size: Ix,
    refcnt: usize,
}
impl<Ix: IndexType, D> Fragment<Ix, D> {
    pub fn size(&self) -> usize {
        self.size.index()
    }
    pub fn tracked(&self) -> bool {
        self.seen || self.size() >= 3
    }
}
impl<Ix: IndexType, D: PartialEq> PartialEq for Fragment<Ix, D> {
    fn eq(&self, other: &Self) -> bool {
        self.size == other.size && self.repr == other.repr && self.custom == other.custom
    }
}

struct AMatch<'a, Ix: IndexType> {
    matched: &'a [Cell<(Ix, usize, bool, u8)>],
    exact: bool,
}
impl<G0, G1: GraphBase, Ix: IndexType> NodeMatcher<G0, &GraphCompactor<G1>> for AMatch<'_, Ix>
where
    G0: DataValueMap<NodeWeight = Atom>,
    G1: NodeIndexable,
    GraphCompactor<G1>: DataValueMap<NodeId = G1::NodeId, NodeWeight = Atom>,
{
    fn enabled() -> bool {
        true
    }
    fn eq(&mut self, g0: &G0, g1: &&GraphCompactor<G1>, n0: G0::NodeId, n1: G1::NodeId) -> bool {
        let Some(l) = g0.node_weight(n0) else {
            return false;
        };
        let Some(mut r) = g1.node_weight(n1) else {
            return false;
        };
        let mat = self.matched[g1.graph.to_index(n1)].get();
        r.single_to_unknown(mat.3)
            .expect("Too many unknown groups would exist on this atom!");
        if self.exact {
            l == r
        } else {
            l.matches(&r)
        }
    }
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Clone)]
pub struct Arena<Ix: IndexType = DefaultIx, D = ()> {
    graph: Graph<Ix>,
    pub(crate) frags: SmallVec<Fragment<Ix, D>, 8>,
    next_vacant_frag: Ix,
    ordered: bool,
    contained: bool,
    require_contained: bool,
}
impl<Ix: IndexType, D> Arena<Ix, D> {
    pub fn new() -> Self {
        Self {
            graph: Graph::default(),
            frags: SmallVec::new(),
            next_vacant_frag: IndexType::max(),
            contained: true,
            ordered: true,
            require_contained: true,
        }
    }

    #[inline(always)]
    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
    }

    /// Check if the fragments are ordered by containment, i.e. a fragment contained in another has a smaller index.
    #[inline(always)]
    pub fn is_ordered(&self) -> bool {
        self.ordered
    }

    /// Check if every (seen) molecule that's been looked up will be tracked as contained in a molecule if it's subgraph isomorphic.
    ///
    /// When false, `contained_group` may return false even when it should return true based on the structure.
    #[inline(always)]
    pub fn is_contained(&self) -> bool {
        self.contained
    }

    /// Get a reference to the fragments.
    #[inline(always)]
    pub fn expose_frags(&self) -> &[Fragment<Ix, D>] {
        &self.frags
    }

    #[inline]
    fn push_frag(&mut self, frag: Fragment<Ix, D>) -> Ix {
        let max = <Ix as IndexType>::max().index() - 1;
        if self.next_vacant_frag == IndexType::max() {
            let idx = self.frags.len();
            assert!(
                idx < max,
                "too many fragments in molecule: limit of {idx} reached!"
            );
            self.frags.push(frag);
            Ix::new(idx)
        } else {
            let ret = self.next_vacant_frag;
            let ni = ret.index();
            self.frags[ni] = frag;
            self.next_vacant_frag = self.frags[(ni + 1)..]
                .iter()
                .position(|f| f.repr == MolRepr::Empty)
                .map_or_else(IndexType::max, |n| {
                    let ret = n + ni;
                    if ret < max {
                        Ix::new(ret)
                    } else {
                        IndexType::max()
                    }
                });
            ret
        }
    }
    fn push_or_set_frag(&mut self, frag: Fragment<Ix, D>, to: Option<Ix>) -> Ix {
        if let Some(idx) = to {
            let old = &mut self.frags[idx.index()];
            old.size = frag.size;
            #[cfg(debug_assertions)]
            {
                let old = std::mem::replace(&mut old.repr, frag.repr);
                debug_assert_eq!(old, MolRepr::TempEmpty);
            }
            #[cfg(not(debug_assertions))]
            {
                old.repr = frag.repr;
            }
            idx
        } else {
            self.push_frag(frag)
        }
    }

    /// Check if `mol` contains `group`
    #[instrument(level = "trace", skip(self))]
    pub fn contains_group(&self, mol: MolIndex<Ix>, group: MolIndex<Ix>) -> bool {
        let mut stack = SmallVec::<_, 3>::new();
        let mut seen = BSType::new();
        self.contains_group_impl(mol.0, group.0, &mut stack, &mut seen)
    }

    fn contains_group_impl(
        &self,
        mol: Ix,
        group: Ix,
        stack: &mut SmallVec<Ix, 3>,
        seen: &mut BSType,
    ) -> bool {
        stack.clear();
        seen.clear();
        if mol.index() < group.index() && self.ordered {
            return false; // quirk of insertion order
        }
        if mol.index() >= self.frags.len() || group.index() >= self.frags.len() {
            return false;
        }
        stack.push(mol);
        while let Some(i) = stack.pop() {
            let idx = i.index();
            if i == group {
                return true;
            }
            if seen.get(idx) {
                continue;
            }
            seen.set(idx, true);
            match self.frags.get(idx).map(|s| &s.repr) {
                None => warn!(idx, "OOB node found when checking for membership"),
                Some(MolRepr::Modify(ModdedMol { base: i, .. })) => stack.push(*i),
                Some(MolRepr::Broken(BrokenMol { frags, .. })) => stack.extend_from_slice(frags),
                Some(MolRepr::Atomic(_) | MolRepr::TempEmpty | MolRepr::Empty) => {}
            }
        }
        false
    }

    /// Get a graph of the molecule at the given index. Note that `Molecule::from_arena` could give
    /// better results as it can borrow from `RefCell`s and `RwLock`s.
    pub fn molecule(&self, mol: MolIndex<Ix>) -> Molecule<Ix, access::RefAcc<Ix, D>> {
        Molecule::from_arena(self, mol)
    }

    /// Try to find a fragment in the arena. Unlike `insert_mol`, this is guaranteed to not (visibly) mutate the arena.
    ///
    /// Note that this can still break invariants in some cases since we kinda ignore fragments that aren't seen.
    pub fn lookup_mol<G>(&self, mol: G) -> Option<MolIndex<Ix>>
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + misc::DataValueMap
            + GraphProp<EdgeType = Undirected>
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
    {
        (0..self.frags.len()).find_map(|i| {
            let frag = &self.frags[i];
            (frag.seen && frag.size() == mol.node_count()).then_some(())?;
            let idx = MolIndex(Ix::new(i));
            is_isomorphic_matching(self.molecule(idx), mol, PartialEq::eq, PartialEq::eq, false)
                .then_some(idx)
        })
    }

    /// Assert that all invariants requested hold.
    ///
    /// This can be a very expensive operation, and should only be used for debugging.
    pub fn integrity_check(&self) {
        use crate::disp::smiles::generate_smiles_quick as gsq;
        for i in 0..self.frags.len() {
            for j in 0..self.frags.len() {
                if i == j {
                    continue;
                }
                let ii = Ix::new(i).into();
                let ji = Ix::new(j).into();
                let imol = self.molecule(ii);
                let jmol = self.molecule(ji);
                if is_isomorphic_matching(imol, jmol, PartialEq::eq, PartialEq::eq, false) {
                    panic!(
                        "fragments {i} and {j} are isomorphic (smiles: {:?})",
                        gsq(imol)
                    );
                }
                let sub_ism =
                    is_isomorphic_matching(jmol, imol, Atom::matches, PartialEq::eq, true);
                let contained = self.contains_group(ii, ji);
                if self.contained && self.frags[j].seen && sub_ism && !contained {
                    panic!("fragment {i} (smiles: {:?}) should contain fragment {j} (smiles: {:?}) but it does not", gsq(imol), gsq(jmol));
                }
                if self.ordered && contained && j > i {
                    panic!("fragment {i} (smiles: {:?}) should have a lower index than fragment {j} (smiles: {:?}) but it does not", gsq(imol), gsq(jmol));
                }
            }
        }
    }

    fn gc_from(&mut self, repr: &MolRepr<Ix>) {
        let mut stack = SmallVec::<_, 3>::new();
        match repr {
            MolRepr::Atomic(b) => self.graph.free_range(b.start(), b.end()),
            MolRepr::Broken(BrokenMol { frags, .. }) => stack.extend_from_slice(frags),
            MolRepr::Modify(ModdedMol { base, .. }) => stack.push(*base),
            _ => {}
        }
        while let Some(i) = stack.pop() {
            let frag = &mut self.frags[i.index()];
            frag.refcnt -= 1;
            if frag.seen || frag.refcnt > 0 || frag.repr == MolRepr::TempEmpty {
                continue;
            }
            if i.index() < self.next_vacant_frag.index() {
                self.next_vacant_frag = i;
            }
            match std::mem::replace(&mut frag.repr, MolRepr::Empty) {
                MolRepr::Atomic(b) => self.graph.free_range(b.start(), b.end()),
                MolRepr::Broken(BrokenMol { frags, .. }) => stack.extend_from_slice(&frags),
                MolRepr::Modify(ModdedMol { base, .. }) => stack.push(base),
                _ => {}
            }
        }
    }

    /// Cleanup unused fragments. Probably doesn't need to be called.
    pub fn gc(&mut self) {
        let mut stack = SmallVec::<_, 3>::new();
        for frag in &mut self.frags {
            if frag.seen || frag.refcnt > 0 {
                continue;
            }
            match std::mem::replace(&mut frag.repr, MolRepr::Empty) {
                MolRepr::Atomic(b) => self.graph.free_range(b.start(), b.end()),
                MolRepr::Broken(BrokenMol { frags, .. }) => stack.extend_from_slice(&frags),
                MolRepr::Modify(ModdedMol { base, .. }) => stack.push(base),
                _ => {}
            }
        }
        while let Some(i) = stack.pop() {
            let frag = &mut self.frags[i.index()];
            frag.refcnt -= 1;
            if frag.seen || frag.refcnt > 0 {
                continue;
            }
            match std::mem::replace(&mut frag.repr, MolRepr::Empty) {
                MolRepr::Atomic(b) => self.graph.free_range(b.start(), b.end()),
                MolRepr::Broken(BrokenMol { frags, .. }) => stack.extend_from_slice(&frags),
                MolRepr::Modify(ModdedMol { base, .. }) => stack.push(base),
                _ => {}
            }
        }
    }
}
impl<Ix: IndexType, D: Default> Arena<Ix, D> {
    /// Enforce that the arena is ordered.
    /// See `is_ordered` for what this means.
    pub fn make_ordered(&mut self) {
        if !self.ordered {
            self.optimize_layout(None);
        }
    }

    /// Enforce that all arena fragments are contained properly.
    /// See `is_contained` for what this means.
    pub fn make_contained(&mut self) {
        if !self.contained {
            self.optimize_layout(None);
        }
        self.require_contained = true;
    }

    fn fix_frags_containing(&mut self, idx: Ix) {
        let _guard = debug_span!("fixing fragments", fixing = idx.index()).entered();
        let mut queue = Vec::new();
        let mut graph = UnGraph::<Atom, Bond, Ix>::default();
        let mut base = 0;
        let node_count = self.frags[idx.index()].size();
        let frags = {
            let mut scratch = Vec::new();
            let (mut src, mut frags) = self
                .frags
                .iter()
                .enumerate()
                .filter_map(|(n, frag)| {
                    let children = match &frag.repr {
                        MolRepr::TempEmpty | MolRepr::Empty => None?,
                        MolRepr::Atomic(_) => &[] as &[_],
                        MolRepr::Broken(b) => &b.frags,
                        MolRepr::Modify(m) => std::slice::from_ref(&m.base),
                    };
                    (frag.size() >= node_count && n != idx.index())
                        .then_some((Ix::new(n), unsafe { &*(children as *const [Ix]) }))
                })
                .partition::<Vec<_>, _>(|v| {
                    v.1.iter().any(|n| {
                        let frag = &self.frags[n.index()];
                        frag.size() >= node_count && *n != idx
                    })
                });
            let mut edge = frags.iter().map(|i| i.0).collect::<Vec<_>>();
            frags.reserve(src.len());

            while !edge.is_empty() {
                for (i, subs) in &mut src {
                    if *i == IndexType::max() {
                        continue;
                    }
                    {
                        if !subs.iter().any(|n| edge.contains(n)) {
                            continue;
                        }
                    }
                    scratch.push(*i);
                    frags.push((*i, *subs));
                    *i = IndexType::max();
                }
                std::mem::swap(&mut edge, &mut scratch);
                if !edge.is_empty() {
                    scratch.clear();
                }
            }
            debug_assert!(
                src.iter().all(|i| i.0 == IndexType::max()),
                "cycle detected in arena fragments?"
            );
            trace!(?frags, "topo-sorted fragments");
            frags
        };
        for (frag, children) in frags {
            let mol = self.molecule(idx.into());
            let m2 = self.molecule(frag.into());
            if !(queue.iter().rev().any(|(q, _)| children.contains(q))
                || (is_isomorphic_matching(mol, m2, Atom::matches, PartialEq::eq, true)
                    && !self.contains_group(frag.into(), idx.into())))
            {
                trace!(frag = frag.index(), "pruning fragment");
                continue;
            }
            trace!(frag = frag.index(), "adding fragment to graph");
            for atom in m2.node_references() {
                graph.add_node(*atom.weight());
            }
            for edge in m2.edge_references() {
                graph.add_edge(
                    IndexType::new(edge.source().0.index() + base),
                    IndexType::new(edge.target().0.index() + base),
                    *edge.weight(),
                );
            }
            base += m2.node_count();
            queue.push((frag, base));
        }
        self.ordered &= !queue.is_empty();
        for frag in &queue {
            let old = std::mem::replace(&mut self.frags[frag.0.index()].repr, MolRepr::TempEmpty);
            self.gc_from(&old);
        }
        base = 0;
        for (f, end) in queue {
            debug!(frag = f.index(), "re-inserting fragment");
            let new_mol = RangeFiltered::new(&graph, base, end);
            self.insert_mol_impl(&new_mol, 0, None, Some(f), 0);
            base = end;
        }
    }

    /// Insert a molecule into the arena, deduplicating common parts.
    fn insert_mol_impl<G>(
        &mut self,
        mol: G,
        isms_from: usize,
        bits: Option<(&mut BSType, usize)>,
        to: Option<Ix>,
        depth: usize,
    ) -> (Ix, Vec<usize>)
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + misc::DataValueMap
            + GraphProp<EdgeType = Undirected>
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
        G::NodeId: Hash + Eq,
    {
        let _guard = (depth == 0).then(|| {
            info_span!(
                "inserting molecule",
                id = INSERT_ID.fetch_add(1, Ordering::Relaxed),
            )
            .entered()
        });

        let node_count = if let Some((_, count)) = &bits {
            *count
        } else {
            mol.node_count()
        };
        if node_count == 0 {
            debug!("handling special case for empty molecule");
            if let Some(f) = self.frags.iter().position(|f| f.size() == 0) {
                return (Ix::new(f), Vec::new());
            } else {
                let idx = self.push_or_set_frag(
                    Fragment {
                        custom: D::default(),
                        repr: MolRepr::TempEmpty,
                        size: Ix::new(0),
                        seen: false,
                        refcnt: 0,
                    },
                    to,
                );
                return (idx, Vec::new());
            }
        }

        debug!(
            depth,
            bits = ?bits.as_ref().map(|i| &i.0),
            smiles = event_enabled!(Level::TRACE, "entering impl")
                .then(|| {
                    use crate::disp::smiles::generate_smiles_quick as smiles;
                    if let Some((bits, _)) = &bits {
                        smiles(&BitFiltered::<_, _, ATOM_BIT_STORAGE, false>::new(
                            mol,
                            (**bits).clone(),
                        ))
                    } else {
                        smiles(mol)
                    }
                }),
            "entering impl"
        );
        struct OnExit(usize);
        impl Drop for OnExit {
            fn drop(&mut self) {
                debug!(depth = self.0, "exiting impl");
            }
        }
        let _guard = OnExit(depth);
        let max = <Ix as IndexType>::max().index();
        assert!(
            node_count < max,
            "molecule has too many atoms: {node_count}, max is {max}"
        );
        let mut scratch = Vec::new();
        let mut frags = {
            let (mut src, frags) = self
                .frags
                .iter()
                .enumerate()
                .filter_map(|(n, frag)| {
                    let children = match &frag.repr {
                        MolRepr::TempEmpty | MolRepr::Empty => None?,
                        MolRepr::Atomic(_) => &[] as &[_],
                        MolRepr::Broken(b) => &b.frags,
                        MolRepr::Modify(m) => std::slice::from_ref(&m.base),
                    };
                    (frag.size() == node_count
                        || (frag.size() < node_count && n >= isms_from && frag.tracked()))
                    .then_some((n, children))
                })
                .partition::<Vec<_>, _>(|v| {
                    v.1.iter().any(|n| {
                        let frag = &self.frags[n.index()];
                        frag.size() == node_count
                            || (frag.size() < node_count
                                && n.index() >= isms_from
                                && frag.tracked())
                    })
                });
            let mut frags = VecDeque::from(frags);
            let mut edge = frags.iter().map(|i| i.0).collect::<Vec<_>>();
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
                src.iter().all(|i| i.0 == usize::MAX),
                "cycle detected in arena fragments?\nsrc: {src:?}"
            );
            trace!(?frags, "topo-sorted fragments");
            frags
        };
        let default = (<Ix as IndexType>::max(), usize::MAX, false, 0u8);
        // vec of (index in frags, index in parts, perfect match)
        let mut matched = vec![default; mol.node_count()]; // this may be inefficient for matching small fragments
        if let Some((bits, _)) = &bits {
            for (n, (_, _, _, i)) in matched.iter_mut().enumerate() {
                if !bits.get(n) {
                    continue;
                }
                let ni = mol.from_index(n);
                *i = mol
                    .neighbors(ni)
                    .filter(|&n2| !bits.get(mol.to_index(n2)))
                    .count()
                    .try_into()
                    .expect("More than 255 pruned groups isn't possible!");
            }
        }
        let seen_buf = UnsafeCell::new(BSType::new());
        let pred_buf = UnsafeCell::new(SmallVec::new());
        let mut prune_buf = Vec::new();
        let mut found = Slab::new();
        while let Some((i_, _)) = frags.pop_front() {
            let i = Ix::new(i_);
            let matched = Cell::from_mut(&mut *matched).as_slice_of_cells();
            let filtered = NodeFilter::new(mol, |n| {
                let idx = mol.to_index(n);
                if let Some(b) = &bits {
                    if !b.0.get(idx) {
                        return false;
                    }
                }
                let v = matched[idx].get();
                if v.0 == IndexType::max() {
                    return true;
                }
                let f1 = &self.frags[i.index()];
                let f2 = &self.frags[v.0.index()];
                // SAFETY: the UnsafeCells' data doesn't escape the call
                f1.seen.cmp(&f2.seen).then(f1.size().cmp(&f2.size())) == std::cmp::Ordering::Greater
                    || unsafe {
                        i != v.0
                            && self.contains_group_impl(
                                i,
                                v.0,
                                &mut *pred_buf.get(),
                                &mut *seen_buf.get(),
                            )
                    }
            });
            trace!(
                i = i.index(),
                filtered = ?from_fn(|f| f.debug_list().entries(filtered.node_identifiers().map(|i| mol.to_index(i))).finish()),
                "checking for isomorphisms"
            );
            let compacted = GraphCompactor::<NodeFilter<G, _>>::new(filtered);
            let frag = self.molecule(i.into());
            let mut amatch = AMatch {
                matched,
                exact: i_ < isms_from.index() || !self.frags[i_].tracked(),
            };
            let mut bmatch = PartialEq::eq;
            let mut found_any = false;
            'isms: for mut ism in
                isomorphisms_iter(&frag, &&compacted, &mut amatch, &mut bmatch, true)
            {
                trace!(frag = i.index(), ?ism, "found fragment");
                prune_buf.clear();
                prune_buf.reserve(ism.len());
                for (n, to) in ism.iter_mut().enumerate() {
                    let to_id = compacted.node_map[*to];
                    let new = mol.to_index(to_id);
                    let mut mol_a = mol.node_weight(to_id).unwrap();
                    let mat = matched[new].get();
                    mol_a
                        .single_to_unknown(mat.3)
                        .expect("Too many unknown groups would exist on this atom!");
                    let frag_a = self.molecule(i.into()).get_atom(Ix::new(n)).unwrap();
                    let matches = if mat.2 {
                        if mol_a != frag_a {
                            continue 'isms;
                        }
                        true
                    } else {
                        mol_a == frag_a
                    };
                    *to = new;
                    prune_buf.push((Ix::new(new), matches));
                }
                let idx = found.insert((i, ism));
                scratch.clear();
                for &(p, m) in &prune_buf {
                    let c = &matched[p.index()];
                    let u = c.get().3;
                    let old_idx = c.replace((i, idx, m, u)).1;
                    let res = found.try_remove(old_idx).is_some();
                    if res {
                        scratch.push(old_idx);
                    }
                }
                for c in matched {
                    if scratch.contains(&c.get().1) {
                        c.set(default);
                    }
                }
                found_any = true;
            }
            if !found_any {
                prune_buf.clear();
                prune_buf.push((i, false));
                frags.retain(|&(n, p)| {
                    !prune_buf.iter().rev().any(|(i, _)| p.contains(i)) || {
                        prune_buf.push((Ix::new(n), false));
                        false
                    }
                });
            }
        }
        for (_, (i, ism)) in &found {
            let cmp = self.molecule((*i).into());
            for (from, &to) in ism.iter().enumerate() {
                let mi = mol.from_index(to);
                let new = mol
                    .neighbors(mi)
                    .filter(|&n| !ism.contains(&mol.to_index(n)))
                    .count()
                    .try_into()
                    .expect("More than 255 pruned groups isn't possible!");
                matched[to].3 = new;
                let mut mol_a = mol.node_weight(mi).unwrap();
                let _ = mol_a.single_to_unknown(new);
                let cmp_a = cmp.get_atom(Ix::new(from)).unwrap();
                matched[to].2 = mol_a == cmp_a;
            }
        }
        {
            let matched = Cell::from_mut(&mut *matched).as_slice_of_cells();
            for (n, c) in matched.iter().enumerate() {
                if c.get().1 != usize::MAX {
                    continue;
                }
                let pruned = mol
                    .neighbors(mol.from_index(n))
                    .filter(|&n2| {
                        let i2 = mol.to_index(n2);
                        matched[i2].get().1 != usize::MAX
                            || bits.as_ref().map_or(false, |b| !b.0.get(i2))
                    })
                    .count()
                    .try_into()
                    .expect("More than 255 pruned groups isn't possible!");
                c.set((IndexType::max(), usize::MAX, false, pruned))
            }
        }
        for (_, (i, ism)) in &mut found {
            let patch: SmallVec<_, 4> = ism
                .iter()
                .enumerate()
                .filter_map(|(from, &to)| {
                    let mat = matched[to];
                    (!mat.2).then(|| {
                        let mut atom = mol.node_weight(mol.from_index(to)).unwrap();
                        let _ = atom.single_to_unknown(mat.3);
                        (Ix::new(from), atom)
                    })
                })
                .collect();
            if patch.is_empty() {
                continue;
            }
            self.frags[i.index()].refcnt += 1;
            *i = self.push_frag(Fragment {
                custom: D::default(),
                repr: MolRepr::Modify(ModdedMol { base: *i, patch }),
                size: Ix::new(ism.len()),
                seen: ism.len() >= 3,
                refcnt: 0,
            });
        }
        if found.is_empty() {
            let mut mapping = vec![(IndexType::max(), 0); mol.node_count()];
            let mut count = 0;
            let mut it = mol.node_references();
            let start = self.graph.allocate_range(node_count, || {
                it.find_map(|n| {
                    let mi = mol.to_index(n.id());
                    if let Some(b) = &bits {
                        if !b.0.get(mi) {
                            return None;
                        }
                    }
                    let mut a = *n.weight();
                    let _ = a.single_to_unknown(matched[mi].3);
                    mapping[mi] = (count, mi);
                    count += 1;
                    Some(a)
                })
                .unwrap()
            });
            for e in mol.edge_references() {
                let a = mapping[mol.to_index(e.source())].0;
                let b = mapping[mol.to_index(e.target())].0;
                if a == IndexType::max() || b == IndexType::max() {
                    continue;
                }
                self.graph.add_edge(
                    NodeIndex::new(start + a),
                    NodeIndex::new(start + b),
                    *e.weight(),
                );
            }
            mapping.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            let new_map = mapping
                .into_iter()
                .filter_map(|x| (x.0 != IndexType::max()).then_some(x.1))
                .collect();
            let idx = self.push_or_set_frag(
                Fragment {
                    custom: D::default(),
                    repr: MolRepr::Atomic(Range::from_start_len(start, node_count)),
                    size: Ix::new(count),
                    seen: false,
                    refcnt: 0,
                },
                to,
            );
            debug!(idx = idx.index(), "no fragments found, inserting molecule");
            return (idx, new_map);
        } else if found.len() == 1 {
            let (n, (_, ism)) = found.iter().next().unwrap();
            if ism.len() == mol.node_count() {
                let res = found.remove(n);
                debug!(idx = res.0.index(), "found a single fragment");
                return res;
            }
        }
        let old_start = self.frags.len();
        let gen_mapping = bits.is_some();
        {
            let matched = Cell::from_mut(&mut *matched).as_slice_of_cells();
            let filtered = NodeFilter::new(mol, |n| {
                let n = mol.to_index(n);
                matched[n].get().0 == IndexType::max() && bits.as_ref().map_or(true, |b| b.0.get(n))
            });
            let mut cgi = ConnectedGraphIter::<u16, ATOM_BIT_STORAGE>::new(&filtered);
            let mut bits = BSType::new();
            let mut ism_buf = Vec::new();
            while let Some(count) = cgi.step(&filtered, &mut bits) {
                trace!(size = count.get(), "inserting fragment");
                let (i, ism) = self.insert_mol_impl(
                    mol,
                    old_start,
                    Some((&mut bits, count.get())),
                    None,
                    depth + 1,
                );
                ism_buf.clone_from(&ism);
                let idx = found.insert((i, ism));
                for ni in ism_buf.drain(..) {
                    let c = &matched[ni];
                    let u = c.get().3;
                    c.set((i, idx, true, u));
                }
            }
        }
        found.compact(|_, from, to| {
            for m in &mut matched {
                if m.1 == from {
                    m.1 = to;
                }
            }
            true
        });
        if found.len() == 1 {
            let frag = found.into_iter().next().unwrap().1;
            debug!(idx = frag.0.index(), "inserted a single fragment");
            if let Some(idx) = to {
                if idx == frag.0 {
                    return frag;
                }
                let [old, new] = self
                    .frags
                    .get_many_mut([frag.0.index(), idx.index()])
                    .unwrap();
                std::mem::swap(&mut old.repr, &mut new.repr);
                std::mem::swap(&mut old.size, &mut new.size);
                if frag.0.index() + 1 == self.frags.len() {
                    self.frags.pop();
                }
                return (idx, frag.1);
            } else {
                return frag;
            }
        }
        let mut bonds = SmallVec::with_capacity(found.len().saturating_sub(1));
        for e in mol.edge_references() {
            let mut ami = mol.to_index(e.source());
            let mut bmi = mol.to_index(e.target());
            if let Some((bits, _)) = &bits {
                if !bits.get(ami) || !bits.get(bmi) {
                    continue;
                }
            }
            let mut an = matched[ami].1;
            let mut bn = matched[bmi].1;
            if an == bn {
                continue;
            }
            if an > bn {
                std::mem::swap(&mut an, &mut bn);
                std::mem::swap(&mut ami, &mut bmi);
            }
            let ai = found[an].1.iter().position(|&to| to == ami).unwrap();
            let bi = found[bn].1.iter().position(|&to| to == bmi).unwrap();
            bonds.push(InterFragBond {
                an: Ix::new(an),
                bn: Ix::new(bn),
                ai: Ix::new(ai),
                bi: Ix::new(bi),
            });
        }
        let out_map = if gen_mapping {
            found.iter().flat_map(|i| &i.1 .1).copied().collect()
        } else {
            Vec::new()
        };
        let frags = found
            .into_iter()
            .map(|f| {
                let res = f.1 .0;
                self.frags[res.index()].refcnt += 1;
                res
            })
            .collect();
        let frag = Fragment {
            custom: D::default(),
            repr: MolRepr::Broken(BrokenMol { frags, bonds }),
            size: Ix::new(node_count),
            seen: false,
            refcnt: 0,
        };
        let idx = self.push_or_set_frag(frag, to);
        debug!(idx = idx.index(), "returning a new broken fragment");
        (idx, out_map)
    }

    pub fn insert_mol<G>(&mut self, mol: G) -> MolIndex<Ix>
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + misc::DataValueMap
            + GraphProp<EdgeType = Undirected>
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
        G::NodeId: Hash + Eq,
    {
        debug!(
            size = mol.node_count(),
            smiles = event_enabled!(Level::TRACE, "entering impl")
                .then(|| crate::disp::smiles::generate_smiles_quick(mol)),
            "top-level insertion"
        );
        let res = self.insert_mol_impl(mol, 0, None, None, 0).0;
        self.frags[res.index()].seen = true;
        if self.require_contained {
            self.fix_frags_containing(res);
        } else {
            self.contained = false;
        }
        res.into()
    }

    /// Optimize the layout, ensuring that all elements are inserted optimally.
    /// A permutation can be output through a mutable reference.
    #[instrument(skip_all, fields(perm = perm.is_some()))]
    pub fn optimize_layout(&mut self, mut perm: Option<&mut Vec<MolIndex<Ix>>>) {
        if self.frags.is_empty() {
            self.ordered = true;
            self.contained = true;
            return;
        }
        if let Some(perm) = &mut perm {
            perm.clear();
            perm.resize(self.frags.len(), <Ix as IndexType>::max().into());
        }
        let mut graph = UnGraph::with_capacity(self.frags.iter().map(|i| i.size()).sum(), 0);
        let mut graph_indices = Vec::with_capacity(self.frags.len());
        let mut base = 0;
        for i in 0..self.frags.len() {
            if self.frags[i].tracked() {
                trace!(i, "adding fragment to graph");
                let mol = self.molecule(Ix::new(i).into());
                mol.node_references().for_each(|n| {
                    graph.add_node(n.atom);
                });
                mol.edge_references().for_each(|e| {
                    graph.add_edge(
                        Ix::new(e.source().0.index() + base).into(),
                        Ix::new(e.target().0.index() + base).into(),
                        *e.weight(),
                    );
                });
                base += mol.node_count();
            }
            graph_indices.push(base);
        }
        debug!(
            nodes = graph.node_count(),
            edges = graph.edge_count(),
            frags = self.frags.len(),
            "built graph"
        );
        let mut frags = self
            .frags
            .iter()
            .enumerate()
            .map(|(n, s)| {
                let seen = s.tracked();
                (Ix::new(n), if seen { s.size() } else { 0 }, seen)
            })
            .collect::<Vec<_>>();
        frags.sort_by_key(|x| x.1);
        let mut start_ = Some(0);
        let mut scratch = Vec::new();
        let mut prune_buf = Vec::new();
        while let Some(start) = start_ {
            let base = frags[start].1;
            let len = frags[start..].iter().position(|x| x.1 > base);
            trace!(size = base, len, "arranging fragments with same size");
            let slice = if let Some(len) = len {
                &mut frags[start..(start + len)]
            } else {
                &mut frags[start..]
            };
            start_ = len.map(|l| start + l);
            if base == 0 {
                continue;
            }
            scratch.clear();
            scratch.extend(slice.iter().map(|s| (s.0, Vec::new(), s.2)));
            for i in 0..scratch.len() {
                for j in 0..scratch.len() {
                    if i == j {
                        continue;
                    }
                    let frag_0 = slice[i].0.index();
                    let frag_1 = slice[j].0.index();
                    if i > j && scratch[i].1.contains(&Ix::new(frag_1)) {
                        continue;
                    }
                    let range_0 = RangeFiltered::new(
                        &graph,
                        if frag_0 == 0 {
                            0
                        } else {
                            graph_indices[frag_0 - 1]
                        },
                        graph_indices[frag_0],
                    );
                    let range_1 = RangeFiltered::new(
                        &graph,
                        if frag_1 == 0 {
                            0
                        } else {
                            graph_indices[frag_1 - 1]
                        },
                        graph_indices[frag_1],
                    );
                    let res = is_isomorphic_matching(
                        &range_0,
                        &range_1,
                        Atom::matches,
                        PartialEq::eq,
                        true,
                    );
                    if res {
                        scratch[j].1.push(Ix::new(frag_0));
                    }
                }
            }
            let mut i = 0;
            while !scratch.is_empty() {
                prune_buf.clear();
                scratch.retain(|(n, c, s)| {
                    !c.is_empty() || {
                        slice[i] = (*n, 0, *s);
                        i += 1;
                        prune_buf.push(*n);
                        false
                    }
                });
                for (_, c, _) in &mut scratch {
                    c.retain(|n| !prune_buf.contains(n));
                }
            }
        }
        self.graph.clear();
        self.frags.clear();
        for &(n, _, s) in &frags {
            if !s {
                continue;
            }
            let frag = n.index();
            trace!(frag, "re-inserting fragment");
            let range = RangeFiltered::new(
                &graph,
                if frag == 0 {
                    0
                } else {
                    graph_indices[frag - 1]
                },
                graph_indices[frag],
            );
            let idx = self.insert_mol_impl(&range, 0, None, None, 0).0;
            debug!(old = frag, new = idx.index(), "re-inserting fragment");
            if let Some(perm) = &mut perm {
                perm[frag] = idx.into();
            }
        }
        self.ordered = true;
        self.contained = true;
    }
}
impl<Ix: IndexType, D> Default for Arena<Ix, D> {
    fn default() -> Self {
        Self::new()
    }
}

static INSERT_ID: AtomicUsize = AtomicUsize::new(0);
