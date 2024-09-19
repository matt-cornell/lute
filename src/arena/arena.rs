//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use super::*;
use crate::graph::*;
use petgraph::graph::DefaultIx;
use petgraph::prelude::*;
use petgraph::visit::*;
use slab::Slab;
use smallvec::SmallVec;
use std::cell::{Cell, UnsafeCell};
use std::collections::VecDeque;
use std::fmt::{self, Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, AtomicU8, Ordering};

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
type BSType = crate::utils::bitset::BitSet<u16, ATOM_BIT_STORAGE>;

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ModdedMol<Ix> {
    pub base: Ix,
    pub patch: SmallVec<(Ix, Atom), 4>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MolRepr<Ix: IndexType> {
    Atomic(BSType),
    Broken(BrokenMol<Ix>),
    Modify(ModdedMol<Ix>),
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
    pub fn in_mut_arena<'a, 'b: 'a, R: ArenaAccessibleMut<Ix = Ix>>(
        &self,
        arena: &'b R,
    ) -> Molecule<Ix, R::AccessMut<'a>> {
        Molecule::from_mut_arena(arena, *self)
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

#[derive(Debug)]
pub struct Fragment<Ix: IndexType, D> {
    pub custom: D,
    pub(crate) repr: MolRepr<Ix>,
    pub(crate) size: Ix,
    pub(crate) seen: AtomicBool,
}
impl<Ix: IndexType, D> Fragment<Ix, D> {
    pub fn size(&self) -> usize {
        self.size.index()
    }
    pub fn seen(&self) -> bool {
        self.seen.load(Ordering::Relaxed)
    }
    pub fn set_seen(&self, seen: bool) -> bool {
        self.seen.swap(seen, Ordering::Relaxed)
    }
}
impl<Ix: IndexType, D: Clone> Clone for Fragment<Ix, D> {
    fn clone(&self) -> Self {
        Self {
            custom: self.custom.clone(),
            repr: self.repr.clone(),
            size: self.size,
            seen: AtomicBool::new(false),
        }
    }
}
impl<Ix: IndexType, D: PartialEq> PartialEq for Fragment<Ix, D> {
    fn eq(&self, other: &Self) -> bool {
        self.size == other.size && self.repr == other.repr && self.custom == other.custom
    }
}

/// Invariant state of the arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum InvariantState {
    /// Core invariants hold:
    /// ```ignore
    /// let ix1 = arena.insert_mol(mol);
    /// assert!(is_isomoprhic_matching(ix1.in_arena(arena), &mol, PartialEq::eq, PartialEq::eq, false)); // isomorphism
    /// let ix2 = arena.insert_mol(mol);
    /// assert_eq!(ix1, ix2); // inempotency
    /// ```
    CoreInvariants,
    /// Core invariants, plus:
    /// ```ignore
    /// let small_frag = arena.expose_frags()[i.index()];
    /// if small_frag.seen() || small_frag.size() >= 3 {
    ///     let small_mol = arena.molecule(i);
    ///     let big_mol = arena.molecule(j);
    ///     if is_isomoprhic_matching(small_mol, big_mol, Atom::matches, PartialEq::eq, true) {
    ///         assert!(arena.contains_group(j, i)); // full group structure tracked
    ///     }
    /// }
    /// ```
    AllFragments,
    /// Core invariants, plus:
    /// ```ignore
    /// if arena.contains(i, j) {
    ///     assert!(i >= j);
    /// }
    /// ```
    OrderedFragments,
    /// `AllFragments` + `OrderedFragments`
    AllOrdered,
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug)]
pub struct Arena<Ix: IndexType = DefaultIx, D = ()> {
    graph: Graph<Ix>,
    pub(crate) frags: SmallVec<Fragment<Ix, D>, 16>,
    invariants: AtomicU8,
}
impl<Ix: IndexType, D> Arena<Ix, D> {
    pub fn new() -> Self {
        Self {
            graph: Graph::default(),
            frags: SmallVec::new(),
            invariants: AtomicU8::new(InvariantState::AllOrdered as _),
        }
    }

    #[inline(always)]
    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
    }

    pub fn current_invariant(&self) -> InvariantState {
        match self.invariants.load(Ordering::Relaxed) {
            0 => InvariantState::CoreInvariants,
            1 => InvariantState::AllFragments,
            2 => InvariantState::OrderedFragments,
            3 => InvariantState::AllOrdered,
            _ => unreachable!(),
        }
    }

    /// Get a reference to the fragments. For debugging purposes only.
    #[inline(always)]
    pub fn expose_frags(&self) -> &[Fragment<Ix, D>] {
        &self.frags
    }

    #[inline]
    fn push_frag(&mut self, frag: Fragment<Ix, D>) -> Ix {
        let max = <Ix as IndexType>::max().index() - 1;
        let idx = self.frags.len();
        assert!(
            idx < max,
            "too many fragments in molecule: limit of {idx} reached!"
        );
        self.frags.push(frag);
        Ix::new(idx)
    }
    fn push_or_set_frag(&mut self, frag: Fragment<Ix, D>, to: Option<Ix>) -> Ix {
        if let Some(idx) = to {
            #[cfg(debug_assertions)]
            {
                let old = std::mem::replace(&mut self.frags[idx.index()], frag);
                debug_assert_eq!(old.repr, MolRepr::Empty);
            }
            #[cfg(not(debug_assertions))]
            {
                self.frags[idx.index()] = frag;
            }
            idx
        } else {
            self.push_frag(frag)
        }
    }

    #[instrument(level = "trace", skip(self))]
    pub fn contains_group(&self, mol: MolIndex<Ix>, group: MolIndex<Ix>) -> bool {
        let mut stack = SmallVec::<_, 3>::new();
        let mut seen = BSType::new();
        self.contains_group_impl(mol.0, group.0, &mut stack, &mut seen)
    }

    /// Check if `mol` contains `group`
    fn contains_group_impl(
        &self,
        mol: Ix,
        group: Ix,
        stack: &mut SmallVec<Ix, 3>,
        seen: &mut BSType,
    ) -> bool {
        stack.clear();
        seen.clear();
        if mol.index() < group.index()
            && self.invariants.load(Ordering::Relaxed) >= InvariantState::AllOrdered as _
        {
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
                Some(MolRepr::Atomic(_) | MolRepr::Empty) => {}
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
            (self.frags[i].size() == mol.node_count()).then_some(())?;
            let idx = MolIndex(Ix::new(i));
            is_isomorphic_matching(self.molecule(idx), mol, PartialEq::eq, PartialEq::eq).then(
                || {
                    self.frags[i].seen.store(true, Ordering::Relaxed);
                    idx
                },
            )
        })
    }
}
impl<Ix: IndexType, D: Default> Arena<Ix, D> {
    pub fn request_invariant(&mut self, state: InvariantState) {
        let inv = self.invariants.get_mut();
        if *inv >= state as _ {
            *inv = state as _;
        } else {
            self.optimize_layout(None);
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
                        repr: MolRepr::Empty,
                        size: Ix::new(0),
                        seen: AtomicBool::new(false),
                    },
                    to,
                );
                return (idx, Vec::new());
            }
        }
        trace!(?bits, depth, "entering impl");
        struct OnExit(usize);
        impl Drop for OnExit {
            fn drop(&mut self) {
                trace!(depth = self.0, "exiting impl");
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
            let (frags, mut src) = self.frags[isms_from..]
                .iter()
                .enumerate()
                .filter_map(|(n, frag)| {
                    let children = match &frag.repr {
                        MolRepr::Atomic(_) | MolRepr::Empty => &[] as &[_],
                        MolRepr::Broken(b) => &b.frags,
                        MolRepr::Modify(m) => std::slice::from_ref(&m.base),
                    };
                    (frag.size() <= node_count && (frag.size() >= 3 || frag.seen()))
                        .then_some((n + isms_from, children))
                })
                .partition::<Vec<_>, _>(|v| v.1.is_empty());
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
                "cycle detected in arena fragments?"
            );
            debug!(?frags, "topo-sorted fragments");
            frags
        };
        // vec of (index in frags, index in parts, perfect match)
        let mut matched = vec![(IndexType::max(), usize::MAX, false, 0u8); mol.node_count()]; // this may be inefficient for matching small fragments
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
        while let Some((i, _)) = frags.pop_front() {
            trace!(i, "checking for isomorphisms");
            let i = Ix::new(i);
            let matched = Cell::from_mut(&mut *matched).as_slice_of_cells();
            let filtered = NodeFilter::new(mol, |n| {
                let idx = mol.to_index(n);
                if let Some(b) = &bits {
                    if !b.0.get(idx) {
                        return false;
                    }
                }
                let v = matched[idx].get();
                // SAFETY: the UnsafeCells' data doesn't escape the call
                unsafe {
                    v.0 == IndexType::max()
                        || (i != v.0
                            && self.contains_group_impl(
                                i,
                                v.0,
                                &mut *pred_buf.get(),
                                &mut *seen_buf.get(),
                            ))
                }
            });
            let compacted = GraphCompactor::<NodeFilter<G, _>>::new(filtered);
            let frag = self.molecule(i.into());
            let mut amatch = Atom::matches;
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
                for &(p, m) in &prune_buf {
                    let c = &matched[p.index()];
                    let u = c.get().3;
                    let old_idx = c.replace((i, idx, m, u)).1;
                    found.try_remove(old_idx);
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
            *i = self.push_frag(Fragment {
                custom: D::default(),
                repr: MolRepr::Modify(ModdedMol { base: *i, patch }),
                size: Ix::new(ism.len()),
                seen: AtomicBool::new(false),
            });
        }
        if found.is_empty() {
            debug!("no fragments found, inserting molecule");
            let mut mapping = vec![(IndexType::max(), 0); mol.node_count()];
            let mut count = 0;
            for n in mol.node_references() {
                let mi = mol.to_index(n.id());
                if let Some(b) = &bits {
                    if !b.0.get(mi) {
                        continue;
                    }
                }
                let mut a = *n.weight();
                let _ = a.single_to_unknown(matched[mi].3);
                let gi = self.graph.add_node(a);
                mapping[mi] = (gi, mi);
                count += 1;
            }
            for e in mol.edge_references() {
                let a = mapping[mol.to_index(e.source())].0;
                let b = mapping[mol.to_index(e.target())].0;
                if a == IndexType::max() || b == IndexType::max() {
                    continue;
                }
                self.graph.add_edge(a, b, *e.weight());
            }
            mapping.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            let (new_bits, new_map) = mapping
                .into_iter()
                .filter_map(|x| (x.0 != IndexType::max()).then_some((x.0.index(), x.1)))
                .unzip();
            let idx = self.push_or_set_frag(
                Fragment {
                    custom: D::default(),
                    repr: MolRepr::Atomic(new_bits),
                    size: Ix::new(count),
                    seen: AtomicBool::new(false),
                },
                to,
            );
            return (idx, new_map);
        } else if found.len() == 1 {
            let (n, (_, ism)) = found.iter().next().unwrap();
            if ism.len() == mol.node_count() {
                return found.remove(n);
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
            return found.into_iter().next().unwrap().1;
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
        let frags = found.into_iter().map(|f| f.1 .0).collect();
        let frag = Fragment {
            custom: D::default(),
            repr: MolRepr::Broken(BrokenMol { frags, bonds }),
            size: Ix::new(node_count),
            seen: AtomicBool::new(false),
        };
        let idx = self.push_or_set_frag(frag, to);
        (idx, out_map)
    }
    #[inline(always)]
    #[instrument(skip(self, mol), fields(size = mol.node_count()))]
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
        let res = self.insert_mol_impl(mol, 0, None, None, 0).0;
        self.frags[res.index()].seen.store(true, Ordering::Relaxed);
        res.into()
    }

    /// Optimize the layout, ensuring that all elements are inserted optimally.
    /// A permutation can be output through a mutable reference.
    #[instrument(skip_all, fields(perm = perm.is_some()))]
    pub fn optimize_layout(&mut self, mut perm: Option<&mut Vec<MolIndex<Ix>>>) {
        if self.frags.is_empty() {
            *self.invariants.get_mut() = InvariantState::AllOrdered as _;
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
            if self.frags[i].seen() {
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
                let seen = s.seen();
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
                    let res =
                        is_isomorphic_matching(&range_0, &range_1, Atom::matches, PartialEq::eq);
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
        *self.invariants.get_mut() = InvariantState::AllOrdered as _;
    }
}
impl<Ix: IndexType, D> Default for Arena<Ix, D> {
    fn default() -> Self {
        Self::new()
    }
}
impl<Ix: IndexType, D: Clone> Clone for Arena<Ix, D> {
    fn clone(&self) -> Self {
        Self {
            graph: self.graph.clone(),
            frags: self.frags.clone(),
            invariants: AtomicU8::new(self.invariants.load(Ordering::Relaxed)),
        }
    }
}
