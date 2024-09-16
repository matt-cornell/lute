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

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MolRepr<Ix: IndexType> {
    Atomic(BSType),
    Broken(BrokenMol<Ix>),
    Modify(ModdedMol<Ix>),
}

/// The `Arena` is the backing storage for everything. It tracks all molecules and handles
/// deduplication.
#[derive(Debug, Default, Clone)]
pub struct Arena<Ix: IndexType = DefaultIx> {
    graph: Graph<Ix>,
    pub(crate) parts: SmallVec<(MolRepr<Ix>, Ix), 16>,
}
impl<Ix: IndexType> Arena<Ix> {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
    }

    /// Get a reference to the `parts` field. For debugging purposes only.
    #[inline(always)]
    pub fn expose_parts(&self) -> &(impl Debug + Clone + PartialEq) {
        &self.parts
    }
    /// Get a reference to a single fragment. For debugging purposes only.
    #[inline(always)]
    pub fn expose_part(&self, index: usize) -> &(impl Debug + Clone + PartialEq) {
        &self.parts[index]
    }

    #[inline]
    fn push_frag(&mut self, frag: (MolRepr<Ix>, Ix)) -> Ix {
        let max = <Ix as IndexType>::max().index() - 1;
        let idx = self.parts.len();
        assert!(
            idx < max,
            "too many fragments in molecule: limit of {idx} reached!"
        );
        self.parts.push(frag);
        Ix::new(idx)
    }

    pub fn contains_group(&self, mol: Ix, group: Ix) -> bool {
        let mut stack = SmallVec::<_, 3>::new();
        let mut seen = BSType::new();
        self.contains_group_impl(mol, group, &mut stack, &mut seen)
    }

    /// Check if `mol` contains `group`
    #[instrument(level = "trace", name = "contains_group", skip(self))]
    fn contains_group_impl(
        &self,
        mol: Ix,
        group: Ix,
        stack: &mut SmallVec<Ix, 3>,
        seen: &mut BSType,
    ) -> bool {
        stack.clear();
        seen.clear();
        // if mol.index() < group.index() {
        //     return false; // quirk of insertion order
        // }
        if mol.index() >= self.parts.len() || group.index() >= self.parts.len() {
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
            match self.parts.get(idx).map(|s| &s.0) {
                None => warn!(idx, "OOB node found when checking for membership"),
                Some(MolRepr::Modify(ModdedMol { base: i, .. })) => stack.push(*i),
                Some(MolRepr::Broken(BrokenMol { frags, .. })) => stack.extend_from_slice(frags),
                Some(MolRepr::Atomic(_)) => {}
            }
        }
        false
    }

    /// Get a graph of the molecule at the given index. Note that `Molecule::from_arena` could give
    /// better results as it can borrow from `RefCell`s and `RwLock`s.
    pub fn molecule(&self, mol: Ix) -> Molecule<Ix, access::RefAcc<Ix>> {
        Molecule::from_arena(self, mol)
    }

    /// Insert a molecule into the arena, deduplicating common parts.
    fn insert_mol_impl<G>(
        &mut self,
        mol: G,
        isms_from: usize,
        bits: Option<(&mut BSType, usize)>,
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
            let (frags, mut src) = self.parts[isms_from..]
                .iter()
                .enumerate()
                .filter_map(|(n, frag)| {
                    let children = match &frag.0 {
                        MolRepr::Atomic(_) => &[] as &[_],
                        MolRepr::Broken(b) => &b.frags,
                        MolRepr::Modify(m) => std::slice::from_ref(&m.base),
                    };
                    (frag.1.index() <= node_count).then_some((n + isms_from, children))
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
                    self.contains_group_impl(i, v.0, &mut *pred_buf.get(), &mut *seen_buf.get())
                }
            });
            let compacted = GraphCompactor::<NodeFilter<G, _>>::new(filtered);
            let frag = self.molecule(i);
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
                    let frag_a = self.molecule(i).get_atom(Ix::new(n)).unwrap();
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
            let cmp = self.molecule(*i);
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
            *i = self.push_frag((
                MolRepr::Modify(ModdedMol { base: *i, patch }),
                Ix::new(ism.len()),
            ));
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
            let idx = self.push_frag((MolRepr::Atomic(new_bits), Ix::new(count)));
            return (idx, new_map);
        } else if found.len() == 1 {
            let (n, (_, ism)) = found.iter().next().unwrap();
            if ism.len() == mol.node_count() {
                return found.remove(n);
            }
        }
        let old_start = self.parts.len();
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
                let (i, ism) =
                    self.insert_mol_impl(mol, old_start, Some((&mut bits, count.get())), depth + 1);
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
        let idx = self.push_frag((
            MolRepr::Broken(BrokenMol { frags, bonds }),
            Ix::new(node_count),
        ));
        (idx, out_map)
    }
    #[inline(always)]
    #[instrument(skip(self, mol), fields(size = mol.node_count()))]
    pub fn insert_mol<G>(&mut self, mol: G) -> Ix
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
        self.insert_mol_impl(mol, 0, None, 0).0
    }
}
