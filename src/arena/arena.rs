//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use super::*;
use crate::graph::*;
use hybridmap::HybridMap;
use itertools::Itertools;
use petgraph::graph::DefaultIx;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallvec::SmallVec;
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
const MODDED_GRAPH_LEN: usize = 4;
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

    /// Simpler, faster version of `insert_mol` for when we know there are no subgraphs
    /// Returns index and mapping where forall `i`:
    /// `self.molecule(ix).get_atom(i) == mol.get_node(mapping[i])`
    ///
    /// Checks for isomorpisms iff `find_isms`
    #[instrument(skip(self, mol), fields(size = mol.node_count()))]
    fn insert_mol_atomic<G>(&mut self, mol: G, find_isms: bool) -> (Ix, Vec<usize>)
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
        if find_isms {
            let frags = self.parts.iter().positions(|p| {
                p.1.index() == mol.node_count() && matches!(p.0, MolRepr::Atomic(_))
            });
            // avoid an allocation if we can
            let frags = if span_enabled!(Level::DEBUG, frags, "found available fragments") {
                let frags = frags.collect::<SmallVec<_, 4>>();
                debug!(?frags, "found available fragments");
                either::Right(frags.into_iter())
            } else {
                either::Left(frags)
            };
            let mut mods = SmallVec::<(Ix, Atom), 4>::with_capacity(mol.node_count());
            let mut node_map = vec![usize::MAX; mol.node_count()];
            let mut amatch = Atom::matches;
            let mut bmatch = PartialEq::eq;
            let mut frag = None;

            'main: for i in frags {
                debug!(idx = i, "matching fragment");
                let cmp = self.molecule(Ix::new(i));
                let mut it =
                    isomorphisms_iter(&cmp, &mol, &mut amatch, &mut bmatch, false).peekable();
                while let Some(ism) = it.next() {
                    debug_assert_eq!(ism.len(), mol.node_count());
                    mods.clear();

                    for (cmp_i, &mol_i) in ism.iter().enumerate() {
                        let graph_atom = cmp.get_atom(Ix::new(cmp_i)).unwrap();
                        let mol_atom = mol.node_weight(mol.from_index(mol_i)).unwrap();
                        if graph_atom != mol_atom {
                            let mi = Ix::new(mol_i);
                            if let Err(idx) = mods.binary_search_by_key(&mi, |m| m.0) {
                                mods.insert(idx, (mi, mol_atom));
                            }
                        }
                        node_map[mol_i] = cmp_i;
                    }

                    // perfect match!
                    if mods.is_empty() {
                        info!(idx = i, ?node_map, "perfect isomorphism");
                        return (Ix::new(i), node_map);
                    }

                    // not quite perfect, let's see if there's already another modded mol
                    for (idx, frag) in self.parts.iter().enumerate() {
                        if let MolRepr::Modify(m) = &frag.0 {
                            if m.base.index() == i && m.patch == mods {
                                info!(idx, base = i, ?node_map, "matched isomorphism");
                                return (Ix::new(idx), node_map);
                            }
                        }
                    }

                    // no similar moddeds exist, this is the last ism
                    if it.peek().is_none() {
                        info!(
                            idx = self.parts.len(),
                            base = i,
                            ?node_map,
                            "created isomorphism"
                        );
                        frag = Some((
                            (
                                MolRepr::Modify(ModdedMol {
                                    base: Ix::new(i),
                                    patch: mods,
                                }),
                                Ix::new(mol.node_count()),
                            ),
                            node_map,
                        ));
                        break 'main;
                    }
                }
            }
            if let Some((frag, map)) = frag {
                return (self.push_frag(frag), map);
            }
        }
        let end = petgraph::graph::NodeIndex::<Ix>::end();
        let mut bits = BSType::new();
        let mut node_map = vec![end; mol.node_bound()];
        let mut atom_map = Vec::with_capacity(mol.node_bound());
        for aref in mol.node_references() {
            let b = self.graph.add_node(*aref.weight());
            let i = mol.to_index(aref.id());
            node_map[i] = b;
            atom_map.push(i);
            bits.set(b.index(), true);
        }
        for eref in mol.edge_references() {
            let s = node_map[mol.to_index(eref.source())];
            let t = node_map[mol.to_index(eref.target())];
            debug_assert_ne!(s, end);
            debug_assert_ne!(t, end);
            self.graph.add_edge(s, t, *eref.weight());
        }
        let idx = self.push_frag((MolRepr::Atomic(bits), Ix::new(mol.node_count())));
        info!(idx = idx.index(), "inserting atomic");
        (idx, atom_map)
    }

    /// Insert a molecule into the arena, deduplicating common parts.
    #[instrument(skip(self, mol), fields(size = mol.node_count()))]
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
                .partition::<Vec<_>, _>(|v| v.1.is_empty());
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
                    frags.push((*i, *subs));
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
            debug!(?frags, "topo-sorted fragments");
            frags
        };

        let mut found = Vec::new();
        let mut matched = vec![(usize::MAX, 0); mol.node_bound()];
        let mut mods = SmallVec::<(Ix, Atom), 4>::with_capacity(mol.node_count());
        let mut additional_rs = SmallVec::<(G::NodeId, Atom, Bond), MODDED_GRAPH_LEN>::new();
        let mut rbonds = HybridMap::<G::NodeId, Atom, MODDED_GRAPH_LEN>::new();
        let mut idx = 0;
        let mut amatch = Atom::matches;
        let mut bmatch = PartialEq::eq;

        let mut search_stack = SmallVec::<_, 2>::new();
        let mut preds_found = SmallVec::<_, 3>::new();
        let mut push_list = SmallVec::<_, 8>::new();
        let mut frag = None;
        'main: while idx < frags.len() {
            let (i, subs) = frags[idx];
            let cmp = self.molecule(Ix::new(i));
            if span_enabled!(Level::TRACE, atoms) {
                trace!(atoms = ?cmp.node_references().map(|n| n.weight().to_string()).collect::<Vec<_>>());
            }
            let mut found_any = false;
            preds_found.clear();

            let mut it = isomorphisms_iter(&cmp, &mol, &mut amatch, &mut bmatch, true).peekable();
            trace!("before first ism check");
            if it.peek().is_some() {
                debug!(idx, "isomorphisms exist");
            }
            'isms: while let Some(ism) = it.next() {
                if ism.len() == mol.node_count() {
                    mods.clear();

                    for (cmp_i, &mol_i) in ism.iter().enumerate() {
                        let graph_atom = cmp.get_atom(Ix::new(cmp_i)).unwrap();
                        let mol_atom = mol.node_weight(mol.from_index(mol_i)).unwrap();
                        if graph_atom != mol_atom {
                            let mi = Ix::new(mol_i);
                            if let Err(idx) = mods.binary_search_by_key(&mi, |m| m.0) {
                                mods.insert(idx, (mi, mol_atom));
                            }
                        }
                    }

                    // perfect match!
                    if mods.is_empty() {
                        info!(idx = i, "perfect isomorphism");
                        return Ix::new(i);
                    }

                    // not quite perfect, let's see if there's already another modded mol
                    for (idx, frag) in self.parts.iter().enumerate() {
                        if let MolRepr::Modify(m) = &frag.0 {
                            if m.base.index() == i && m.patch == mods {
                                info!(idx, base = i.index(), "matched isomorphism");
                                return Ix::new(idx);
                            }
                        }
                    }

                    // no similar moddeds exist, this is the last ism
                    if it.peek().is_none() {
                        info!(idx = self.parts.len(), base = i, "created isomorphism");
                        frag = Some((
                            MolRepr::Modify(ModdedMol {
                                base: Ix::new(i),
                                patch: mods,
                            }),
                            Ix::new(mol.node_count()),
                        ));
                        break 'main;
                    }
                }
                debug!(idx = i, ?ism, "found subgraph isomorphism");
                found_any = true;
                push_list.clear();
                for (cmp_i, &mol_i) in ism.iter().enumerate() {
                    let graph_id = Ix::new(cmp_i).into();
                    let mol_id = mol.from_index(mol_i);
                    let graph_atom = cmp.get_atom(graph_id).unwrap();
                    let mol_atom = mol.node_weight(mol_id).unwrap();
                    trace!(graph_atom.protons, mol_atom.protons);

                    // find the unmatched neighbors
                    let mut neighbors = mol
                        .edges(mol_id)
                        .map(|e| {
                            let i = if e.source() == mol_id {
                                e.target()
                            } else {
                                e.source()
                            };
                            (*e.weight(), i)
                        })
                        .collect::<SmallVec<_, 4>>();
                    cmp.neighbors(graph_id)
                        .for_each(|n| neighbors.retain(|e| ism[mol.to_index(e.1)] != n.0.index()));

                    // graph atom is an R-group, matched on mol. no need to track it as being owned
                    if graph_atom.protons == 0 && !mol_atom.protons == 0 {
                        additional_rs.reserve(neighbors.len());
                        for (b, n) in neighbors {
                            let idx = additional_rs
                                .binary_search_by_key(&mol.to_index(n), |r| mol.to_index(r.0))
                                .unwrap_or_else(|x| x);
                            trace!(neighbor = mol.to_index(n), bond = ?b, "adding additional R-group");
                            additional_rs.insert(idx, (n, Atom::new(0), b));
                        }
                        continue;
                    }

                    if matched[mol_i].0 != usize::MAX {
                        trace!(
                            mol_i,
                            frag = matched[mol_i].0,
                            cmp_i = matched[mol_i].1,
                            "atom already matched"
                        );
                        // search through subgraphs
                        if preds_found.is_empty() && !subs.is_empty() {
                            let _span =
                                trace_span!("searching predecessors", frag = mol_i).entered();
                            search_stack.clear();
                            search_stack.extend_from_slice(subs);
                            while let Some(pred) = search_stack.pop() {
                                trace!(pred = pred.index(), "found predecessor");
                                if !preds_found.contains(&pred) {
                                    preds_found.push(pred);
                                }
                                match self.parts[pred.index()].0 {
                                    MolRepr::Broken(ref b) => {
                                        search_stack.extend_from_slice(&b.frags)
                                    }
                                    MolRepr::Atomic(_) => {}
                                    MolRepr::Redirect(to)
                                    | MolRepr::Modify(ModdedMol { base: to, .. }) => {
                                        search_stack.push(to)
                                    }
                                }
                            }
                            // sort it just for that lookup speed boost
                            preds_found.sort();
                        }
                        // predecessor not found, this atom is already accounted for
                        if preds_found.binary_search(&Ix::new(mol_i)).is_err() {
                            trace!("predecessor not found");
                            continue 'isms;
                        }
                    }
                    let diff = mol_atom.data.other() - cmp.neighbors(graph_id).count() as u8;
                    assert!(
                        neighbors.len() < 2,
                        "too many neighbors: {}",
                        neighbors.len()
                    );
                    for (b, n) in neighbors {
                        assert_eq!(b, Bond::Single);
                        debug!(
                            neighbor = mol.to_index(n),
                            "adding an additional suppressed R-group"
                        );
                        if let Some(a) = rbonds.get_mut(&n) {
                            a.data.set_other(a.data.other() - diff);
                            a.data.set_unknown(a.data.unknown() + diff);
                        } else {
                            let mut a = mol.node_weight(n).unwrap();
                            a.data.set_other(a.data.other() - diff);
                            a.data.set_unknown(a.data.unknown() + diff);
                            rbonds.insert(n, a);
                        }
                    }
                    push_list.push((mol_i, i, cmp_i));
                }
                for &(mol_i, i, cmp_i) in &push_list {
                    matched[mol_i] = (i, cmp_i);
                }
                found.push((i, ism));
            }
            if found_any {
                idx += 1;
            } else {
                let mut index = 0;
                frags.retain(|(_, subs)| {
                    let res = index < idx || (index != idx && !subs.contains(&Ix::new(i)));
                    index += 1;
                    res
                });
                if idx + 1 == frags.len() {
                    break;
                }
            }
        }
        if let Some(frag) = frag {
            return self.push_frag(frag);
        }
        if found.is_empty() {
            info!("no subgraphs found, delegating to atomic");
            return self.insert_mol_atomic(mol, false).0;
        }
        let modded = ModdedGraph {
            graph: mol,
            mods: rbonds,
            additional: additional_rs,
        };
        trace!(
            n_mods = modded.mods.len(),
            n_adds = modded.additional.len(),
            "tracked modifications"
        );
        let filtered = NodeFilter::new(&modded, |i| matched[modded.to_index(i)].0 == usize::MAX);
        trace!(?matched, "matched nodes");
        let ext_start = found.len();
        found.extend(ConnectedGraphIter::new(&filtered).iter(mol).map(|bits| {
            let graph = GraphCompactor::<
                BitFiltered<
                    &NodeFilter<&ModdedGraph<G, MODDED_GRAPH_LEN>, _>,
                    usize,
                    ATOM_BIT_STORAGE,
                >,
            >::new(BitFiltered::new(&filtered, bits));
            let (ix, map) = self.insert_mol_atomic(&graph, true);
            let out = (0..map.len())
                .map(|i| modded.to_index(graph.node_map[map[i]]))
                .collect();
            (ix.index(), out)
        }));
        debug!(count = found.len() - ext_start, "unmatched sections split");
        for (i, ism) in &found[ext_start..] {
            let cmp = self.molecule(Ix::new(*i));
            for (cmp_i, &mol_i) in ism.iter().enumerate() {
                let graph_atom = cmp.get_atom(Ix::new(cmp_i)).unwrap();
                let mol_atom = mol.node_weight(mol.from_index(mol_i)).unwrap();
                if mol_atom.protons == 0 || graph_atom.protons != 0 {
                    matched[mol_i] = (*i, cmp_i);
                }
            }
        }
        found.sort();
        error!("BrokenMol linking not implemented");
        <Ix as IndexType>::max()
    }
}
