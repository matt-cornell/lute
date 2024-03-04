//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

use super::*;
use crate::graph::*;
use petgraph::data::DataMap;
use petgraph::graph::DefaultIx;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::collections::HashMap;
use std::hash::Hash;

const ATOM_BIT_STORAGE: usize = 2;

type Graph<Ix> = StableGraph<Atom, Bond, Undirected, Ix>;
type BSType = crate::utils::bitset::BitSet<usize, ATOM_BIT_STORAGE>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
struct InterFragBond<Ix> {
    an: Ix,
    ai: Ix,
    bn: Ix,
    bi: Ix,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BrokenMol<Ix> {
    frags: SmallVec<Ix, 8>,
    bonds: SmallVec<InterFragBond<Ix>, 8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ModdedMol<Ix> {
    base: Ix,
}

#[allow(clippy::large_enum_variant, dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
enum MolRepr<Ix> {
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
    parts: SmallVec<(MolRepr<Ix>, Ix), 16>,
}
impl<Ix: IndexType> Arena<Ix> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn graph(&self) -> &Graph<Ix> {
        &self.graph
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
    #[allow(clippy::needless_range_loop)]
    pub fn insert_mol<G>(&mut self, mol: G) -> Ix
    where
        G: Data<NodeWeight = Atom, EdgeWeight = Bond>
            + DataMap
            + GraphProp<EdgeType = Undirected>
            + GraphRef
            + GetAdjacencyMatrix
            + NodeCompactIndexable
            + IntoEdgesDirected
            + IntoNodeReferences,
        G::NodeId: Hash + Eq,
    {
        let compacted = (0..self.parts.len())
            .filter_map(|i| {
                if let (MolRepr::Atomic(a), _) = &self.parts[i] {
                    Some((
                        i,
                        GraphCompactor::<BitFiltered<&Graph<Ix>, usize, ATOM_BIT_STORAGE>>::new(
                            BitFiltered::<&Graph<Ix>, usize, ATOM_BIT_STORAGE>::new(
                                unsafe { &*std::ptr::addr_of!(self.graph) },
                                a.clone(),
                            ),
                        ),
                    ))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        // keep track of matched atoms so there's no overlap
        let mut matched: SmallVec<_, 16> = smallvec![(0, 0); mol.node_bound()];
        // keep track of found isomorphisms, don't try to handle them in the search
        let mut found = SmallVec::<_, 8>::new();
        for (k, (n, cmp)) in compacted.iter().enumerate() {
            for ism in
                subgraph_isomorphisms_iter(&cmp, &mol, &mut Atom::eq_or_r, &mut PartialEq::eq)
            {
                if ism.len() == mol.node_count() {
                    // simplest case, this molecule already exists
                    return Ix::new(*n);
                }
                if ism.iter().any(|&i| matched[i].0.index() != 0) {
                    continue;
                }
                ism.iter().enumerate().for_each(|(n, &i)| {
                    if self.graph[cmp.node_map[n]].protons != 0 {
                        matched[i] = (k + 1, n);
                    }
                });
                found.push((*n, ism));
            }
        }

        if found.is_empty() {
            // simple case: no subgraph isomorhpisms found
            let mut map = vec![NodeIndex::end(); mol.node_bound()];
            let mut bits = BSType::new();
            mol.node_references().for_each(|n| {
                let idx = self.graph.add_node(*n.weight());
                map[mol.to_index(n.id())] = idx;
                bits.set(idx.index(), true);
            });
            mol.edge_references().for_each(|e| {
                self.graph.add_edge(
                    map[mol.to_index(e.source())],
                    map[mol.to_index(e.target())],
                    *e.weight(),
                );
            });
            let out = self.parts.len();
            self.parts
                .push((MolRepr::Atomic(bits), Ix::new(mol.node_count())));
            Ix::new(out)
        } else {
            let mut nb = BSType::with_capacity(mol.node_bound());
            let mut needs_clear = false;
            let mut frags = SmallVec::with_capacity(found.len() + 1);
            frags.extend(found.iter().map(|(i, _)| Ix::new(*i)));
            let last = frags.len();
            frags.push(Ix::new(self.parts.len()));
            let mut bonds = SmallVec::new();

            for (n, (i, ism)) in found.iter().enumerate() {
                let cmp = &compacted[*i].1;
                for (j, &mi) in ism.iter().enumerate() {
                    let idx = cmp.node_map[j];
                    let atom = self.graph[idx];
                    let ma = mol.node_weight(mol.from_index(mi)).unwrap();
                    if atom.protons == 0 && ma.protons != 0 {
                        let bn = Ix::new(last);
                        let bi = Ix::new(mi);
                        let an = Ix::new(n);
                        let ai = Ix::new(idx.index());
                        bonds.push(InterFragBond { an, ai, bn, bi });
                    }
                    if atom.data.unknown() > 0 {
                        if needs_clear {
                            nb.clear();
                        }
                        needs_clear = true;
                        for c in self.graph.neighbors(idx) {
                            nb.set(c.index(), true);
                        }

                        // look for atoms bonded to this atom that don't have a corresponding node
                        // in the fragment, those *must* be an R-group
                        bonds.extend(mol.neighbors(mol.from_index(mi)).filter_map(|mi| {
                            let mix = mol.to_index(mi);
                            if nb.get(mix) {
                                None // already accounted for
                            } else {
                                Some(InterFragBond {
                                    bn: Ix::new(last),
                                    bi: Ix::new(mix),
                                    an: Ix::new(n),
                                    ai: Ix::new(idx.index()),
                                })
                            }
                        }));
                    }
                }
            }

            if matched.iter().all(|m| m.0 != 0) {
                // all atoms make up other molecules
                for bond in &mut bonds {
                    let (bn, bi) = matched[bond.bi.index()];
                    // keep an < bn
                    if bn >= bond.an.index() {
                        bond.bn = Ix::new(bn);
                        bond.bi = Ix::new(bi);
                    } else {
                        bond.bn = std::mem::replace(&mut bond.an, Ix::new(bn));
                        bond.bi = std::mem::replace(&mut bond.ai, Ix::new(bi));
                    }
                }
            } else {
                let mut map = vec![NodeIndex::end(); mol.node_bound()];
                let mut bits = BSType::new();
                let mut count = 0;
                mol.node_references().for_each(|n| {
                    let ni = mol.to_index(n.id());
                    if matched[ni].0.index() == 0 {
                        count += 1;
                        let idx = self.graph.add_node(*n.weight());
                        map[ni] = idx;
                        bits.set(idx.index(), true);
                    }
                });
                mol.edge_references().for_each(|e| {
                    let si = mol.to_index(e.source());
                    let ti = mol.to_index(e.target());
                    if matched[ti].0.index() != 0 || matched[si].0.index() != 0 {
                        return;
                    }
                    self.graph.add_edge(
                        map[mol.to_index(e.source())],
                        map[mol.to_index(e.target())],
                        *e.weight(),
                    );
                });

                for bond in &mut bonds {
                    let ix = map[bond.bi.index()].index();
                    if ix == <Ix as IndexType>::max().index() {
                        // handle the case of a node that's bonded directly to another atom
                        let (bn, bi) = matched[bond.bi.index()];
                        if bn >= bond.an.index() {
                            bond.bn = Ix::new(bn);
                            bond.bi = Ix::new(bi);
                        } else {
                            bond.bn = std::mem::replace(&mut bond.an, Ix::new(bn));
                            bond.bi = std::mem::replace(&mut bond.ai, Ix::new(bi));
                        }
                    } else {
                        bond.bi = Ix::new(ix);
                    }
                }

                let filtered = NodeFilter::new(mol, |n| matched[mol.to_index(n)].0 == 0);

                let mut iter = DisjointGraphIter::new(&filtered).iter(&filtered).peekable();
                let g0 = iter.by_ref().next().unwrap();
                if iter.peek().is_some() {
                    // bits in a `usize`
                    let usbits = std::mem::size_of::<usize>() * 8;
                    // mapping from bit index to the index in `out`
                    let mut map = HashMap::with_capacity(g0.count_ones());
                    let mut insert_bits = |bits: BSType, idx| {
                        let mut count = 0;
                        for (byte, &w) in bits.as_slice().iter().enumerate() {
                            for bit in 0..usbits {
                                if w & (1 << bit) != 0 {
                                    map.insert(byte * usbits + bit, idx);
                                    count += 1;
                                }
                            }
                        }
                        // graph, node count, initial index, inclusion marker
                        (
                            GraphCompactor::<
                                BitFiltered<&NodeFilter<G, _>, usize, ATOM_BIT_STORAGE>,
                            >::new(BitFiltered::new(&filtered, bits)),
                            count,
                            idx,
                            None,
                        )
                    };
                    let mut chunks: SmallVec<_, 4> = smallvec![insert_bits(g0, 0)];
                    chunks.extend(
                        iter.by_ref()
                            .enumerate()
                            .map(|(n, bits)| insert_bits(bits, n + 1)),
                    );
                    chunks.sort_by_key(|p| p.1);
                    // mapping from the old fragment index to the one after the sort. fragments
                    // isomorphic to a previous fragment will not be referred to
                    let mut fmap: SmallVec<_, 4> =
                        smallvec![(usize::MAX, usize::MAX); chunks.len()];
                    fmap[chunks[0].2] = (0, 0);
                    // removal list, we can't do it inline because `self.graph` is borrowed
                    let mut rm = vec![];
                    // last index
                    let mut li = 0;
                    for i in 1..chunks.len() {
                        if let Some(ism) = find_isomorphism_matching(
                            &&chunks[last].0,
                            &&chunks[i].0,
                            &mut Atom::eq_match_r,
                            &mut PartialEq::eq,
                        ) {
                            chunks[i].3 = Some(ism);
                            fmap[chunks[i].2] = (li, i);
                            for (byte, &w) in chunks[i].0.graph.filter.as_slice().iter().enumerate()
                            {
                                for bit in 0..usbits {
                                    if w & (1 << bit) != 0 {
                                        rm.push(Ix::new(byte * usbits + bit));
                                    }
                                }
                            }
                        } else {
                            li = i;
                            fmap[chunks[i].2] = (i, i);
                        }
                    }

                    for bond in &mut bonds {
                        if bond.bn.index() == last {
                            let old = map[&bond.bi.index()];
                            let new = fmap[old];
                            bond.bn = Ix::new(last + new.1);
                            bond.bi = if let Some(ism) = chunks[old].3.as_deref() {
                                // chunks[old].0.inv_map[...] -- get compacted index from old chunk
                                // ism[x] -- run through isomorphism
                                // chunks[new].0.node_map[x] -- get new index
                                // Ix::new(x.index()) -- convert the NodeIndex into an Ix
                                Ix::new(mol.to_index(
                                    chunks[new.0].0.node_map[ism
                                        [chunks[old].0.inv_map[&mol.from_index(bond.bi.index())]]],
                                ))
                            } else {
                                bond.bi
                            };
                        }
                    }

                    let mut frag_count = 0;
                    self.parts
                        .extend(chunks.into_iter().filter_map(|(g, c, _, i)| {
                            i.is_none().then(|| {
                                frag_count += 1;
                                (MolRepr::Atomic(g.graph.filter), Ix::new(c))
                            })
                        }));
                    frags.extend((last..(last + frag_count)).map(Ix::new));

                    for n in rm {
                        self.graph.remove_node(n.into());
                    }
                } else {
                    self.parts.push((MolRepr::Atomic(g0), Ix::new(count)));
                }
            };

            let out = Ix::new(self.parts.len());
            self.parts.push((
                MolRepr::Broken(BrokenMol { frags, bonds }),
                Ix::new(mol.node_count()),
            ));

            out
        }
    }
}
