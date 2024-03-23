use crate::core::*;
use crate::empirical::*;
use modular_bitfield::prelude::*;
use petgraph::data::*;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::cell::UnsafeCell;
use std::cmp::Ordering;
use std::fmt::{self, Debug, Formatter};

const ELECTRON_COUNTS: &[u8] = &[0, 2, 10, 18, 30, 36, 54, 86, 118];

fn unshared_electrons(protons: u8, charge: i8, bonded: u8) -> u8 {
    let elecs = (protons as i8 - charge) as u8;
    let shell = ELECTRON_COUNTS.binary_search(&elecs).unwrap_or_else(|e| e);
    if shell == 0 {
        return 0;
    }
    let lower = ELECTRON_COUNTS[shell];
    let upper = ELECTRON_COUNTS[shell + 1];
    // noble gas
    if elecs == lower {
        return 8 - bonded;
    }
    // p-block
    if let Some(res) = elecs.checked_sub(upper - 6) {
        return res + 2 - bonded;
    }
    // s-block
    if elecs <= lower + 2 {
        return elecs - lower - bonded;
    }
    // some metal
    2_u8.saturating_sub(bonded)
}

#[bitfield]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomInfoFields {
    pub bonded: B7,
    pub neighbors_pi: bool,
}

/// `AtomInfo` holds some more data about the atoms this atom is bonded to. It can calculate things
/// like unshared electrons and hybridization because of this.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomInfo {
    pub atom: Atom,
    pub fields: AtomInfoFields,
}

impl AtomInfo {
    /// The number of unshared (non-bonding) valence electrons.
    pub fn unshared(&self) -> u8 {
        unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded())
    }
    /// The hybridization of this atom. 1 is s, 2 is sp, 3 is sp2, etc.
    pub fn hybridization(&self) -> u8 {
        let unshared =
            unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded());
        let groups = self.atom.data.total_bonds() + (unshared + 1) / 2;
        match groups {
            4 => {
                if self.fields.neighbors_pi() {
                    3
                } else {
                    4
                }
            }
            0..=6 => groups,
            g => panic!("atom has {g} bonding groups, which shouldn't be possible"),
        }
    }
    /// The number of delocalized pi electrons available for an aromatic system.
    pub fn pi_donor(&self) -> u8 {
        let unshared =
            unshared_electrons(self.atom.protons, self.atom.charge, self.fields.bonded());
        let mut groups = self.atom.data.total_bonds() + unshared.div_ceil(2);
        if self.fields.neighbors_pi() && groups == 4 {
            groups -= 1;
        }
        unshared.saturating_sub(groups * 2)
    }
}

/// It's a lot easier to return an atom by value, since we have to handle borrows and locking, etc.
/// We rely on this instead of `DataMap` so that `Molecule` can implement the full set of functions.
pub trait ValueMolecule: Data<NodeWeight = Atom, EdgeWeight = Bond> {
    fn get_atom(&self, id: Self::NodeId) -> Option<Atom>;
    fn get_bond(&self, id: Self::EdgeId) -> Option<Bond>;
}
impl<T: DataMap<NodeWeight = Atom, EdgeWeight = Bond>> ValueMolecule for T {
    fn get_atom(&self, id: Self::NodeId) -> Option<Atom> {
        self.node_weight(id).copied()
    }
    fn get_bond(&self, id: Self::EdgeId) -> Option<Bond> {
        self.edge_weight(id).copied()
    }
}

/// Because we may need to calculate more, the CIP structure needs interior mutability
#[derive(Debug, Clone)]
struct CipPriorityInner<G: Visitable> {
    pub weights: SmallVec<(u16, u8), 24>,
    pub indices: SmallVec<u16, 4>,
    graph: G,
    edge: SmallVec<(G::NodeId, u8), 4>,
    seen: G::Map,
}
impl<G: Visitable> CipPriorityInner<G> {
    pub fn new(graph: G, node: G::NodeId) -> Self {
        let map = graph.visit_map();
        Self::with_seen(graph, node, map)
    }
    pub fn with_seen(graph: G, node: G::NodeId, seen: G::Map) -> Self {
        Self {
            weights: smallvec![],
            indices: smallvec![],
            edge: smallvec![(node, 1)],
            graph,
            seen,
        }
    }
    pub fn on_branch(graph: G, node: G::NodeId, center: G::NodeId) -> Self {
        let mut map = graph.visit_map();
        map.visit(center);
        Self::with_seen(graph, node, map)
    }
}
impl<G: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>> CipPriorityInner<G> {
    fn compute_step(&mut self) -> bool {
        let edge = std::mem::take(&mut self.edge);
        if edge.is_empty() {
            return false;
        }
        for (n, w) in edge {
            self.seen.visit(n);
            self.edge.extend(self.graph.edges(n).filter_map(|e| {
                let n2 = if e.source() == n {
                    e.target()
                } else {
                    e.source()
                };
                (!self.seen.is_visited(&n2)).then_some((n2, e.weight().bond_count().floor() as u8))
            }));
            let a = self.graph.node_weight(n).unwrap();
            self.weights
                .push((((a.protons as u16) << 9) + a.isotope, w));
        }
        let mut idx = self.indices.last().map_or(0, |x| *x as usize);
        self.weights[idx..].sort();
        // yeesh... worst-case quadratic time here. There's gotta be a better way, but this will probably only be working on like 20 atoms *at most*
        while idx + 1 < self.weights.len() {
            let (w2, n2) = self.weights[idx + 1];
            let (w1, n1) = &mut self.weights[idx];
            if *w1 == w2 {
                *n1 += n2;
                self.weights.remove(idx + 1);
            } else {
                idx += 1;
            }
        }
        self.indices.push(self.weights.len() as _);
        !self.edge.is_empty()
    }
    fn cmp_impl<H: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>>(
        &mut self,
        other: &mut CipPriorityInner<H>,
        at: usize,
        out: &mut Ordering,
    ) -> bool {
        let ret = (self.indices.len() > at || (self.indices.len() == at && self.compute_step()))
            | (other.indices.len() > at || (other.indices.len() == at && other.compute_step()));
        let lhs: &[_] = self.indices.get(at).map_or(&[], |i| {
            let i = *i as usize;
            if at == 0 {
                &self.weights[..i]
            } else {
                &self.weights[(self.indices[at - 1] as _)..i]
            }
        });
        let rhs: &[_] = other.indices.get(at).map_or(&[], |i| {
            let i = *i as usize;
            if at == 0 {
                &other.weights[..i]
            } else {
                &other.weights[(other.indices[at - 1] as _)..i]
            }
        });
        *out = lhs.cmp(&rhs);
        ret && *out == Ordering::Equal
    }
    pub fn cmp<H: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>>(
        &mut self,
        other: &mut CipPriorityInner<H>,
    ) -> Ordering {
        let mut out = Ordering::Equal;
        let mut i = 0;
        while self.cmp_impl(other, i, &mut out) {
            i += 1;
        }
        out
    }
}

/// CIP priority can be partially calculated, this allows only necessary calculations to be done
pub struct CipPriority<G: Visitable>(UnsafeCell<CipPriorityInner<G>>);
impl<G: Visitable> CipPriority<G> {
    pub fn new(graph: G, node: G::NodeId) -> Self {
        Self(UnsafeCell::new(CipPriorityInner::new(graph, node)))
    }
    pub fn with_seen(graph: G, node: G::NodeId, seen: G::Map) -> Self {
        Self(UnsafeCell::new(CipPriorityInner::with_seen(
            graph, node, seen,
        )))
    }
    pub fn on_branch(graph: G, node: G::NodeId, center: G::NodeId) -> Self {
        Self(UnsafeCell::new(CipPriorityInner::on_branch(
            graph, node, center,
        )))
    }
}
impl<G: Visitable + Debug> Debug for CipPriority<G>
where
    G::NodeId: Debug,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        unsafe {
            let inner = &*self.0.get();
            f.debug_struct("CipPriority")
                .field("weights", &inner.weights)
                .field("indices", &inner.indices)
                .field("edge", &inner.edge)
                .finish_non_exhaustive()
        }
    }
}
impl<
        G1: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>,
        G2: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>,
    > PartialEq<CipPriority<G2>> for CipPriority<G1>
{
    fn eq(&self, other: &CipPriority<G2>) -> bool {
        <Self as PartialOrd<_>>::partial_cmp(&self, other) == Some(Ordering::Equal)
    }
}
impl<G: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>> Eq
    for CipPriority<G>
{
}
impl<
        G1: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>,
        G2: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>,
    > PartialOrd<CipPriority<G2>> for CipPriority<G1>
{
    fn partial_cmp(&self, other: &CipPriority<G2>) -> Option<Ordering> {
        unsafe {
            let lhs = &mut *self.0.get();
            let rhs = &mut *other.0.get();
            Some(lhs.cmp(rhs))
        }
    }
}
impl<G: Visitable + IntoEdges + DataMap<NodeWeight = Atom, EdgeWeight = Bond>> Ord
    for CipPriority<G>
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

/// Most of the molecule operations. This should be implemented for whatever is necessary, and
/// impossible to implement for additional types.
pub trait Molecule<'a>
where
    &'a Self: GraphProp<EdgeType = petgraph::Undirected>
        + Data<NodeWeight = Atom, EdgeWeight = Bond>
        + 'a,
{
    /// Find the mass of a molecule.
    fn mass(&'a self) -> f32
    where
        &'a Self: IntoNodeReferences,
    {
        self.node_references()
            .map(|a| {
                a.weight().mass()
                    + a.weight().data.hydrogen() as f32 * crate::atom_info::ATOM_DATA[1].mass
            })
            .sum()
    }

    fn empirical(&'a self) -> EmpiricalFormula
    where
        &'a Self: IntoNodeReferences,
    {
        self.node_references().map(|a| *a.weight()).collect()
    }

    /// Get the `AtomData` for a molecule in the graph.
    fn atom_data(&'a self, id: AtomId<'a, Self>) -> AtomInfo
    where
        &'a Self: ValueMolecule + IntoEdges,
    {
        let atom = self.get_atom(id).unwrap();
        let mut neighbors_pi = false;
        let mut bond_electrons = 0.0;
        for edge in self.edges(id) {
            bond_electrons += edge.weight().bond_count();
            let neighbor = if edge.source() == id {
                if edge.target() == id {
                    panic!("Molecule should not contain self-loops!")
                } else {
                    edge.target()
                }
            } else {
                edge.source()
            };
            if !neighbors_pi {
                neighbors_pi = self.edges(neighbor).any(|e| e.weight().bond_count() > 1.0);
            }
        }
        AtomInfo {
            atom,
            fields: AtomInfoFields::new()
                .with_bonded(bond_electrons.floor() as _)
                .with_neighbors_pi(neighbors_pi),
        }
    }

    /// Get CIP priority for an atom.
    fn cip_priority(&'a self, atom: AtomId<'a, Self>, center: AtomId<'a, Self>) -> CipPriority<&Self>
    where
        Self: Visitable + Sized,
    {
        CipPriority::on_branch(self, atom, center)
    }

    /// Resolve tautomers. This also ensures that aromaticicty is tracked, and rings have the
    /// "correct" conjugation.
    fn resolve_tautomer(&'a mut self)
    where
        Self: IntoEdges + DataMapMut,
    {
        todo!()
    }
}

type AtomId<'a, T> = <&'a T as GraphBase>::NodeId;

impl<'a, T: 'a> Molecule<'a> for T
where
    &'a T: GraphProp<EdgeType = petgraph::Undirected> + Data<NodeWeight = Atom, EdgeWeight = Bond>,
{
    
}
