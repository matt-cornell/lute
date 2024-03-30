use std::collections::BTreeSet;

use super::arena::*;
use super::*;
pub mod graph_traits;
mod node_impls;
use hybridmap::HybridMap;
pub use node_impls::*;

/// A `Molecule` acts like a graph, and can have graph algorithms used on it. It's immutable, with all
/// mutations making (efficient) copies.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct Molecule<Ix, R> {
    arena: R,
    pub index: Ix,
}
impl<Ix, R> Molecule<Ix, R> {
    pub fn from_arena<'a, 'b: 'a, A: ArenaAccessible<Ix = Ix, Access<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor(),
            index,
        }
    }

    pub fn from_mut_arena<'a, 'b: 'a, A: ArenaAccessibleMut<Ix = Ix, AccessMut<'a> = R> + 'a>(
        arena: &'b A,
        index: Ix,
    ) -> Self
    where
        Self: 'a,
    {
        Self {
            arena: arena.get_accessor_mut(),
            index,
        }
    }

    pub fn arena(&self) -> R::Ref<'_>
    where
        R: ArenaAccessor<Ix = Ix>,
    {
        self.arena.get_arena()
    }

    pub fn arena_mut(&self) -> R::RefMut<'_>
    where
        R: ArenaAccessorMut<Ix = Ix>,
    {
        self.arena.get_arena_mut()
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Molecule<Ix, R> {
    /// Check if this molecule contains the underlying group.
    pub fn contains(&self, group: Ix) -> bool {
        self.arena.get_arena().contains_group(self.index, group)
    }

    /// Get an atom, along with the corresponding node index in the underlying graph.
    /// Slower than `get_atom` because it can't short-circuit at `MolRepr::Modify`s.
    #[inline]
    #[instrument(level = "debug", skip_all, fields(mol_idx = self.index.index(), idx))]
    pub fn get_atom_idx(&self, idx: impl Into<NodeIndex<Ix>>) -> Option<(Ix, Atom)> {
        let idx = idx.into();
        Span::current().record("idx", &idx.0.index());
        self.get_atom_idx_impl(idx)
    }
    fn get_atom_idx_impl(&self, idx: NodeIndex<Ix>) -> Option<(Ix, Atom)> {
        let mut idx = idx.0.index();
        let arena = self.arena.get_arena();
        (idx < arena.parts.get(self.index.index())?.1.index()).then_some(())?;
        let mut oride = None;

        let mut ix = self.index;
        loop {
            trace!(mol_idx = ix.index(), atom_idx = idx, "searching for atom");
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => {
                    if oride.is_none() {
                        if let Ok(a) = m.patch.binary_search_by_key(&Ix::new(idx), |m| m.0) {
                            oride = Some(m.patch[a].1);
                        }
                    }
                    ix = m.base
                }
                MolRepr::Atomic(b) => {
                    let i = b.nth(idx)?;
                    return Some((
                        Ix::new(i),
                        oride.unwrap_or_else(|| arena.graph()[petgraph::graph::NodeIndex::new(i)]),
                    ));
                }
                MolRepr::Broken(b) => {
                    let empty = BTreeSet::new();
                    let mut skips = HybridMap::<Ix, BTreeSet<Ix>, 4>::new();
                    for i in &b.bonds {
                        if arena
                            .molecule(i.an)
                            .get_atom(NodeIndex(i.ai))
                            .unwrap()
                            .protons
                            == 0
                        {
                            if let Some(s) = skips.get_mut(&i.an) {
                                s.insert(i.ai);
                            } else {
                                skips.insert(i.an, [i.ai].into());
                            }
                        }
                        if arena
                            .molecule(i.bn)
                            .get_atom(NodeIndex(i.bi))
                            .unwrap()
                            .protons
                            == 0
                        {
                            if let Some(s) = skips.get_mut(&i.bn) {
                                s.insert(i.bi);
                            } else {
                                skips.insert(i.bn, [i.bi].into());
                            }
                        }
                    }
                    let mut found = false;
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        let mut ic = 0;
                        for i in r {
                            if i.index() < ic {
                                ic += 1;
                            }
                            if i.index() < s {
                                s -= 1;
                            } else {
                                break;
                            }
                        }
                        if idx > s {
                            idx -= s;
                        } else {
                            ix = *p;
                            idx += ic;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return None;
                    }
                }
            }
        }
    }

    /// Get an atom in this molecule.
    #[inline]
    #[instrument(level = "debug", skip_all, fields(mol_idx = self.index.index(), idx))]
    pub fn get_atom(&self, idx: impl Into<NodeIndex<Ix>>) -> Option<Atom> {
        let idx = idx.into();
        Span::current().record("idx", &idx.0.index());
        self.get_atom_impl(idx)
    }
    fn get_atom_impl(&self, idx: NodeIndex<Ix>) -> Option<Atom> {
        let mut idx = idx.0.index();
        let arena = self.arena.get_arena();
        (idx < arena.parts.get(self.index.index())?.1.index()).then_some(())?;

        let mut ix = self.index;
        loop {
            trace!(mol_idx = ix.index(), atom_idx = idx, "searching for atom");
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => {
                    if let Ok(a) = m.patch.binary_search_by_key(&Ix::new(idx), |m| m.0) {
                        return Some(m.patch[a].1);
                    } else {
                        ix = m.base
                    }
                }
                MolRepr::Atomic(b) => {
                    let i = b.nth(idx)?;
                    return Some(arena.graph()[petgraph::graph::NodeIndex::new(i)]);
                }
                MolRepr::Broken(b) => {
                    let empty = BTreeSet::new();
                    let mut skips = HybridMap::<Ix, BTreeSet<Ix>, 4>::new();
                    for i in &b.bonds {
                        if arena
                            .molecule(i.an)
                            .get_atom(NodeIndex(i.ai))
                            .unwrap()
                            .protons
                            == 0
                        {
                            if let Some(s) = skips.get_mut(&i.an) {
                                s.insert(i.ai);
                            } else {
                                skips.insert(i.an, [i.ai].into());
                            }
                        }
                        if arena
                            .molecule(i.bn)
                            .get_atom(NodeIndex(i.bi))
                            .unwrap()
                            .protons
                            == 0
                        {
                            if let Some(s) = skips.get_mut(&i.bn) {
                                s.insert(i.bi);
                            } else {
                                skips.insert(i.bn, [i.bi].into());
                            }
                        }
                    }
                    let mut found = false;
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        let mut ic = 0;
                        for i in r {
                            if i.index() < ic {
                                ic += 1;
                            }
                            if i.index() < s {
                                s -= 1;
                            } else {
                                break;
                            }
                        }
                        if idx > s {
                            idx -= s;
                        } else {
                            ix = *p;
                            idx += ic;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return None;
                    }
                }
            }
        }
    }

    /// Get a bond between two atoms in this molecule.
    #[inline]
    #[instrument(level = "debug", skip_all, fields(mol_idx = self.index.index(), idx))]
    pub fn get_bond(&self, idx: impl Into<EdgeIndex<Ix>>) -> Option<Bond> {
        let idx = idx.into();
        let span = Span::current();
        span.record("idx.source", &idx.source().index());
        span.record("idx.target", &idx.target().index());
        self.get_bond_impl(idx.into())
    }
    fn get_bond_impl(&self, idx: EdgeIndex<Ix>) -> Option<Bond> {
        let mut idx0 = idx.source().index();
        let mut idx1 = idx.target().index();
        let arena = self.arena.get_arena();
        let s = arena.parts.get(self.index.index())?.1.index();
        (idx0 < s && idx1 < s).then_some(())?;

        let mut ix = self.index;
        loop {
            trace!(
                mol_idx = ix.index(),
                atom_idx0 = idx0,
                atom_idx1 = idx1,
                "searching for bond"
            );
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => ix = m.base,
                MolRepr::Atomic(b) => {
                    // inefficient, but works
                    let i0 = b.nth(idx0)?;
                    let i1 = b.nth(idx1)?;
                    // let [i0, i1] = b.nth_many_short([idx0, idx1])?;
                    let g = arena.graph();
                    return Some(g[g.find_edge(Ix::new(i0).into(), Ix::new(i1).into())?]);
                }
                MolRepr::Broken(b) => {
                    let empty = BTreeSet::new();
                    let mut skips = HybridMap::<Ix, BTreeSet<Ix>, 4>::new();
                    let mut ibs = HybridMap::<InterFragBond<Ix>, Bond, 4>::new();
                    for i in &b.bonds {
                        let (ix, a) = arena.molecule(i.an).get_atom_idx(NodeIndex(i.ai)).unwrap();
                        let mut weight = None;
                        if a.protons == 0 {
                            weight = weight.or_else(|| {
                                arena.graph().edges(ix.into()).next().map(|e| *e.weight())
                            });
                            if let Some(s) = skips.get_mut(&i.an) {
                                s.insert(i.ai);
                            } else {
                                skips.insert(i.an, [i.ai].into());
                            }
                        }
                        let (ix, a) = arena.molecule(i.bn).get_atom_idx(NodeIndex(i.bi)).unwrap();
                        if a.protons == 0 {
                            weight = weight.or_else(|| {
                                arena.graph().edges(ix.into()).next().map(|e| *e.weight())
                            });
                            if let Some(s) = skips.get_mut(&i.bn) {
                                s.insert(i.bi);
                            } else {
                                skips.insert(i.bn, [i.bi].into());
                            }
                        }
                        ibs.insert(*i, weight.unwrap_or(Bond::Single));
                    }
                    let mut first = None;
                    let mut found = false;
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        let (mut i0c, mut i1c) = (0, 0);
                        for i in r {
                            if i.index() < i0c {
                                i0c += 1;
                            }
                            if i.index() < i1c {
                                i1c += 1;
                            }
                            if i.index() < s {
                                s -= 1;
                            } else {
                                break;
                            }
                        }
                        if let Some(first) = first {
                            if idx1 > s {
                                idx1 -= s;
                            } else {
                                idx1 += i1c;
                                return ibs
                                    .get(&InterFragBond {
                                        an: first,
                                        ai: Ix::new(idx0),
                                        bn: *p,
                                        bi: Ix::new(idx1),
                                    })
                                    .copied();
                            }
                        } else if idx0 > s {
                            idx0 -= s;
                            idx1 -= s;
                        } else if idx1 > s {
                            idx1 -= s;
                            idx0 += i0c;
                            first = Some(*p);
                        } else {
                            ix = *p;
                            idx0 += i0c;
                            idx1 += i1c;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return None;
                    }
                }
            }
        }
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>, S: ArenaAccessor<Ix = Ix>> PartialEq<Molecule<Ix, S>>
    for Molecule<Ix, R>
{
    fn eq(&self, other: &Molecule<Ix, S>) -> bool {
        self.index == other.index && std::ptr::eq(&*self.arena(), &*other.arena())
    }
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> crate::graph::misc::DataValueMap
    for Molecule<Ix, R>
{
    fn node_weight(&self, idx: NodeIndex<Ix>) -> Option<Atom> {
        self.get_atom(idx)
    }
    fn edge_weight(&self, idx: EdgeIndex<Ix>) -> Option<Bond> {
        self.get_bond(idx)
    }
}
