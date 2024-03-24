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
    index: Ix,
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
    pub fn get_atom_idx(&self, idx: NodeIndex<Ix>) -> Option<(Ix, Atom)> {
        let mut idx = idx.0.index();
        let arena = self.arena.get_arena();
        (idx < arena.parts.get(self.index.index())?.1.index()).then_some(())?;
        let mut oride = None;

        let mut ix = self.index;
        loop {
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => {
                    if oride.is_none() {
                        if let Some(a) = m.patch.get(&Ix::new(idx)) {
                            oride = Some(*a);
                        }
                    }
                    ix = m.base
                }
                MolRepr::Atomic(b) => {
                    let mut i = 0;
                    for w in b.as_slice() {
                        let o = w.count_ones() as usize;
                        if idx > o {
                            idx -= o;
                            i += usize::BITS as usize;
                        } else {
                            let ix = i + idx;
                            return Some((
                                Ix::new(ix),
                                oride.unwrap_or_else(|| {
                                    arena.graph()[petgraph::graph::NodeIndex::new(ix)]
                                }),
                            ));
                        }
                    }
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
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        for i in r {
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
                            break;
                        }
                    }
                }
            }
        }
    }

    /// Get an atom in this molecule.
    pub fn get_atom(&self, idx: NodeIndex<Ix>) -> Option<Atom> {
        let mut idx = idx.0.index();
        let arena = self.arena.get_arena();
        (idx < arena.parts.get(self.index.index())?.1.index()).then_some(())?;

        let mut ix = self.index;
        loop {
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => {
                    if let Some(a) = m.patch.get(&Ix::new(idx)) {
                        return Some(*a);
                    } else {
                        ix = m.base
                    }
                }
                MolRepr::Atomic(b) => {
                    let mut i = 0;
                    for w in b.as_slice() {
                        let o = w.count_ones() as usize;
                        if idx > o {
                            idx -= o;
                            i += usize::BITS as usize;
                        } else {
                            return Some(arena.graph()[petgraph::graph::NodeIndex::new(i + idx)]);
                        }
                    }
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
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        for i in r {
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
                            break;
                        }
                    }
                }
            }
        }
    }

    /// Get a bond between two atoms in this molecule.
    pub fn get_bond(&self, idx: EdgeIndex<Ix>) -> Option<Bond> {
        let mut idx0 = idx.source().index();
        let mut idx1 = idx.target().index();
        let arena = self.arena.get_arena();
        let s = arena.parts.get(self.index.index())?.1.index();
        (idx0 < s && idx1 < s).then_some(())?;

        let mut ix = self.index;
        loop {
            match &arena.parts[ix.index()].0 {
                MolRepr::Redirect(i) => ix = *i,
                MolRepr::Modify(m) => ix = m.base,
                MolRepr::Atomic(b) => {
                    let (mut i0, mut i1) = (0, 0);
                    let mut both = false;
                    for w in b.as_slice() {
                        let o = w.count_ones() as usize;
                        if both {
                            if idx0 > o {
                                idx0 -= o;
                                idx1 -= o;
                                i0 += usize::BITS as usize;
                                i1 += usize::BITS as usize;
                            } else if idx1 > o {
                                idx1 -= o;
                                i1 += usize::BITS as usize;
                                both = false;
                            } else {
                                break;
                            }
                        } else if idx1 > o {
                            idx1 -= o;
                            i1 += usize::BITS as usize;
                        } else {
                            break;
                        }
                    }
                    let g = arena.graph();
                    return Some(
                        g[g.find_edge(Ix::new(i0 + idx0).into(), Ix::new(i1 + idx1).into())?],
                    );
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
                    for p in &b.frags {
                        let mut s = arena.parts[p.index()].1.index();
                        let r = skips.get(p).unwrap_or(&empty);
                        for i in r {
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
                            first = Some(*p);
                        } else {
                            ix = *p;
                            break;
                        }
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

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> crate::molecule::ValueMolecule for Molecule<Ix, R> {
    fn get_atom(&self, idx: NodeIndex<Ix>) -> Option<Atom> {
        self.get_atom(idx)
    }
    fn get_bond(&self, idx: EdgeIndex<Ix>) -> Option<Bond> {
        self.get_bond(idx)
    }
}
