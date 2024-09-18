use crate::utils::bitset::BitSet;

use super::*;
use petgraph::{visit::*, Direction};
use std::iter::Map;
use std::ops::Range;

impl<Ix: Copy + Ord, R> GraphBase for Molecule<Ix, R> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
}

impl<Ix: Copy + Ord, R: Copy> GraphRef for Molecule<Ix, R> {}

impl<Ix: Copy + Ord, R> GraphProp for Molecule<Ix, R> {
    type EdgeType = petgraph::Undirected;
}

impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> NodeCount for Molecule<Ix, R> {
    fn node_count(&self) -> usize {
        self.arena.get_arena().frags[self.index.0.index()].size()
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> NodeIndexable for Molecule<Ix, R> {
    fn node_bound(&self) -> usize {
        self.arena.get_arena().frags[self.index.0.index()].size()
    }
    fn to_index(&self, a: NodeIndex<Ix>) -> usize {
        a.0.index()
    }
    fn from_index(&self, i: usize) -> NodeIndex<Ix> {
        NodeIndex(Ix::new(i))
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> NodeCompactIndexable for Molecule<Ix, R> {}

impl<Ix: IndexType + Ord, R> Data for Molecule<Ix, R> {
    type NodeWeight = Atom;
    type EdgeWeight = Bond;
}

impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeIdentifiers
    for Molecule<Ix, R>
{
    type NodeIdentifiers = Map<Range<usize>, fn(usize) -> NodeIndex<Ix>>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        (0..self.node_count()).map(|ix| NodeIndex(Ix::new(ix)))
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeReferences for Molecule<Ix, R> {
    type NodeRef = NodeReference<Ix>;
    type NodeReferences = iter::NodeReferences<Ix, R>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeReferences {
            range: 0..self.node_bound(),
            mol_idx: self.index.0,
            arena: self.arena,
        }
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdgeReferences for Molecule<Ix, R> {
    type EdgeRef = EdgeReference<Ix>;
    type EdgeReferences = iter::EdgeReferences<Ix, R>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeReferences::new(self.index.0, self.arena)
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdges for Molecule<Ix, R> {
    type Edges = iter::EdgesDirected<Ix, R>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::EdgesDirected(
            iter::Edges::new(self.index.0, a.0, self.arena),
            Direction::Outgoing,
        )
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdgesDirected for Molecule<Ix, R> {
    type EdgesDirected = iter::EdgesDirected<Ix, R>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgesDirected(iter::Edges::new(self.index.0, a.0, self.arena), dir)
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoNeighbors for Molecule<Ix, R> {
    type Neighbors = iter::Neighbors<Ix, R>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::Neighbors(iter::Edges::new(self.index.0, a.0, self.arena))
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> IntoNeighborsDirected
    for Molecule<Ix, R>
{
    type NeighborsDirected = iter::Neighbors<Ix, R>;

    fn neighbors_directed(self, n: Self::NodeId, _dir: Direction) -> Self::NeighborsDirected {
        iter::Neighbors(iter::Edges::new(self.index.0, n.0, self.arena))
    }
}
impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix> + Copy> GetAdjacencyMatrix for Molecule<Ix, R> {
    type AdjMatrix = crate::utils::bitset::BitSet<usize, 1>;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        let mut out = Self::AdjMatrix::new();
        for e in self.edge_references() {
            let a = e.source().0.index();
            let b = e.target().0.index();
            if a == b {
                continue;
            }
            let (s, t) = if a < b { (a, b) } else { (b, a) };
            let i = t * (t - 1) / 2 + s;
            out.set(i, true);
        }
        out
    }
    fn is_adjacent(&self, matrix: &Self::AdjMatrix, a: Self::NodeId, b: Self::NodeId) -> bool {
        if a == b {
            return false;
        }
        let a = a.0.index();
        let b = b.0.index();
        let (s, t) = if a < b { (a, b) } else { (b, a) };
        let i = t * (t - 1) / 2 + s;
        matrix.get(i)
    }
}
impl<Ix: IndexType + Ord, R> Visitable for Molecule<Ix, R> {
    type Map = BitSet<usize, 1>;

    fn reset_map(&self, map: &mut Self::Map) {
        map.clear();
    }
    fn visit_map(&self) -> Self::Map {
        BitSet::new()
    }
}

pub mod iter {
    use super::*;
    use petgraph::stable_graph::WalkNeighbors;
    use smallvec::{smallvec, SmallVec};
    use std::fmt::{self, Debug, Formatter};

    // TODO: make this more efficient-- lookup is cheap during traversal
    #[derive(Debug, Clone)]
    pub struct NodeReferences<Ix, R> {
        pub(crate) range: Range<usize>,
        pub(crate) mol_idx: Ix,
        pub(crate) arena: R,
    }
    impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> Iterator for NodeReferences<Ix, R> {
        type Item = NodeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            self.range
                .next()
                .map(|i| NodeReference::new(self.mol_idx.into(), NodeIndex(Ix::new(i)), self.arena))
        }
    }

    #[derive(Clone)]
    pub struct EdgeReferences<Ix, R> {
        stack: SmallVec<(Ix, Ix), 3>,
        // `StableGraph` doesn't have walkers, so we need to buffer our references
        buffer: SmallVec<EdgeReference<Ix>, 3>,
        arena: R,
    }
    impl<Ix: IndexType, R> EdgeReferences<Ix, R> {
        pub fn new(mol_idx: Ix, arena: R) -> Self {
            Self {
                stack: smallvec![(mol_idx, Ix::new(0))],
                buffer: smallvec![],
                arena,
            }
        }
    }
    impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> Iterator for EdgeReferences<Ix, R> {
        type Item = EdgeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(eref) = self.buffer.pop() {
                return Some(eref);
            }
            let arena = self.arena.get_arena();
            while let Some((idx, off)) = self.stack.pop() {
                match arena.frags[idx.index()].repr {
                    MolRepr::Modify(ModdedMol { base, .. }) => self.stack.push((base, off)),
                    MolRepr::Broken(BrokenMol {
                        ref frags,
                        ref bonds,
                    }) => {
                        self.stack.reserve(frags.len());
                        let offsets = frags
                            .iter()
                            .scan(off.index(), |count, idx| {
                                let old = *count;
                                *count += arena.frags[idx.index()].size();
                                self.stack.push((*idx, Ix::new(old)));
                                Some(old)
                            })
                            .collect::<SmallVec<_, 3>>();
                        let mut it = bonds.iter().map(|b| {
                            EdgeReference::with_weight(
                                Ix::new(offsets[b.an.index()] + b.ai.index()),
                                Ix::new(offsets[b.bn.index()] + b.bi.index()),
                                Bond::Single,
                            )
                        });
                        if let Some(ret) = it.next() {
                            self.buffer.extend(it);
                            return Some(ret);
                        }
                    }
                    MolRepr::Atomic(ref b) => {
                        let offset = off.index();
                        let mut it = arena.graph().edge_references().filter_map(|e| {
                            let s = b.index(e.source().index())?;
                            let t = b.index(e.target().index())?;
                            Some(EdgeReference::with_weight(
                                Ix::new(s + offset),
                                Ix::new(t + offset),
                                *e.weight(),
                            ))
                        });
                        if let Some(ret) = it.next() {
                            self.buffer.extend(it);
                            return Some(ret);
                        }
                    }
                }
            }
            None
        }
    }

    #[derive(Clone)]
    enum State<Ix: IndexType> {
        Uninit,
        Broken(SmallVec<Ix, 4>),
        Atomic(WalkNeighbors<Ix>),
    }
    impl<Ix: IndexType> Debug for State<Ix> {
        fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
            match self {
                Self::Uninit => f.write_str("Uninit"),
                Self::Broken(_) => f.write_str("Broken"),
                Self::Atomic(_) => f.write_str("Atomic"),
            }
        }
    }

    /// This is probably *not* what you want to use! It returns all neighboring edges, and is more
    /// for implementation than direct use!
    #[derive(Debug, Clone)]
    pub struct Edges<Ix: IndexType + Ord, R> {
        orig_idx: Ix,
        mol_idx: Ix,
        atom_idx: Ix,
        offset: Ix,
        arena: R,
        state: State<Ix>,
    }
    impl<Ix: IndexType + Ord, R> Edges<Ix, R> {
        pub fn new(mol_idx: Ix, atom_idx: Ix, arena: R) -> Self {
            Self {
                mol_idx,
                atom_idx,
                arena,
                orig_idx: atom_idx,
                state: State::Uninit,
                offset: Ix::new(0),
            }
        }
    }
    impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> Iterator for Edges<Ix, R> {
        type Item = EdgeReference<Ix>;

        #[instrument(level = "debug", name = "edges_next", skip(self), fields(mol_idx = self.mol_idx.index(), atom_idx = self.atom_idx.index()))]
        fn next(&mut self) -> Option<Self::Item> {
            loop {
                trace!(
                    mol_idx = self.mol_idx.index(),
                    atom_idx = self.atom_idx.index(),
                    "searching for atom"
                );
                if self.mol_idx == <Ix as IndexType>::max() {
                    return None;
                }
                let arena = self.arena.get_arena();
                match arena.frags[self.mol_idx.index()].repr {
                    MolRepr::Modify(ModdedMol { base, .. }) => self.mol_idx = base,
                    MolRepr::Atomic(ref b) => {
                        if !matches!(self.state, State::Atomic(_)) {
                            let i = b.nth(self.atom_idx.index())?;
                            trace!(id = i, "initializing walker at atomic level");
                            let w = arena.graph().neighbors(Ix::new(i).into()).detach();
                            self.state = State::Atomic(w);
                        }
                        let State::Atomic(w) = &mut self.state else {
                            unreachable!()
                        };

                        while let Some((e, n)) = w.next(arena.graph()) {
                            let Some(idx) = b.index(n.index()) else {
                                continue;
                            };

                            // get the correct offset index
                            let offset = idx + self.offset.index();
                            return Some(EdgeReference::with_weight(
                                self.orig_idx,
                                Ix::new(offset),
                                arena.graph()[e],
                            ));
                        }
                        return None;
                    }
                    MolRepr::Broken(ref b) => {
                        let next_frag = if let State::Broken(next) = &mut self.state {
                            if let Some(id) = next.pop() {
                                return Some(EdgeReference::with_weight(
                                    self.orig_idx,
                                    id,
                                    Bond::Single,
                                ));
                            }
                            true
                        } else {
                            false
                        };
                        debug!("generating state for broken molecule");
                        if next_frag {
                            for &i in &b.frags {
                                let size = arena.frags[i.index()].size();
                                if let Some(new) = self.atom_idx.index().checked_sub(size) {
                                    self.atom_idx = Ix::new(new);
                                    self.offset = Ix::new(self.offset.index() + size);
                                } else {
                                    self.mol_idx = i;
                                    break;
                                }
                            }
                            self.state = State::Uninit;
                        } else {
                            let offsets = b
                                .frags
                                .iter()
                                .scan(self.offset.index(), |count, idx| {
                                    let old = *count;
                                    *count += arena.frags[idx.index()].size();
                                    Some(old)
                                })
                                .collect::<SmallVec<_, 3>>();
                            let mut n = 0;
                            let mut i = self.atom_idx.index();
                            for &f in &b.frags {
                                if let Some(new) = i.checked_sub(arena.frags[f.index()].size()) {
                                    n += 1;
                                    i = new;
                                }
                            }
                            let mut it = b.bonds.iter().filter_map(|ifb| {
                                if ifb.an.index() == n && ifb.ai.index() == i {
                                    return Some(Ix::new(offsets[ifb.bn.index()] + ifb.bi.index()));
                                }
                                if ifb.bn.index() == n && ifb.bi.index() == i {
                                    return Some(Ix::new(offsets[ifb.an.index()] + ifb.ai.index()));
                                }
                                None
                            });
                            let ret = it.next();
                            let cont = it.collect::<SmallVec<_, 4>>();
                            if cont.is_empty() {
                                self.mol_idx = b.frags[n];
                                self.atom_idx = Ix::new(i);
                                self.offset = Ix::new(offsets[n]);
                                self.state = State::Uninit;
                            } else {
                                self.state = State::Broken(cont);
                            }
                            if let Some(ret) = ret {
                                return Some(EdgeReference::with_weight(
                                    self.orig_idx,
                                    ret,
                                    Bond::Single,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    pub struct Neighbors<Ix: IndexType + Ord, R>(pub Edges<Ix, R>);

    impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> Iterator for Neighbors<Ix, R> {
        type Item = NodeIndex<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            self.0.next().map(|n| n.target())
        }
    }

    pub struct EdgesDirected<Ix: IndexType + Ord, R>(pub Edges<Ix, R>, pub Direction);

    impl<Ix: IndexType + Ord, R: ArenaAccessor<Ix = Ix>> Iterator for EdgesDirected<Ix, R> {
        type Item = EdgeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            let idx = self.0.orig_idx;
            self.0.next().map(|e| {
                if (e.source().0 == idx) ^ (self.1 == Direction::Outgoing) {
                    e.rev()
                } else {
                    e
                }
            })
        }
    }
}
