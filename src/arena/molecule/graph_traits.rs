use super::*;
use petgraph::{visit::*, Direction};
use std::iter::Map;
use std::ops::Range;

impl<Ix: Copy + PartialEq, R> GraphBase for Molecule<Ix, R> {
    type NodeId = NodeIndex<Ix>;
    type EdgeId = EdgeIndex<Ix>;
}

impl<Ix: Copy + PartialEq, R: Copy> GraphRef for Molecule<Ix, R> {}

impl<Ix: Copy + PartialEq, R> GraphProp for Molecule<Ix, R> {
    type EdgeType = petgraph::Undirected;
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> NodeCount for Molecule<Ix, R> {
    fn node_count(&self) -> usize {
        self.arena.get_arena().parts[self.index.index()].1.index()
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> NodeIndexable for Molecule<Ix, R> {
    fn node_bound(&self) -> usize {
        self.arena.get_arena().parts[self.index.index()].1.index()
    }
    fn to_index(&self, a: NodeIndex<Ix>) -> usize {
        a.0.index()
    }
    fn from_index(&self, i: usize) -> NodeIndex<Ix> {
        NodeIndex(Ix::new(i))
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> NodeCompactIndexable for Molecule<Ix, R> {}

impl<Ix: IndexType, R> Data for Molecule<Ix, R> {
    type NodeWeight = Atom;
    type EdgeWeight = Bond;
}

impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeIdentifiers for Molecule<Ix, R> {
    type NodeIdentifiers = Map<Range<usize>, fn(usize) -> NodeIndex<Ix>>;

    fn node_identifiers(self) -> Self::NodeIdentifiers {
        (0..self.node_count()).map(|ix| NodeIndex(Ix::new(ix)))
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNodeReferences for Molecule<Ix, R> {
    type NodeRef = NodeReference<Ix>;
    type NodeReferences = iter::NodeReferences<Ix, R>;

    fn node_references(self) -> Self::NodeReferences {
        iter::NodeReferences {
            range: 0..self.node_bound(),
            mol_idx: self.index,
            arena: self.arena,
        }
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdgeReferences for Molecule<Ix, R> {
    type EdgeRef = EdgeReference<Ix>;
    type EdgeReferences = iter::EdgeReferences<Ix, R>;

    fn edge_references(self) -> Self::EdgeReferences {
        iter::EdgeReferences::new(self.index, self.arena)
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdges for Molecule<Ix, R> {
    type Edges = iter::Edges<Ix, R>;

    fn edges(self, a: Self::NodeId) -> Self::Edges {
        iter::Edges::new(self.index, a.0, self.arena)
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoEdgesDirected for Molecule<Ix, R> {
    type EdgesDirected = iter::EdgesDirected<Ix, R>;

    fn edges_directed(self, a: Self::NodeId, dir: Direction) -> Self::EdgesDirected {
        iter::EdgesDirected(iter::Edges::new(self.index, a.0, self.arena), dir)
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNeighbors for Molecule<Ix, R> {
    type Neighbors = iter::Neighbors<Ix, R>;

    fn neighbors(self, a: Self::NodeId) -> Self::Neighbors {
        iter::Neighbors(iter::Edges::new(self.index, a.0, self.arena))
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> IntoNeighborsDirected for Molecule<Ix, R> {
    type NeighborsDirected = iter::Neighbors<Ix, R>;

    fn neighbors_directed(
        self,
        n: Self::NodeId,
        _: petgraph::prelude::Direction,
    ) -> Self::NeighborsDirected {
        iter::Neighbors(iter::Edges::new(self.index, n.0, self.arena))
    }
}
impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix> + Copy> GetAdjacencyMatrix for Molecule<Ix, R> {
    type AdjMatrix = Self;

    fn adjacency_matrix(&self) -> Self::AdjMatrix {
        *self
    }
    fn is_adjacent(&self, matrix: &Self, a: Self::NodeId, b: Self::NodeId) -> bool {
        matrix.get_bond(EdgeIndex::new(a.0, b.0)).is_some()
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
    impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for NodeReferences<Ix, R> {
        type Item = NodeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            self.range
                .next()
                .map(|i| NodeReference::new(self.mol_idx, NodeIndex(Ix::new(i)), self.arena))
        }
    }

    #[allow(dead_code)]
    #[derive(Clone)]
    pub struct EdgeReferences<Ix, R> {
        stack: SmallVec<Ix, 3>,
        arena: R,
    }
    impl<Ix, R> EdgeReferences<Ix, R> {
        pub fn new(mol_idx: Ix, arena: R) -> Self {
            Self {
                stack: smallvec![mol_idx],
                arena,
            }
        }
    }
    impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for EdgeReferences<Ix, R> {
        type Item = EdgeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            todo!()
        }
    }

    #[derive(Clone)]
    enum State<Ix: IndexType> {
        Uninit,
        Broken(SmallVec<(Ix, Bond), 4>),
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

    #[derive(Debug, Clone)]
    pub struct Edges<Ix: IndexType, R> {
        orig_idx: Ix,
        mol_idx: Ix,
        atom_idx: Ix,
        arena: R,
        state: State<Ix>,
        offset: Ix,
        skip: BTreeSet<Ix>,
    }
    impl<Ix: IndexType, R> Edges<Ix, R> {
        pub fn new(mol_idx: Ix, atom_idx: Ix, arena: R) -> Self {
            Self {
                mol_idx,
                atom_idx,
                arena,
                orig_idx: atom_idx,
                state: State::Uninit,
                offset: Ix::new(0),
                skip: BTreeSet::new(),
            }
        }
    }
    impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for Edges<Ix, R> {
        type Item = EdgeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            loop {
                if self.mol_idx == <Ix as IndexType>::max() {
                    return None;
                }
                let arena = self.arena.get_arena();
                match arena.parts[self.mol_idx.index()].0 {
                    MolRepr::Redirect(ix) | MolRepr::Modify(ModdedMol { base: ix, .. }) => {
                        self.mol_idx = ix
                    }
                    MolRepr::Atomic(ref b) => {
                        if !matches!(self.state, State::Atomic(_)) {
                            let w = arena.graph().neighbors(self.atom_idx.into()).detach();
                            self.state = State::Atomic(w);
                        }
                        let State::Atomic(w) = &mut self.state else {
                            unreachable!()
                        };

                        while let Some((e, n)) = w.next(arena.graph()) {
                            if !b.get(n.index()) {
                                continue;
                            }

                            // get the correct offset index
                            let offset = self.skip.range(..Ix::new(n.index())).count().index()
                                + self.offset.index();

                            return Some(EdgeReference::with_weight(
                                EdgeIndex::new(self.orig_idx, Ix::new(offset)),
                                arena.graph()[e],
                            ));
                        }
                        return None;
                    }
                    MolRepr::Broken(ref b) => {
                        if let State::Broken(next) = &mut self.state {
                            if let Some((id, w)) = next.pop() {
                                return Some(EdgeReference::with_weight(
                                    EdgeIndex::new(self.orig_idx, id),
                                    w,
                                ));
                            }
                        }
                        let empty = BTreeSet::new();
                        let mut skips = HybridMap::<Ix, BTreeSet<Ix>, 4>::new();
                        let mut ibs = HybridMap::<_, SmallVec<_, 4>, 4>::new();
                        let mut idx = self.atom_idx.index();
                        for i in &b.bonds {
                            let (ix, a) =
                                arena.molecule(i.an).get_atom_idx(NodeIndex(i.ai)).unwrap();
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
                            let (ix, a) =
                                arena.molecule(i.bn).get_atom_idx(NodeIndex(i.bi)).unwrap();
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
                            let weight = weight.unwrap_or(Bond::Single);
                            eprintln!("bond: {i:#?}");
                            if let Some(s) = ibs.get_mut(&(i.an, i.ai)) {
                                s.push((i.bn, i.bi, weight));
                            } else {
                                ibs.insert((i.an, i.ai), smallvec![(i.bn, i.bi, weight)]);
                            }
                            if let Some(s) = ibs.get_mut(&(i.bn, i.bi)) {
                                s.push((i.an, i.ai, weight));
                            } else {
                                ibs.insert((i.bn, i.bi), smallvec![(i.an, i.ai, weight)]);
                            }
                        }
                        let mut new_idx = None;
                        let mut f = SmallVec::<usize, 4>::with_capacity(b.frags.len());
                        let mut counter = 0;
                        let mut up_idx = true;
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
                            f.push(counter);
                            counter += s;
                            if up_idx {
                                if idx > s {
                                    idx -= s;
                                } else {
                                    idx += ic;
                                    new_idx = Some(*p);
                                    up_idx = false;
                                }
                            }
                        }
                        self.mol_idx = new_idx.unwrap_or(<Ix as IndexType>::max());
                        self.atom_idx = Ix::new(idx);
                        eprintln!("ibs: {ibs:#?}");
                        let next = ibs.get(&(self.mol_idx, self.atom_idx)).map_or_else(
                            SmallVec::new,
                            |v| {
                                v.iter()
                                    .map(|&(n, i, b)| {
                                        let id = f[n.index()]
                                            + i.index()
                                            + skips.get(&n).map_or(0, |s| s.range(..i).count());
                                        (Ix::new(id), b)
                                    })
                                    .collect()
                            },
                        );
                        self.state = State::Broken(next);
                    }
                }
            }
        }
    }

    pub struct Neighbors<Ix: IndexType, R>(pub Edges<Ix, R>);

    impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for Neighbors<Ix, R> {
        type Item = NodeIndex<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            self.0.next().map(|n| {
                if n.source().0 == self.0.atom_idx {
                    n.target()
                } else {
                    n.source()
                }
            })
        }
    }

    pub struct EdgesDirected<Ix: IndexType, R>(pub Edges<Ix, R>, pub Direction);

    impl<Ix: IndexType, R: ArenaAccessor<Ix = Ix>> Iterator for EdgesDirected<Ix, R> {
        type Item = EdgeReference<Ix>;

        fn next(&mut self) -> Option<Self::Item> {
            let idx = self.0.orig_idx;
            self.0
                .find(|i| (self.1 == Direction::Incoming) ^ (i.source().0 == idx))
        }
    }
}
