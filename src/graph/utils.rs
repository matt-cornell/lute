use crate::utils::bitset::BitSet;
use num_traits::PrimInt;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::fmt::{self, Binary, Debug, Formatter};

/// Iterate over disjoint graphs in a graph, returning their bits.
#[derive(Clone)]
pub struct DisjointGraphIter<T = usize, const N: usize = 8> {
    pub full: BitSet<T, N>,
    pub seen: BitSet<T, N>,
}

impl<T: PrimInt, const N: usize> DisjointGraphIter<T, N> {
    pub fn new<G: IntoNodeIdentifiers + NodeIndexable>(graph: G) -> Self {
        let mut full = BitSet::with_capacity(graph.node_bound());
        let mut max = 0;
        for id in graph.node_identifiers() {
            let idx = graph.to_index(id);
            full.set(idx, true);
            if idx > max {
                max = idx;
            }
        }
        let seen = BitSet::with_capacity(max);
        Self { full, seen }
    }
}

impl<T: Binary, const N: usize> Debug for DisjointGraphIter<T, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("DisjointGraphIter")
            .field("full", &self.full)
            .field("seen", &self.seen)
            .finish()
    }
}

impl<G: IntoNeighbors + NodeIndexable, T: PrimInt, const N: usize> Walker<G>
    for DisjointGraphIter<T, N>
{
    type Item = BitSet<T, N>;
    fn walk_next(&mut self, graph: G) -> Option<BitSet<T, N>> {
        let mut start = None;
        let bits = T::zero().count_zeros() as usize;
        for (n, (&f, &s)) in self
            .full
            .as_slice()
            .iter()
            .zip(self.seen.as_slice())
            .enumerate()
        {
            let bit = (f & !s).trailing_zeros() as usize;
            if bit < bits {
                start = Some(n * bits + bit);
                break;
            }
        }
        let mut stack: SmallVec<G::NodeId, 8> = smallvec![graph.from_index(start?)];
        let mut out = BitSet::with_capacity(self.seen.as_slice().len() * bits);
        while let Some(id) = stack.pop() {
            let idx = graph.to_index(id);
            self.seen.set(idx, true);
            out.set(idx, true);
            for n in graph.neighbors(id) {
                if !self.seen.get(graph.to_index(n)) {
                    stack.push(n);
                }
            }
        }
        Some(out)
    }
}
