use crate::utils::bitset::BitSet;
use num_traits::PrimInt;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::fmt::{self, Binary, Debug, Formatter};

/// Iterate over connected graphs in a graph, returning their bits.
#[derive(Clone)]
pub struct ConnectedGraphIter<T = usize, const N: usize = 8, F = fn(usize) -> usize> {
    pub full: BitSet<T, N>,
    pub seen: BitSet<T, N>,
    pub cvt: F,
}

impl<T: PrimInt, const N: usize> ConnectedGraphIter<T, N, fn(usize) -> usize> {
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
        Self {
            full,
            seen,
            cvt: std::convert::identity,
        }
    }
    pub fn from_full(full: BitSet<T, N>) -> Self {
        Self {
            seen: BitSet::with_capacity(std::mem::size_of_val(full.as_slice()) * 8),
            full,
            cvt: std::convert::identity,
        }
    }
}
impl<T: PrimInt, const N: usize, F: FnMut(usize) -> usize> ConnectedGraphIter<T, N, F> {
    pub fn with_cvt<G: IntoNodeIdentifiers + NodeIndexable>(graph: G, cvt: F) -> Self {
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
        Self { full, seen, cvt }
    }
}

impl<T: Binary, const N: usize, F> Debug for ConnectedGraphIter<T, N, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("DisjointGraphIter")
            .field("full", &self.full)
            .field("seen", &self.seen)
            .finish()
    }
}

impl<G: IntoNeighbors + NodeIndexable, T: PrimInt, const N: usize, F: FnMut(usize) -> usize>
    Walker<G> for ConnectedGraphIter<T, N, F>
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
            out.set((self.cvt)(idx), true);
            stack.extend(graph.neighbors(id).filter(|&id| {
                self.full.get(graph.to_index(id)) && !self.seen.get(graph.to_index(id))
            }))
        }
        Some(out)
    }
}
