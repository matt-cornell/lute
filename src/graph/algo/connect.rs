use crate::utils::bitset::BitSet;
use num_traits::PrimInt;
use petgraph::visit::*;
use smallvec::{smallvec, SmallVec};
use std::fmt::{self, Binary, Debug, Formatter};
use std::num::NonZeroUsize;

/// Iterate over connected graphs in a graph, returning their bits.
#[derive(Clone)]
pub struct ConnectedGraphIter<T = usize, const N: usize = 8, F = fn(usize) -> usize> {
    pub full: BitSet<T, N>,
    pub cvt: F,
}

impl<T: PrimInt, const N: usize> ConnectedGraphIter<T, N, fn(usize) -> usize> {
    pub fn new<G: IntoNodeIdentifiers + NodeIndexable>(graph: G) -> Self {
        let mut full = BitSet::new();
        let mut max = 0;
        for id in graph.node_identifiers() {
            let idx = graph.to_index(id);
            full.set(idx, true);
            if idx > max {
                max = idx;
            }
        }
        Self {
            full,
            cvt: std::convert::identity,
        }
    }
    pub fn from_full(full: BitSet<T, N>) -> Self {
        Self {
            full,
            cvt: std::convert::identity,
        }
    }
}
impl<T: PrimInt, const N: usize, F: FnMut(usize) -> usize> ConnectedGraphIter<T, N, F> {
    pub fn with_cvt<G: IntoNodeIdentifiers + NodeIndexable>(graph: G, cvt: F) -> Self {
        let mut full = BitSet::new();
        let mut max = 0;
        for id in graph.node_identifiers() {
            let idx = graph.to_index(id);
            full.set(idx, true);
            if idx > max {
                max = idx;
            }
        }
        Self { full, cvt }
    }
    pub fn step<G: IntoNeighbors + NodeIndexable, U: PrimInt, const O: usize>(
        &mut self,
        graph: G,
        out: &mut BitSet<U, O>,
    ) -> Option<NonZeroUsize> {
        let start = self.full.nth(0)?;
        let start = graph.from_index(start);
        out.clear();
        let mut stack: SmallVec<G::NodeId, 8> = smallvec![start];
        let mut count = 0;
        while let Some(id) = stack.pop() {
            count += 1;
            let idx = graph.to_index(id);
            self.full.set(idx, false);
            let cvt = (self.cvt)(idx);
            out.set(cvt, true);
            stack.extend(graph.neighbors(id).filter(|&id| {
                let idx = graph.to_index(id);
                self.full.get(idx)
            }))
        }
        NonZeroUsize::new(count)
    }
}

impl<T: Binary, const N: usize, F> Debug for ConnectedGraphIter<T, N, F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("ConnectedGraphIter")
            .field("full", &self.full)
            .finish()
    }
}

impl<G: IntoNeighbors + NodeIndexable, T: PrimInt, const N: usize, F: FnMut(usize) -> usize>
    Walker<G> for ConnectedGraphIter<T, N, F>
{
    type Item = BitSet<T, N>;
    fn walk_next(&mut self, graph: G) -> Option<BitSet<T, N>> {
        let start = self.full.nth(0);
        let start = graph.from_index(start?);
        let mut stack: SmallVec<G::NodeId, 8> = smallvec![start];
        let mut out = BitSet::new();
        while let Some(id) = stack.pop() {
            let idx = graph.to_index(id);
            self.full.set(idx, false);
            let cvt = (self.cvt)(idx);
            out.set(cvt, true);
            stack.extend(graph.neighbors(id).filter(|&id| {
                let idx = graph.to_index(id);
                self.full.get(idx)
            }))
        }
        Some(out)
    }
}
