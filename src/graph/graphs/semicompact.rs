use itertools::Itertools;
use petgraph::graph::{EdgeIndex, Graph, IndexType, NodeIndex};
use petgraph::{Direction, EdgeType};

/// A trait for types acting like `Option`. Mostly here for its special impl for `Atom`, which uses a scratch bit.
pub trait Optional {
    type Inner;

    fn is_some(&self) -> bool;
    /// Create a version of this type that's empty
    fn none() -> Self;
    /// Create a version of this type
    fn some(val: Self::Inner) -> Self;

    fn unwrap_ref(&self) -> &Self::Inner;

    fn unwrap_mut(&mut self) -> &mut Self::Inner;
}
impl Optional for crate::core::Atom {
    type Inner = Self;

    fn is_some(&self) -> bool {
        self.data.scratch() & (1 << 5) == 0
    }
    fn none() -> Self {
        Self::new_scratch(0, 1 << 5)
    }
    fn some(mut val: Self::Inner) -> Self {
        val.map_scratch(|s| s & ((1 << 5) - 1));
        val
    }
    fn unwrap_ref(&self) -> &Self::Inner {
        self
    }
    fn unwrap_mut(&mut self) -> &mut Self::Inner {
        self
    }
}
impl<T> Optional for Option<T> {
    type Inner = T;

    fn is_some(&self) -> bool {
        Option::is_some(self)
    }
    fn none() -> Self {
        None
    }
    fn some(val: T) -> Self {
        Some(val)
    }
    fn unwrap_ref(&self) -> &Self::Inner {
        self.as_ref().unwrap()
    }
    fn unwrap_mut(&mut self) -> &mut Self::Inner {
        self.as_mut().unwrap()
    }
}

/// A type that "allocates" runs of contiguous nodes.
/// Might make this a full graph type later, for now it's just a thin wrapper that you could mess up
#[derive(Debug, Clone)]
pub struct GraphNodeAlloc<N, E, Ty: EdgeType, Ix: IndexType> {
    /// The underlying graph we use for storage
    pub inner: Graph<N, E, Ty, Ix>,
    /// Holes we've created in the graph
    holes: [[Ix; 2]; 8],
}
impl<N, E, Ty: EdgeType, Ix: IndexType> GraphNodeAlloc<N, E, Ty, Ix> {
    pub fn new() -> Self {
        Self::with_capacity(0, 0)
    }

    /// Create a new semi-connected graph with a given capacity
    pub fn with_capacity(nodes: usize, edges: usize) -> Self {
        Self {
            inner: Graph::with_capacity(nodes, edges),
            holes: [[IndexType::max(); 2]; 8],
        }
    }

    /// Create a new semi-connected graph from an underlying graph, assuming that all elements are in the "filled" state.
    pub fn from_compact(graph: Graph<N, E, Ty, Ix>) -> Self {
        Self {
            inner: graph,
            holes: [[IndexType::max(); 2]; 8],
        }
    }

    pub fn clear(&mut self) {
        self.inner.clear();
        self.holes = [[IndexType::max(); 2]; 8];
    }
}
impl<N: Optional, E, Ty: EdgeType, Ix: IndexType> GraphNodeAlloc<N, E, Ty, Ix> {
    /// Allocate a contiguous range of nodes, returns the index of the first one.
    pub fn allocate_range(&mut self, size: usize, mut fill: impl FnMut() -> N::Inner) -> usize {
        if size == 0 {
            return 0;
        }
        let best = self
            .holes
            .iter()
            .enumerate()
            .filter_map(|(n, &r @ [from, to])| {
                let sz = to.index() - from.index();
                (sz >= size).then_some((n, r, sz))
            })
            .min_by_key(|x| size.abs_diff(x.2 * 2));
        if let Some((i, [from, _to], sz)) = best {
            let start = from.index();
            for i in start..(start + size) {
                let elem = self.inner.node_weight_mut(Ix::new(i).into()).unwrap();
                debug_assert!(!elem.is_some(), "overwrite of element detected!");
                *elem = N::some(fill());
            }
            if sz > size {
                self.holes[i][0] = Ix::new(start + size);
            } else {
                self.holes[i..].rotate_left(1);
                self.update_holes(1);
            }
            for i in start..(start + size) {
                self.inner[NodeIndex::new(i)] = N::some(fill());
            }
            start
        } else {
            let start = self.inner.node_count();
            self.inner.reserve_nodes(size);
            for _ in 0..size {
                self.inner.add_node(N::some(fill()));
            }
            start
        }
    }

    pub fn free_range(&mut self, start: usize, end: usize) {
        for i in start..end {
            self.inner[NodeIndex::new(i)] = N::none();
            self.detach_node(i);
        }
        let mut merge_to: Option<&mut Ix> = None;
        let mut rotate = None;
        for (i, [s, e]) in self.holes.iter_mut().enumerate() {
            if s.index() == end {
                if let Some(s2) = merge_to {
                    *s2 = *s;
                    rotate = Some(i);
                    break;
                } else {
                    *s = Ix::new(start);
                    merge_to = Some(e);
                }
            } else if e.index() == start {
                if let Some(e2) = merge_to {
                    *e2 = *e;
                    rotate = Some(i);
                    break;
                } else {
                    *e = Ix::new(end);
                    merge_to = Some(s);
                }
            }
        }
        if let Some(i) = rotate {
            self.holes[i..].rotate_left(1);
        }
        self.update_holes(1);
    }

    /// Remove all edges from a node.
    fn detach_node(&mut self, idx: usize) {
        let mut next = self.inner.raw_nodes()[idx].next_edge(Direction::Outgoing);
        while next != EdgeIndex::end() {
            next = self.inner.raw_nodes()[idx].next_edge(Direction::Outgoing);
            self.inner.remove_edge(next);
        }
        next = self.inner.raw_nodes()[idx].next_edge(Direction::Incoming);
        while next != EdgeIndex::end() {
            next = self.inner.raw_nodes()[idx].next_edge(Direction::Incoming);
            self.inner.remove_edge(next);
        }
    }

    /// Update our hole list, with last being the number of elements from the end that we need to check
    fn update_holes(&mut self, last: usize) {
        let hlen = self.holes.len();
        let start_ix = self.holes[hlen - last - 1][0].index();
        let empty = [IndexType::max(); 2];
        if start_ix >= self.inner.node_count() {
            self.holes[(hlen - last)..].fill(empty);
            return;
        }
        let chunks = self.inner.raw_nodes()[start_ix..]
            .iter()
            .enumerate()
            .chunk_by(|n| !n.1.weight.is_some());
        let mut iter = chunks.into_iter().filter_map(|(p, mut g)| {
            p.then(|| -> Option<_> {
                let start = g.next()?.0;
                let end = g.last()?.0;
                Some([start + start_ix, end + start_ix])
            })
            .flatten()
        });
        for hole in &mut self.holes[(hlen - last)..] {
            *hole = iter.next().map_or(empty, |arr| arr.map(Ix::new));
        }
    }
}

impl<N, E, Ty: EdgeType, Ix: IndexType> Default for GraphNodeAlloc<N, E, Ty, Ix> {
    fn default() -> Self {
        Self::new()
    }
}
