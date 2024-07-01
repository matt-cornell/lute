use petgraph::algo::Measure;
use petgraph::prelude::*;
use petgraph::unionfind::UnionFind;
use petgraph::visit::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::fmt::{self, Debug, Formatter};
use std::mem::MaybeUninit;

/// Taken from petgraph
#[derive(Debug, Clone, Copy)]
struct MinScored<K, T>(K, T);
impl<K: PartialOrd, T> PartialEq for MinScored<K, T> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl<K: PartialOrd, T> Eq for MinScored<K, T> {}
impl<K: PartialOrd, T> PartialOrd for MinScored<K, T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
#[allow(clippy::eq_op)] // we have to check for PartialEq
impl<K: PartialOrd, T> Ord for MinScored<K, T> {
    fn cmp(&self, other: &Self) -> Ordering {
        let a = &self.0;
        let b = &other.0;
        a.partial_cmp(b).unwrap_or_else(|| match (a != a, b != b) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => Ordering::Equal,
        })
    }
}

/// Find the number of cycles in a graph by counting the edges that aren't in a spanning tree.
pub fn num_cycles<G: NodeCompactIndexable + IntoEdgeReferences>(graph: G) -> usize {
    let mut union = UnionFind::new(graph.node_count());
    let mut out = 0;
    for edge in graph.edge_references() {
        let ai = graph.to_index(edge.source());
        let bi = graph.to_index(edge.target());
        if !union.union(ai, bi) {
            out += 1;
        }
    }
    out
}

type CycleCsr<K> = petgraph::csr::Csr<(), (K, usize), Undirected, usize>;

/// An iterator over the cycle basis of a graph.
/// Once created, it doesn't need to borrow from the graph at all, since it keeps its own copy of
/// the data.
#[derive(Clone)]
pub struct CycleBasis<K> {
    visited: <CycleCsr<K> as Visitable>::Map,
    tree: CycleCsr<K>,
    scores: Vec<(Option<K>, Vec<usize>)>,
    visit_next: BinaryHeap<MinScored<K, usize>>,
    cycle_edges: VecDeque<(usize, usize)>,
    depth: usize,
    result: Vec<usize>,
}
impl<K: Debug> Debug for CycleBasis<K> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_struct("CycleBasis")
            .field("scores", &self.scores)
            .field("visit_next", &self.visit_next)
            .field("cycle_edges", &self.cycle_edges)
            .field("depth", &self.depth)
            .finish_non_exhaustive()
    }
}
impl<K: Measure> CycleBasis<K> {
    /// Construct a new `CycleBasis` from a graph given a cost function.
    pub fn new<G: NodeCompactIndexable + IntoEdgeReferences, F: FnMut(G::EdgeRef) -> K>(
        graph: G,
        mut cost: F,
    ) -> Self {
        let nodes = graph.node_count();
        let mut edges = BinaryHeap::with_capacity(nodes);
        let mut union = UnionFind::new(nodes);
        edges.extend(graph.edge_references().map(|e| {
            let i = (graph.to_index(e.source()), graph.to_index(e.target()));
            MinScored(cost(e), i)
        }));
        let mut tree = CycleCsr::with_nodes(nodes);
        let mut cycle_edges = VecDeque::new();
        while let Some(MinScored(w, (ai, bi))) = edges.pop() {
            let in_tree = union.union(ai, bi);
            tree.add_edge(ai, bi, (w, if in_tree { 0 } else { cycle_edges.len() + 1 }));
            if !in_tree {
                cycle_edges.push_back((ai, bi));
            }
        }
        let visited = tree.visit_map();
        let scores = vec![(None, Vec::new()); nodes];
        let visit_next = BinaryHeap::new();
        CycleBasis {
            visited,
            tree,
            scores,
            visit_next,
            cycle_edges,
            depth: 0,
            result: Vec::new(),
        }
    }
    /// Construct a new `CycleBasis`, using `PartialOrd` node weights.
    pub fn new_ord<G: NodeCompactIndexable + IntoEdgeReferences<EdgeWeight = K>>(graph: G) -> Self
    where
        K: Copy,
    {
        Self::new(graph, |e| *e.weight())
    }
}
impl CycleBasis<usize> {
    /// Construct a `CycleBasis` that gives every edge a weight of 1.
    pub fn new_struct<G: NodeCompactIndexable + IntoEdgeReferences>(graph: G) -> Self {
        Self::new(graph, |_| 1usize)
    }
}
impl<K: Measure + Copy> CycleBasis<K> {
    /// The function that does all of the work. Takes start and end nodes, along with the allowed
    /// cycle edges to use.
    fn compute(&mut self, ai: usize, bi: usize, n: usize) -> K {
        self.tree.reset_map(&mut self.visited);
        for (s, p) in &mut self.scores {
            *s = None;
            p.clear();
        }
        self.scores[ai].0 = Some(K::default());
        self.scores[ai].1.push(ai);
        self.visit_next.clear();
        self.visit_next.push(MinScored(K::default(), ai));
        while let Some(MinScored(score, node)) = self.visit_next.pop() {
            if self.visited.is_visited(&node) {
                continue;
            }
            if node == bi {
                break;
            }
            for edge in self.tree.edges(node) {
                let next = edge.target();
                if self.visited.is_visited(&next) {
                    continue;
                }
                let &(weight, d) = edge.weight();
                if d > n {
                    continue;
                }
                let next_score = score + weight;
                let Ok([(_, new_path), (entry, path)]) = self.scores.get_many_mut([node, next])
                else {
                    continue;
                };
                if let Some(old_score) = entry {
                    if next_score < *old_score {
                        *old_score = next_score;
                        self.visit_next.push(MinScored(next_score, next));

                        path.clear();

                        {
                            path.reserve(new_path.len());
                            MaybeUninit::copy_from_slice(
                                &mut path.spare_capacity_mut()[..new_path.len()],
                                new_path,
                            );
                            unsafe { path.set_len(new_path.len()) }
                        }

                        path.push(next);
                    }
                } else {
                    *entry = Some(next_score);
                    self.visit_next.push(MinScored(next_score, next));

                    {
                        path.reserve(new_path.len());
                        MaybeUninit::copy_from_slice(
                            &mut path.spare_capacity_mut()[..new_path.len()],
                            new_path,
                        );
                        unsafe { path.set_len(new_path.len()) }
                    }

                    path.push(next);
                }
            }
            self.visited.visit(node);
        }
        std::mem::swap(&mut self.scores[bi].1, &mut self.result);
        self.scores[bi]
            .0
            .take()
            .expect("a path should've been found!")
    }
    /// Compute the next cycle, but hold on to the storage to allow it to be reclaimed.
    pub fn step(&mut self) -> Option<K> {
        let n = self.depth;
        let (ai, bi) = self.cycle_edges.pop_front()?;
        self.depth += 1;
        Some(self.compute(ai, bi, n))
    }
    /// Peek at the nth cycle.
    pub fn peek_nth(&mut self, mut n: usize) -> Option<K> {
        let &(ai, bi) = self.cycle_edges.get(n)?;
        n += self.depth;
        Some(self.compute(ai, bi, n))
    }
    /// Get the current cycle, if computed with `step`. If used as an iterator, this returns an
    /// empty slice.
    pub fn current(&self) -> &[usize] {
        &self.result
    }
    /// See `current`.
    pub fn current_mut(&mut self) -> &mut Vec<usize> {
        &mut self.result
    }
    pub fn len(&self) -> usize {
        self.cycle_edges.len()
    }
    pub fn is_empty(&self) -> bool {
        self.cycle_edges.is_empty()
    }
}
impl<K: Measure + Copy> Iterator for CycleBasis<K> {
    type Item = (K, Vec<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        let n = self.depth;
        let (ai, bi) = self.cycle_edges.pop_front()?;
        self.depth += 1;
        let weight = self.compute(ai, bi, n);
        Some((weight, std::mem::take(&mut self.result)))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.cycle_edges.len();
        (len, Some(len))
    }
}
impl<K: Measure + Copy> DoubleEndedIterator for CycleBasis<K> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let (ai, bi) = self.cycle_edges.pop_back()?;
        let n = self.depth + self.cycle_edges.len();
        let weight = self.compute(ai, bi, n);
        Some((weight, std::mem::take(&mut self.result)))
    }
}
impl<K: Measure + Copy> ExactSizeIterator for CycleBasis<K> {}
