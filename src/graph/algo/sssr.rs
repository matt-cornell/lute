use petgraph::visit::*;
use std::collections::hash_map::Entry;
use std::collections::VecDeque;
use std::hash::Hash;
use ahash::*;

#[derive(Default, Clone)]
pub struct RingsIterator<N> {
    pub stack: Vec<N>,
    /// Map of node to predecessors
    pub paths: HashMap<N, (N, usize)>,
    /// We can find more cycles in a single iteration, store the leftovers here
    pub to_ret: Vec<VecDeque<N>>,
}
impl<N> RingsIterator<N> {
    pub fn new<G: IntoNodeIdentifiers<NodeId = N>>(graph: G) -> Self {
        let stack = graph.node_identifiers().next();
        Self {
            stack: stack.into_iter().collect(),
            paths: HashMap::new(),
            to_ret: Vec::new(),
        }
    }
}
impl<N: Hash + Eq + Copy, G: IntoNeighbors<NodeId = N>> Walker<G> for RingsIterator<N> {
    type Item = VecDeque<N>;
    fn walk_next(&mut self, context: G) -> Option<VecDeque<N>> {
        if let Some(ret) = self.to_ret.pop() {
            return Some(ret);
        }
        while let Some(n) = self.stack.pop() {
            let mut path = VecDeque::new();
            let d = {
                if let Some(&(mut p, n)) = self.paths.get(&n) {
                    if n > 0 {
                        path.reserve_exact(n);
                        path.push_back(p);
                        while let Some(&(x, n)) = self.paths.get(&p) {
                            if n == 0 {
                                break;
                            }
                            path.push_back(x);
                            p = x;
                        }
                        n + 1
                    } else {
                        1
                    }
                } else {
                    1
                }
            };
            'succs: for succ in context.neighbors(n) {
                match self.paths.entry(succ) {
                    Entry::Occupied(e) => {
                        let mut pred = e.into_mut().0;
                        let mut ret = path.clone();
                        while !path.contains(&pred) {
                            ret.push_front(pred);
                            let Some(&(p, n)) = self.paths.get(&pred) else {
                                continue 'succs;
                            };
                            if n == 0 {
                                continue 'succs;
                            }
                            pred = p;
                        }
                        self.to_ret.push(ret);
                    }
                    Entry::Vacant(e) => {
                        e.insert((n, d));
                    }
                }
            }
            if let Some(ret) = self.to_ret.pop() {
                return Some(ret);
            }
        }
        None
    }
}
