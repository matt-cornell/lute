pub use super::color::*;
use crate::core::*;
use crate::graph::algo::CycleBasis;
use ahash::AHashSet;
use petgraph::visit::*;
use petgraph::Undirected;
use std::f64::consts::{PI, TAU};
use std::fmt::{self, Display, Formatter};

pub fn fmt_as_svg<G>(graph: G) -> String
where
    G: NodeCompactIndexable + IntoEdgeReferences + Visitable,
{
    let cycles = CycleBasis::new_struct(graph)
        .map(|x| x.1)
        .collect::<Vec<_>>();
    let mut edges = AHashSet::with_capacity(cycles.len());
    for cycle in &cycles {
        for &[low, hi] in cycle.array_windows() {
            edges.insert(std::cmp::minmax(low, hi));
        }
        if let Some(&last) = cycle.last() {
            let first = cycle[0];
            if first != last {
                edges.insert(std::cmp::minmax(first, last));
            }
        }
    }
    todo!()
}
