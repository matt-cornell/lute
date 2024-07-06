use crate::atom_info::ATOM_DATA;
use crate::core::*;
use ahash::AHashMap;
use petgraph::visit::*;
use std::fmt::Write;
use std::hash::Hash;

fn bond2str(bond: Bond) -> &'static str {
    match bond {
        Bond::Non => ".",
        Bond::Single => "",
        Bond::Double | Bond::DoubleE | Bond::DoubleZ => "=",
        Bond::Triple => "#",
        Bond::Quad => "$",
        Bond::Left => "/",
        Bond::Right => "\\",
        Bond::Aromatic => ":",
        _ => "!",
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SmilesConfig {
    pub isomers: bool,
    pub isotopes: bool,
    pub canon: bool,
}
impl SmilesConfig {
    /// Default config, gives pretty output
    pub const fn new() -> Self {
        Self {
            isomers: true,
            isotopes: true,
            canon: true,
        }
    }
    /// Reuse rings and don't canonicalize, better for serialization.
    pub const fn fast_roundtrip() -> Self {
        Self {
            isomers: true,
            isotopes: true,
            canon: false,
        }
    }
}

pub fn generate_smiles<G>(graph: G, cfg: SmilesConfig) -> String
where
    G: Data<NodeWeight = Atom, EdgeWeight = Bond> + IntoNodeReferences + IntoEdges + NodeCount,
    G::NodeId: Hash + Eq + std::fmt::Debug,
{
    if graph.node_count() == 0 {
        return String::new();
    } // guarantee order is non-empty
    let mut out = String::with_capacity(graph.node_count());
    let mut order = Vec::new();
    if cfg.canon {
        let primes = primal::Sieve::new(
            primal::estimate_nth_prime(graph.node_count() as u64 + 1).1 as usize,
        );
        order.extend(graph.node_references().map(|n| {
            let atom = *n.weight();
            let data = atom.data;
            let protons = atom.protons as usize; // 8 bits
            let hydro = data.hydrogen() as usize; // 4 bits
            let non_h = (data.single() + data.other() + data.unknown()) as usize; // 6 bits
            let count = graph
                .edges(n.id())
                .map(|e| e.weight().bond_count())
                .sum::<f32>() as usize; // 6 bits
            let charge = (atom.charge.signum() + 1) as usize; // 2 bits
            let weight = hydro | (charge << 6) | (protons << 8) | (count << 16) | (non_h << 22);
            // 26 bits in total, but we use a u64 in case other weights get bigger
            (atom, n.id(), weight)
        }));
        order.sort_unstable_by_key(|n| n.2);
        let mut old = Vec::new();
        let mut scratch = Vec::new();
        let mut neighbors = Vec::new();
        let mut doubled = false;
        loop {
            {
                let mut it = primes.primes_from(0);
                let mut last = None;
                let mut prime = it.next().unwrap();
                for weight in &mut order {
                    if last != Some(weight.2) {
                        last = Some(weight.2);
                        prime = it.next().unwrap();
                    }
                    weight.2 = prime;
                }
            }
            scratch.clone_from(&order);
            for (_, id, weight) in &mut order {
                neighbors.clear();
                neighbors.extend(graph.neighbors(*id));
                *weight = scratch
                    .iter()
                    .filter_map(|x| neighbors.contains(&x.1).then_some(x.2))
                    .product();
            }
            order.sort_by_key(|n| n.2);
            if order.iter().map(|i| &i.1).eq(&old) {
                if doubled {
                    break;
                }
                doubled = true;
                let mut last: Option<&mut usize> = None;
                let mut done = true;
                let mut decr = true;
                for i in &mut order {
                    i.2 *= 2;
                    if let Some(last) = &mut last {
                        if decr && **last == i.2 {
                            **last -= 1;
                            decr = false;
                            done = false;
                        } else {
                            decr = true;
                        }
                        *last = &mut i.2;
                    } else {
                        last = Some(&mut i.2);
                    }
                }
                if done {
                    break;
                }
            } else {
                doubled = false;
            }
            old.clear();
            old.extend(order.iter().map(|i| i.1));
        }
    } else {
        order.extend(graph.node_references().map(|i| (*i.weight(), i.id(), 0)))
    }
    let lookup = order
        .iter()
        .enumerate()
        .map(|(n, o)| (o.1, n))
        .collect::<AHashMap<_, _>>();
    let mut visited = AHashMap::<G::NodeId, _>::with_capacity(order.len());
    let mut last_idx = 0;
    let mut queue = Vec::new();
    let mut counter = 0usize;
    let mut buf = String::new();
    let mut edge_buf = Vec::new();
    loop {
        if let Some((node, prev, edge, branch)) = queue.pop() {
            if let Some(node) = node {
                if let Some(idx) = visited.get_mut(&node) {
                    counter += 1;
                    if let Some(ch) = u32::try_from(counter)
                        .ok()
                        .and_then(|d| char::from_digit(d, 10))
                    {
                        out.insert(*idx, ch);
                        out.push_str(bond2str(edge));
                        out.push(ch);
                        *idx += 1;
                    } else {
                        buf.clear();
                        let _ = write!(buf, "%{counter}");
                        out.insert_str(*idx, &buf);
                        out.push_str(bond2str(edge));
                        out.push_str(&buf);
                        *idx += buf.len();
                    }
                } else {
                    if branch {
                        out.push('(');
                        queue.push((None, None, Bond::Non, false));
                    }
                    if !(edge == Bond::Non && out.is_empty()) {
                        out.push_str(bond2str(edge));
                    }
                    let idx = lookup[&node];
                    let atom = order[idx].0;
                    edge_buf.clear();
                    edge_buf.extend(graph.edges(node));
                    let bond_count = (atom.data.unknown() + atom.data.hydrogen()) as isize
                        + edge_buf
                            .iter()
                            .map(|i| i.weight().bond_count())
                            .sum::<f32>() as isize
                        + atom.charge as isize;
                    let ex_bonds = match atom.protons {
                        x @ 6..=9 => Some((10 - (x as i8) + atom.charge) as isize),
                        x @ 14..=17 => Some((18 - (x as i8) + atom.charge) as isize),
                        35 | 53 => Some(1),
                        _ => None,
                    };
                    let force = ex_bonds != Some(bond_count)
                        || (atom.protons > 0 && atom.isotope > 0 && cfg.isotopes)
                        || (atom.data.chirality().is_chiral() && cfg.isomers);
                    'print_atom: {
                        if !force {
                            match atom.protons {
                                0 => match atom.isotope {
                                    0xFFFF => {
                                        out.push('*');
                                        break 'print_atom;
                                    }
                                    0xFFFE => {
                                        out.push('A');
                                        break 'print_atom;
                                    }
                                    0xFFFC => {
                                        out.push('Q');
                                        break 'print_atom;
                                    }
                                    0x0100 => {
                                        out.push('X');
                                        break 'print_atom;
                                    }
                                    0x042B => {
                                        out.push('M');
                                        break 'print_atom;
                                    }
                                    _ => {}
                                },
                                5 => {
                                    out.push('B');
                                    break 'print_atom;
                                }
                                6 => {
                                    out.push('C');
                                    break 'print_atom;
                                }
                                7 => {
                                    out.push('N');
                                    break 'print_atom;
                                }
                                8 => {
                                    out.push('O');
                                    break 'print_atom;
                                }
                                9 => {
                                    out.push('F');
                                    break 'print_atom;
                                }
                                15 => {
                                    out.push('P');
                                    break 'print_atom;
                                }
                                16 => {
                                    out.push('S');
                                    break 'print_atom;
                                }
                                17 => {
                                    out.push_str("Cl");
                                    break 'print_atom;
                                }
                                35 => {
                                    out.push_str("Br");
                                    break 'print_atom;
                                }
                                53 => {
                                    out.push('I');
                                    break 'print_atom;
                                }
                                _ => {}
                            }
                        }
                        out.push('[');
                        if atom.isotope > 0 && cfg.isotopes {
                            let _ = write!(out, "{}", atom.isotope);
                        }
                        out.push_str(ATOM_DATA[atom.protons as usize].sym);
                        if ex_bonds != Some(bond_count) {
                            out.push('H');
                            let h = atom.data.hydrogen();
                            if h != 1 {
                                let _ = write!(out, "{h}");
                            }
                        }
                        match atom.charge {
                            0 => {}
                            1 => out.push('+'),
                            -1 => out.push('-'),
                            n => {
                                let _ = write!(out, "{n:+}");
                            }
                        }
                        out.push(']');
                    }
                    out.extend(std::iter::repeat('&').take(atom.data.unknown() as usize));
                    visited.insert(node, out.len());
                    let idx = queue.len();
                    queue.extend(
                        graph
                            .edges(node)
                            .filter(|e| prev != Some(e.target()))
                            .map(|e| (Some(e.target()), Some(node), *e.weight(), true)),
                    );
                    if cfg.canon {
                        queue[idx..].sort_by_cached_key(|i| usize::MAX - lookup[&i.0.unwrap()]);
                        for (n, i) in queue[idx..].iter().enumerate() {
                            if i.2 != Bond::Single {
                                queue.swap(idx, idx + n);
                                break;
                            }
                        }
                    }
                    if let Some(i) = queue.get_mut(idx) {
                        i.3 = false;
                    }
                }
            } else {
                out.push(')');
            }
        } else if cfg.canon {
            if let Some((n, next)) = order[last_idx..]
                .iter()
                .enumerate()
                .find(|i| !visited.contains_key(&i.1 .1))
            {
                last_idx += n + 1;
                queue.push((Some(next.1), None, Bond::Non, false));
                if !out.is_empty() {
                    out.push('.');
                }
            } else {
                break;
            }
        } else {
            if let Some((n, next)) = graph
                .node_identifiers()
                .skip(last_idx)
                .enumerate()
                .find(|i| !visited.contains_key(&i.1))
            {
                last_idx += n + 1;
                queue.push((Some(next), None, Bond::Non, false));
            } else {
                break;
            }
        }
    }
    out
}
