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
    let mut visited: AHashMap<_, (_, Vec<G::NodeId>)> = AHashMap::with_capacity(order.len());
    let mut last_idx = 0;
    let mut queue: Vec<(Option<G::NodeId>, Option<(_, _)>, _, _)> = Vec::new();
    let mut counter = 0usize;
    let mut buf = String::new();
    let mut edge_buf = Vec::new();
    loop {
        if let Some((node, prev, edge, branch)) = queue.pop() {
            if let Some(node) = node {
                if let Some(&(idx, ref prevs)) = visited.get(&node) {
                    if prevs.contains(&prev.unwrap().0) {
                        match out.pop() {
                            None => {}
                            Some(')') => {
                                let mut depth = 1;
                                for (n, i) in out.bytes().rev().enumerate() {
                                    match i {
                                        b')' => depth += 1,
                                        b'(' => {
                                            depth -= 1;
                                            if depth == 0 {
                                                out.remove(out.len() - n - 1);
                                                break;
                                            }
                                        }
                                        _ => {}
                                    }
                                }
                            }
                            Some(c) => out.push(c),
                        }
                        continue;
                    }
                    counter += 1;
                    let len = if let Some(ch) = u32::try_from(counter)
                        .ok()
                        .and_then(|d| char::from_digit(d, 10))
                    {
                        out.insert(idx, ch);
                        out.push_str(bond2str(edge));
                        out.push(ch);
                        1
                    } else {
                        buf.clear();
                        let _ = write!(buf, "%{counter}");
                        out.insert_str(idx, &buf);
                        out.push_str(bond2str(edge));
                        out.push_str(&buf);
                        buf.len()
                    };
                    for (n, (i, prevs)) in &mut visited {
                        if *i >= idx {
                            *i += len;
                        }
                        if *n == prev.unwrap().0 {
                            prevs.push(node);
                        }
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
                    let aromatic = ([5, 6, 7, 8, 15, 16].contains(&atom.protons)
                        || (atom.protons == 0 && [0xFFFE, 0xFFFC].contains(&atom.isotope)))
                        && {
                            let mut it = edge_buf.iter().filter(|e| *e.weight() == Bond::Aromatic);
                            it.next().is_some() && it.next().is_some()
                        };
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
                    let needs_h = ex_bonds != Some(bond_count)
                        && ex_bonds
                            .map_or(true, |e| e > bond_count - atom.data.hydrogen() as isize);
                    let force = needs_h
                        || (atom.protons > 0 && atom.isotope > 0 && cfg.isotopes)
                        || (atom.data.chirality().is_chiral() && cfg.isomers);
                    if aromatic {
                        match out.pop() {
                            None | Some(':') => {}
                            Some(c) => out.push(c),
                        }
                    }
                    'print_atom: {
                        if !force {
                            match atom.protons {
                                0 => match atom.isotope {
                                    0xFFFF => {
                                        out.push('*');
                                        break 'print_atom;
                                    }
                                    0xFFFE => {
                                        out.push(if aromatic { 'a' } else { 'A' });
                                        break 'print_atom;
                                    }
                                    0xFFFC => {
                                        out.push(if aromatic { 'q' } else { 'Q' });
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
                                    out.push(if aromatic { 'b' } else { 'B' });
                                    break 'print_atom;
                                }
                                6 => {
                                    out.push(if aromatic { 'c' } else { 'C' });
                                    break 'print_atom;
                                }
                                7 => {
                                    out.push(if aromatic { 'n' } else { 'N' });
                                    break 'print_atom;
                                }
                                8 => {
                                    out.push(if aromatic { 'o' } else { 'O' });
                                    break 'print_atom;
                                }
                                9 => {
                                    out.push('F');
                                    break 'print_atom;
                                }
                                15 => {
                                    out.push(if aromatic { 'p' } else { 'P' });
                                    break 'print_atom;
                                }
                                16 => {
                                    out.push(if aromatic { 's' } else { 'S' });
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
                        let idx = out.len();
                        out.push_str(ATOM_DATA[atom.protons as usize].sym);
                        if aromatic {
                            // ascii to ascii is safe
                            unsafe {
                                out.as_bytes_mut()[idx].make_ascii_lowercase();
                            }
                        }
                        if needs_h {
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
                    let idx = queue.len();
                    queue.extend(
                        graph
                            .edges(node)
                            .filter(|e| prev.map(|p| p.0) != Some(e.target()))
                            .map(|e| {
                                let w = *e.weight();
                                let weight = if aromatic && w == Bond::Aromatic {
                                    Bond::Single
                                } else {
                                    w
                                };
                                (Some(e.target()), Some((node, edge)), weight, true)
                            }),
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
                    visited.insert(
                        node,
                        (out.len(), prev.into_iter().map(|i| i.0).collect::<Vec<_>>()),
                    );
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
