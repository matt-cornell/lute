use crate::atom_info::ATOM_DATA;
use crate::core::*;
use petgraph::visit::*;
use petgraph::Undirected;
use std::ascii::Char;
use std::fmt::{self, Debug, Display, Formatter};

/// Convenience for sorting
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum AtomKind {
    Unknown,
    Carbon,
    Symbol(&'static str),
    Hydrogen,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct InChIKeyParseError(pub usize);
impl Display for InChIKeyParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let n = self.0;
        match n {
            0..14 | 15..23 | 26 => write!(f, "expected an ASCII char for character {n}"),
            14 | 25 => write!(f, "expected a hyphen for character {n}"),
            23 => f.write_str("expected 'S' for character 24"),
            24 => f.write_str("expected 'A' for character 25"),
            _ => write!(f, "expected a 27-character input, found one of length {n}"),
        }
    }
}

/// Fixed-size buffer for an InChIKey, assumes use of 1S.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InChIKey {
    pub first: [Char; 14],
    pub second: [Char; 8],
    pub last: Char,
}
impl InChIKey {
    /// Format this key as an array of ASCII characters that can easily be safely converted to a string or bytes.
    pub fn as_bytes(self) -> [Char; 27] {
        let mut out = [Char::HyphenMinus; 27];
        out[0..14].copy_from_slice(&self.first);
        out[15..23].copy_from_slice(&self.second);
        out[23] = Char::CapitalS;
        out[24] = Char::CapitalA;
        out[26] = self.last;
        out
    }
}
/// The `Debug` is overridden to be easily readable.
impl Debug for InChIKey {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_tuple("InChIKey")
            .field(&self.as_bytes().as_str())
            .finish()
    }
}
impl Display for InChIKey {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_bytes().as_str())
    }
}
impl FromStr for InChIKey {
    type Err = InChIKeyParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let first = if let Some(s) = s.get(0..14) {
            if let Some(c) = s.as_ascii() {
                c
            } else {
                for (n, b) in s.bytes().enumerate() {
                    if !b.is_ascii() {
                        return Err(InChIKeyParseError(n));
                    }
                }
                unreachable!()
            }
        } else {
            return Err(InChIKeyParseError(
                s.floor_char_boundary(std::cmp::min(13, s.len())),
            ));
        };
        if s.as_bytes().get(14) != Some(&b'-') {
            return Err(InChIKeyParseError(14));
        }
        let second = if let Some(s) = s.get(15..23) {
            if let Some(c) = s.as_ascii() {
                c
            } else {
                for (n, b) in s.bytes().enumerate() {
                    if !b.is_ascii() {
                        return Err(InChIKeyParseError(n + 15));
                    }
                }
                unreachable!()
            }
        } else {
            return Err(InChIKeyParseError(
                s.floor_char_boundary(std::cmp::min(22, s.len())),
            ));
        };
        if s.as_bytes().get(23) != Some(&b'S') {
            return Err(InChIKeyParseError(23));
        }
        if s.as_bytes().get(24) != Some(&b'A') {
            return Err(InChIKeyParseError(24));
        }
        if s.as_bytes().get(25) != Some(&b'-') {
            return Err(InChIKeyParseError(25));
        }
        let Some(last) = s.as_bytes().get(26).and_then(u8::as_ascii) else {
            return Err(InChIKeyParseError(26));
        };
        if s.len() > 27 {
            return Err(InChIKeyParseError(27));
        }
        Ok(Self {
            first,
            second,
            last,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct InChI {
    pub formula: Vec<(u8, usize)>,
    pub conns: Vec<u8>,
}
impl InChI {
    pub fn new<G>(graph: G) -> Self
    where
        G: NodeCompactIndexable
            + Data<NodeWeight = Atom, EdgeWeight = Bond>
            + GraphProp<EdgeType = Undirected>
            + IntoNodeReferences,
    {
        let mut atoms: Vec<(Atom, G::NodeId, (u8, SmallVec<u8, 1>))> =
            Vec::with_capacity(graph.node_bound());
        atoms.extend(graph.node_references().map(|node| {
            let atom = *node.weight();
            (
                atom,
                node.id(),
                (0u8, smallvec![atom.data.single() + atom.data.other()]),
            )
        }));
        atoms.sort_unstable_by_key(|a| match a.0.protons {
            0 => AtomKind::Unknown,
            1 => AtomKind::Hydrogen,
            6 => AtomKind::Carbon,
            n => AtomKind::Symbol(ATOM_DATA[n as usize].sym),
        });
        let mut color = 0u8;
        let mut old = Vec::new();
        let mut last = None;
        for i in &mut atoms {
            let w = i.2 .1[0];
            if last != Some(w) {
                color += 1;
                last = Some(w);
            }
            i.2 .0 = color;
        }
        let mut ns = Vec::new();
        loop {
            if atoms == old {
                let l = atoms.len();
                let mut done = true;
                for i in 0..(l - 1) {
                    let next = atoms[i + 1].2 .0;
                    let curr = &mut atoms[i].2 .0;
                    if *curr == next {
                        if i == 0 {
                            *curr = 0;
                        } else {
                            *curr = atoms[i - 1].2 .0 + 1;
                        }
                        done = false;
                    }
                }
                if done {
                    break;
                }
            }
            old.clone_from(&atoms);
            color = 0;
            last = None;
            for n in &mut atoms {
                ns.clear();
                ns.extend(graph.neighbors(n.1));
                let weight = old
                    .iter()
                    .filter_map(|n| ns.contains(&n.1).then_some(n.2 .0))
                    .sum();
                n.2 .0 = weight;
            }
            atoms.sort_unstable_by_key(|i| &i.2);
            for i in &mut atoms {
                let w = i.2 .1[0];
                if last != Some(w) {
                    color += 1;
                    last = Some(w);
                }
                i.2 .0 = color;
            }
        }
    }
}
impl Display for InChI {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("InChI=1S/")?;
        for &(p, c) in &self.formula {
            write!(f, "{}{c}", ATOM_DATA[p as usize].sym)?;
        }
        {
            let mut first = true;
            for i in &self.conns {
                if first {
                    first = false;
                } else {
                    f.write_str("-")?;
                }
                write!(f, "{i}")?;
            }
        }
        Ok(())
    }
}
