#![allow(dead_code, unused_variables)]

use crate::core::*;
use crate::graph::misc::DataValueMap;
use crate::graph::CycleBasis;
use crate::utils::bitset::BitSet;
use arrayvec::ArrayString;
use itertools::Itertools;
use petgraph::visit::*;
use smallvec::{smallvec_inline, SmallVec};
use std::cell::Cell;
use std::fmt::Write;
use thiserror::Error;

const MAX_NUM_LEN: usize = 28;

bitflags::bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct IupacConfig: u8 {
        /// No substitutions, just the PIN.
        const PREFERRED = 0;
        /// Use iso-, neo-, and tert- prefixes
        const ISO_NEO_TERT = 1;
        /// Don't generate names like "dimethyl ether", prefer "methoxymethane".
        const NEVER_ETHER = 2;
        /// Give chemicals more common names if possible.
        const COMMON_SUBS = 4;
        /// Give names like "formic acid", forces -aldehyde suffix for aldehydes.
        const FORM_ACET = 8;
        /// Don't have distinct R/S and E/Z stereoisomers.
        const NO_STEREO = 16;
        /// No numbers, e.g. 1-propanol and 2-propanol both become "propanol".
        /// "Isopropanol" could still be generated with this set.
        const NO_NUMBERS = 32;
        /// Generate "hyrochloric acid" for HCl instead of "hydrogen chloride"
        const HYDRO_ACIDS = 64;
        /// Combination of options to give the most common names.
        const COMMON_NAME = 29;
    }
}

/// Kind of a count
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum NumKind {
    /// Count of a group (*di*chloride)
    Count {
        /// If false, mono gives an empty output
        necessary: bool,
    },
    /// Length of a chain (*eth*ane)
    Length,
    /// Alternate length (*form*aldehyde)
    AltLength,
}

/// A halogen
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
enum Halogen {
    Hal,
    At,
    Br,
    Cl,
    F,
    I,
}
impl Halogen {
    pub fn from_u8(p: u8) -> Self {
        match p {
            0 => Self::Hal,
            9 => Self::F,
            17 => Self::Cl,
            35 => Self::Br,
            53 => Self::I,
            85 => Self::At,
            _ => unreachable!(),
        }
    }
    pub fn fragment(self) -> &'static str {
        match self {
            Self::Hal => "hal",
            Self::At => "astat",
            Self::Br => "brom",
            Self::Cl => "chlor",
            Self::F => "fluor",
            Self::I => "iod",
        }
    }
}

/// Substitution kind of an oxygen (for oxyhalides)
#[derive(Debug, Clone, Copy)]
enum SubKind {
    /// Substituted onto another part of a molecule (methyl perchlorate)
    Mol,
    /// Conjugate base (perchlorate anion)
    Base,
    /// Radical (perchlorate radical)
    R,
    /// Acid (perchloric acid)
    H,
}

/// Something's wrong with the molecule passed into `iupac_name`
#[derive(Debug, Clone, Error)]
pub enum InvalidMolecule<N> {
    #[error("Invalid valence on atom: expected one of {expected:?}, found {num}")]
    InvalidValence {
        atom: N,
        num: isize,
        expected: &'static [u8],
    },
    #[error("Invalid charge on atom: expected one of {expected:?}, found {charge:+}")]
    InvalidCharge {
        atom: N,
        charge: i8,
        expected: &'static [i8],
    },
    #[error("Too many of {radical} ({num}), IUPAC doesn't specify how to name more than 10,000")]
    TooManyOfSub { radical: String, num: usize },
    #[error("Unimplemented feature: {0}")]
    Unimplemented(&'static str),
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct RingScore {
    atoms: usize,
    hetero_n: usize,
    hetero_o: usize,
    hetero_s: usize,
    hetero_other: usize,
    bond_score: usize,
}

/// Generate the name for a number. Returns true if NumKind::AltLength was given, and an alternate name was given.
fn num_name(n: usize, kind: NumKind, vowel_end: bool) -> (ArrayString<MAX_NUM_LEN>, bool) {
    let mut out = ArrayString::new();
    assert!(n < 10000, "only chains less than 10000 long are supported!");
    if n < 10 {
        match kind {
            NumKind::Count { necessary } => {
                let names = [
                    "",
                    if !necessary { "" } else { "mon" },
                    "di",
                    "tri",
                    "tetr",
                    "pent",
                    "hex",
                    "hept",
                    "oct",
                    "non",
                ];
                out.push_str(names[n]);
                out.push_str(names[n]);
                if vowel_end && n > 3 {
                    out.push('a');
                }
            }
            NumKind::Length => {
                let names = [
                    "", "meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non",
                ];
                out.push_str(names[n]);
                if vowel_end && n > 3 {
                    out.push('a');
                }
            }
            NumKind::AltLength => {
                let names = [
                    "", "form", "acet", "propion", "butyr", "valer", "capron", "enanthal", "oct",
                    "non",
                ];
                out.push_str(names[n]);
                return (out, n < 8);
            }
        };
        return (out, false);
    }
    let last_two = n % 100;
    const DIGITS: [&str; 10] = [
        "", "hen", "do", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona",
    ];
    match last_two / 10 {
        0 | 1 => {
            out.push_str(DIGITS[last_two]);
            if last_two >= 10 {
                out.push_str("dec");
            }
        }
        2 => {
            let last = last_two % 10;
            if last == 0 {
                out.push_str("icos");
            } else {
                out.push_str(DIGITS[last]);
                out.push_str("cos");
            }
        }
        3 => {
            let last = last_two % 10;
            out.push_str(DIGITS[last]);
            out.push_str("triacont");
        }
        n => {
            let last = last_two % 10;
            out.push_str(DIGITS[last]);
            out.push_str(DIGITS[n]);
            out.push_str("cont");
        }
    }
    out.push_str(
        [
            "", "hect", "dict", "trict", "tetract", "pentact", "hexact", "heptact", "octact",
            "nonact",
        ][n / 100 % 10],
    );
    out.push_str(
        [
            "", "kili", "dili", "trili", "tetrali", "pentali", "hexali", "heptali", "octali",
            "nonali",
        ][n / 1000],
    );
    (out, false)
}

fn name_radical<G>(
    graph: G,
    bits: &mut BitSet<usize, 2>,
    out: &mut String,
    from: Option<[G::NodeId; 2]>,
) -> Result<usize, InvalidMolecule<G::NodeId>>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond> + IntoEdges,
{
    Err(InvalidMolecule::Unimplemented("radical naming"))
}

fn score_ring<G>(graph: G, nodes: &[(usize, G::NodeId)]) -> RingScore
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond> + IntoEdges,
{
    match nodes {
        [] => Default::default(),
        [(_, n)] => {
            let a = graph.node_weight(*n).unwrap();
            RingScore {
                atoms: 1,
                hetero_n: (a.protons == 7) as _,
                hetero_o: (a.protons == 8) as _,
                hetero_s: (a.protons == 16) as _,
                hetero_other: ![6, 7, 8, 16].contains(&a.protons) as _,
                bond_score: 0,
            }
        }
        &[(_, n1), (_, n2)] => {
            let a1 = graph.node_weight(n1).unwrap();
            let a2 = graph.node_weight(n2).unwrap();
            RingScore {
                atoms: 2,
                hetero_n: (a1.protons == 7) as usize + (a2.protons == 7) as usize,
                hetero_o: (a1.protons == 8) as usize + (a2.protons == 8) as usize,
                hetero_s: (a1.protons == 16) as usize + (a2.protons == 16) as usize,
                hetero_other: ![6, 7, 8, 16].contains(&a1.protons) as usize
                    + ![6, 7, 8, 16].contains(&a2.protons) as usize,
                bond_score: graph
                    .edges(n1)
                    .find(|e| e.target() == n2)
                    .unwrap()
                    .weight()
                    .bond_count()
                    .ceil() as _,
            }
        }
        _ => {
            let mut bond_score = 0.0f32;
            let mut hetero_n = 0;
            let mut hetero_o = 0;
            let mut hetero_s = 0;
            let mut hetero_other = 0;
            for (&(_, n1), &(_, n2)) in nodes.iter().circular_tuple_windows() {
                bond_score += graph
                    .edges(n1)
                    .find(|e| e.target() == n2)
                    .unwrap()
                    .weight()
                    .bond_count();
                match graph.node_weight(n1).unwrap().protons {
                    6 => {}
                    7 => hetero_n += 1,
                    8 => hetero_o += 1,
                    16 => hetero_s += 1,
                    _ => hetero_other += 1,
                }
            }
            RingScore {
                atoms: nodes.len(),
                bond_score: bond_score.ceil() as _,
                hetero_n,
                hetero_o,
                hetero_s,
                hetero_other,
            }
        }
    }
}

fn name_ring_system<G>(
    graph: G,
    sys: &[Vec<(usize, G::NodeId)>],
    out: &mut String,
) -> Result<(), InvalidMolecule<G::NodeId>>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond> + IntoEdges,
{
    Err(InvalidMolecule::Unimplemented("ring naming"))
}

fn number_ring_system<G>(
    graph: G,
    sys: &mut [Vec<(usize, G::NodeId)>],
) -> Result<[bool; 2], InvalidMolecule<G::NodeId>>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond> + IntoEdges,
{
    Err(InvalidMolecule::Unimplemented("Ring numbering"))
}

/// Try to generate the IUPAC name for a molecule
pub fn iupac_name<G>(graph: G, cfg: IupacConfig) -> Result<String, InvalidMolecule<G::NodeId>>
where
    G: DataValueMap<NodeWeight = Atom, EdgeWeight = Bond>
        + IntoNodeReferences
        + NodeCompactIndexable
        + IntoEdges,
{
    let mut cycles = CycleBasis::new_struct(graph)
        .map(|i| {
            let ring =
                i.1.into_iter()
                    .map(|n| (0, graph.from_index(n)))
                    .collect::<Vec<_>>();
            (score_ring(graph, &ring), smallvec_inline![ring])
        })
        .collect::<Vec<_>>();
    {
        let mut i = 0;
        let mut atoms = BitSet::<usize, 1>::new();
        let mut bonds = BitSet::<usize, 1>::new();
        while i < cycles.len() {
            let mut common_atoms = 0;
            let mut common_het_n = 0;
            let mut common_het_o = 0;
            let mut common_het_s = 0;
            let mut common_het_other = 0;
            let mut common_bonds = 0.0;
            atoms.clear();
            bonds.clear();
            let mut j = 0;
            while j < i {
                let cy1 = &cycles[i].1[0];
                for cy2 in &cycles[j].1 {
                    for (n, (&(_, n1), &(_, n2))) in cy1.iter().circular_tuple_windows().enumerate()
                    {
                        if let Some(idx) = cy2.iter().position(|&n| n.1 == n1) {
                            common_atoms += 1;
                            if !atoms.set(n, true) {
                                match graph.node_weight(n1).unwrap().protons {
                                    6 => {}
                                    7 => common_het_n += 1,
                                    8 => common_het_o += 1,
                                    16 => common_het_s += 1,
                                    _ => common_het_other += 1,
                                }
                            }
                            let l = cy2.len();
                            if cy2[(idx + 1) % l].1 == n2
                                || cy2[(idx + l - 1) % l].1 == n2 && !bonds.set(n, true)
                            {
                                common_bonds += graph
                                    .edges(n1)
                                    .find(|n| n.target() == n2)
                                    .unwrap()
                                    .weight()
                                    .bond_count();
                            }
                        }
                    }
                }
                if common_atoms > 0 {
                    let (
                        RingScore {
                            atoms,
                            hetero_n,
                            hetero_o,
                            hetero_s,
                            hetero_other,
                            bond_score,
                        },
                        ring,
                    ) = cycles.remove(j);
                    i -= 1;
                    let (score, rings) = &mut cycles[i];
                    score.atoms += atoms;
                    score.hetero_n += hetero_n;
                    score.hetero_o += hetero_o;
                    score.hetero_s += hetero_s;
                    score.hetero_other += hetero_other;
                    score.bond_score += bond_score;
                    rings.extend(ring);
                } else {
                    j += 1;
                }
            }
            let score = &mut cycles[i].0;
            score.atoms -= common_atoms;
            score.hetero_n -= common_het_n;
            score.hetero_o -= common_het_o;
            score.hetero_s -= common_het_s;
            score.hetero_other -= common_het_other;
            score.bond_score -= common_bonds.ceil() as usize;
            i += 1;
        }
    }
    cycles.sort_by_key(|c| c.0);
    let mut frags = SmallVec::<String, 1>::new();
    let mut subs_buf = SmallVec::<_, 2>::new();
    let mut string_reuse = Vec::new();
    let mut unvisited = graph
        .node_identifiers()
        .map(|i| graph.to_index(i))
        .collect::<BitSet<usize, 2>>();

    while !unvisited.all_zero() {
        if let Some((_, mut base)) = cycles.pop() {
            let [can_rot, can_flip] = number_ring_system(graph, &mut base)?;
            for ring in &base {
                for &(num, node) in ring {
                    for ed in graph.edges(node) {
                        let new = ed.target();
                        if ring.iter().any(|n| n.1 == new)
                            || base.iter().any(|r| r.iter().any(|n| n.1 == new))
                        {
                            continue;
                        }
                        let mut s: String = string_reuse.pop().unwrap_or_default();
                        s.clear();
                        let priority =
                            name_radical(graph, &mut unvisited, &mut s, Some([new, node]))?;
                        subs_buf.push((priority, smallvec_inline![num], s, Cell::new(None)));
                    }
                }
            }
            if can_rot || can_flip {
                // TODO: rotate and flip rings for optimal numbering
            }
            subs_buf.sort_unstable_by(|a, b| a.1.cmp(&b.1));
            subs_buf.dedup_by(|a, b| {
                a.2 == b.2 && {
                    string_reuse.push(std::mem::take(&mut a.2));
                    b.1.extend_from_slice(&a.1);
                    true
                }
            });
            subs_buf.sort_unstable_by(|a, b| {
                let an = a.3.get().unwrap_or_else(|| {
                    let name = num_name(a.1.len(), NumKind::Count { necessary: false }, false).0;
                    a.3.set(Some(name));
                    name
                });
                let bn = b.3.get().unwrap_or_else(|| {
                    let name = num_name(b.1.len(), NumKind::Count { necessary: false }, false).0;
                    b.3.set(Some(name));
                    name
                });
                an.cmp(&bn).then_with(|| a.2.cmp(&b.2))
            });
            let mut out = string_reuse.pop().unwrap_or_default();
            out.clear();
            for (_, mut subs, name, prefix) in subs_buf.drain(..) {
                if !cfg.contains(IupacConfig::NO_NUMBERS) {
                    if out.is_empty() {
                        out.push('-');
                    }
                    subs.sort_unstable();
                    let mut first = true;
                    for i in subs {
                        if first {
                            first = false;
                        } else {
                            out.push(',');
                        }
                        let _ = write!(out, "{i}");
                    }
                }
                out.push_str(&prefix.get().unwrap());
                out.push_str(&name);
                string_reuse.push(name);
            }
            name_ring_system(graph, &base, &mut out)?;
            frags.push(out);
        } else {
            let mut out = string_reuse.pop().unwrap_or_default();
            out.clear();
            name_radical(graph, &mut unvisited, &mut out, None)?;
        }
    }
    Ok("".to_string())
}
