use petgraph::prelude::*;
use std::str::FromStr;
use std::collections::HashMap;
use thiserror::Error;
use std::borrow::Cow;
use SmilesErrorKind::*;
use strum::*;
use atoi::FromRadix10;
use std::fmt::{self, Display, Formatter};
use crate::atom_info::ATOM_DATA;

/// An atom in the molecule graph
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Atom {
    pub protons: u8,
    pub isotope: Option<u16>,
    pub charge: i8,
    /// Aromaticity can be found in the bonds, and that should checked. This is only here to
    /// simplify SMILES parsing
    aromatic: bool,
}
impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(f, "{}", ATOM_DATA[self.protons as usize].name)?;
            if let Some(isotope) = self.isotope {
                write!(f, "-{isotope}")?;
            }
        }
        else {
            use fmtastic::Superscript;
            if let Some(isotope) = self.isotope {
                write!(f, "{}", Superscript(isotope))?;
            }
            write!(f, "{}", ATOM_DATA[self.protons as usize].sym)?;
            match self.charge {
                0 => {},
                1 => f.write_str("⁺")?,
                -1 => f.write_str("⁻")?,
                _ => write!(f, "{:+}", Superscript(self.charge))?
            }
        }
        Ok(())
    }
}

/// A bond between atoms in the molecule graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, AsRefStr, Display)]
pub enum Bond {
    /// Non-bond, shouldn't appear in final graph
    #[strum(serialize = "non")]
    Non,
    #[strum(serialize = "single")]
    Single,
    #[strum(serialize = "double")]
    Double,
    #[strum(serialize = "triple")]
    Triple,
    #[strum(serialize = "quadruple")]
    Quad,
    #[strum(serialize = "aromatic")]
    Aromatic,
    #[strum(serialize = "left")]
    Left,
    #[strum(serialize = "right")]
    Right,
}
impl Bond {
    pub fn bond_count(self) -> f32 {
        match self {
            Self::Single | Self::Left | Self::Right => 1f32,
            Self::Double => 2f32,
            Self::Triple => 3f32,
            Self::Quad => 4f32,
            Self::Aromatic => 1.5f32,
            Self::Non => 0f32
        }
    }
}

/// A molecule graph is an undirected graph between atoms, connected with bonds
pub type MoleculeGraph = UnGraph<Atom, Bond>;

/// Newtype around `petgraph`'s graph so some traits can be implemented

#[derive(Debug, Clone)]
pub struct Molecule(pub MoleculeGraph);
impl Molecule {
    /// Parse the string as SMILES, but with the addition that `R` acts as a wildcard
    pub fn from_smiles(s: &str) -> Result<Self, SmilesError> {
        SmilesParser::new(s).parse().map(Self)
    }
}
/// It's better to use the `from_smiles` function, because it allows for borrowing from the source
/// string in the error. Using `from_str` will force an allocation.
impl FromStr for Molecule {
    type Err = SmilesError<'static>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_smiles(s).map_err(SmilesError::into_owned)
    }
}

/// Inner enum for `SmilesError`
#[derive(Debug, Clone, Error)]
pub enum SmilesErrorKind<'a> {
    #[error("{0} is not a recognized element")]
    UnknownElement(Cow<'a, bstr::BStr>),
    #[error("a loop was opened without an atom")]
    LoopWithoutAtom,
    #[error("loop {0} was not closed by the end of the formula")]
    UnclosedLoop(usize),
    #[error("expected a ring number")]
    ExpectedRingNumber,
    #[error("expected an atom")]
    ExpectedAtom,
    #[error("expected a closing bracket")]
    ExpectedClosingBracket,
    #[error("bonds in loop don't match: {0} vs {1}")]
    LoopBondMismatch(Bond, Bond),
    #[error("duplicate bonds between atoms")]
    DuplicateBond,
}
impl SmilesErrorKind<'_> {
    pub fn into_owned(self) -> SmilesErrorKind<'static> {
        use SmilesErrorKind::*;
        match self {
            UnknownElement(e) => UnknownElement(e.into_owned().into()),
            LoopWithoutAtom => LoopWithoutAtom,
            UnclosedLoop(i) => UnclosedLoop(i),
            ExpectedRingNumber => ExpectedRingNumber,
            ExpectedAtom => ExpectedAtom,
            ExpectedClosingBracket => ExpectedClosingBracket,
            LoopBondMismatch(b1, b2) => LoopBondMismatch(b1, b2),
            DuplicateBond => DuplicateBond
        }
    }
}

/// Something went wrong trying to parse a SMILES string
#[derive(Debug, Clone, Error)]
#[error("an error occured at byte {index} in the SMILES string: {kind}")]
pub struct SmilesError<'a> {
    pub index: usize,
    pub kind: SmilesErrorKind<'a>
}
impl<'a> SmilesError<'a> {
    /// Convenience method
    pub const fn new(index: usize, kind: SmilesErrorKind<'a>) -> Self {
        Self {
            index,
            kind,
        }
    }
    /// Convert self into a static lifetime by heap-allocating the string
    pub fn into_owned(self) -> SmilesError<'static> {
        SmilesError {
            index: self.index,
            kind: self.kind.into_owned()
        }
    }
}

/// Parser for a SMILES string
pub struct SmilesParser<'a> {
    // use byte slice because SMILES shouldn't have non-ASCII data
    input: &'a [u8],
    index: usize,
    rings: HashMap<usize, (NodeIndex, Option<Bond>)>,
    graph: MoleculeGraph
}
impl<'a> SmilesParser<'a> {
    pub fn new<I: AsRef<[u8]> + ?Sized>(input: &'a I) -> Self {
        let input = input.as_ref();
        debug_assert!(input.is_ascii());
        Self {
            input,
            index: 0,
            rings: Default::default(),
            graph: Default::default()
        }
    }
    /// Parse a "chain". This can really be anything, though, it just returns the first atom in the
    /// group so it can be bonded to something else
    fn parse_chain(&mut self, nested: bool) -> Result<Option<(NodeIndex, Bond, bool)>, SmilesError<'a>> {
        if self.index >= self.input.len() {
            return Ok(None);
        }
        let first_bond = if nested {
            self.get_bond()
        } else {
            (Bond::Non, false)
        };
        
        let first_atom = self.get_atom()?.unwrap();
        let mut last_atom = first_atom;
        while self.index < self.input.len() {
            match self.input[self.index] {
                b'(' => {
                    self.index += 1;
                    let start_idx = self.index;
                    let (atom, bond, ex) = self.parse_chain(true)?.ok_or(SmilesError::new(self.index, ExpectedAtom))?;
                    if self.graph.contains_edge(last_atom, atom) {
                        Err(SmilesError::new(start_idx, DuplicateBond))?
                    }
                    self.graph.add_edge(last_atom, atom, if !ex && self.graph[last_atom].aromatic && self.graph[atom].aromatic {Bond::Aromatic} else {bond});
                }
                b')' if nested => {
                    self.index += 1;
                    break
                }
                _ => {
                    let (bond, ex) = self.handle_loops(last_atom)?;
                    let atom = self.get_atom()?;
                    if let Some(atom) = atom {
                        if bond != Bond::Non {
                            self.graph.add_edge(last_atom, atom, if !ex && self.graph[last_atom].aromatic && self.graph[atom].aromatic {Bond::Aromatic} else {bond});
                        }
                        last_atom = atom;
                    }
                    else if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom))?
                    }
                    else {
                        break
                    }
                }
            }
        }
        Ok(Some((first_atom, first_bond.0, first_bond.1)))
    }
    /// Parse an atom or an atom with hydrogens attached (in brackets)
    fn get_atom(&mut self) -> Result<Option<NodeIndex>, SmilesError<'a>> {
        match self.input.get(self.index) {
            None => Ok(None),
            Some(&b'B') => {
                self.index += 1;
                if self.input.get(self.index) ==  Some(&b'r') {
                    self.index += 1;
                    Ok(Some(self.graph.add_node(Atom {protons: 35, isotope: None, charge: 0, aromatic: false})))
                } else {
                    Ok(Some(self.graph.add_node(Atom {protons: 5, isotope: None, charge: 0, aromatic: false})))
                }
            }
            Some(&b'C') => {
                self.index += 1;
                if self.input.get(self.index) ==  Some(&b'l') {
                    self.index += 1;
                    Ok(Some(self.graph.add_node(Atom {protons: 17, isotope: None, charge: 0, aromatic: false})))
                } else {
                    Ok(Some(self.graph.add_node(Atom {protons: 6, isotope: None, charge: 0, aromatic: false})))
                }
            }
            Some(&b'N') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 7, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'O') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 8, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'P') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 15, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'S') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 16, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'F') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 9, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'I') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 53, isotope: None, charge: 0, aromatic: false})))
            }
            Some(&b'R') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 0, isotope: Some(0), charge: 0, aromatic: false})))
            }
            Some(&b'n') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 7, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b'o') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 8, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b'p') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 15, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b's') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 16, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b'b') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 5, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b'c') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 6, isotope: None, charge: 0, aromatic: true})))
            }
            Some(&b'r') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom {protons: 0, isotope: Some(0), charge: 0, aromatic: true})))
            }
            Some(&b'[') => {
                self.index += 1;
                let (isotope, used) = u16::from_radix_10(&self.input[self.index..]);
                self.index += used;
                let mut isotope = (used != 0).then_some(isotope);
                let atom = match self.input.get(self.index) {
                    None => Err(SmilesError::new(self.index, ExpectedAtom))?,
                    Some(&b'b') => self.graph.add_node(Atom {protons: 5, isotope, charge: 0, aromatic: true}),
                    Some(&b'c') => self.graph.add_node(Atom {protons: 6, isotope, charge: 0, aromatic: true}),
                    Some(&b'n') => self.graph.add_node(Atom {protons: 7, isotope, charge: 0, aromatic: true}),
                    Some(&b'o') => self.graph.add_node(Atom {protons: 8, isotope, charge: 0, aromatic: true}),
                    Some(&b'p') => self.graph.add_node(Atom {protons: 15, isotope, charge: 0, aromatic: true}),
                    Some(&b's') => self.graph.add_node(Atom {protons: 16, isotope, charge: 0, aromatic: true}),
                    Some(&b'r') => self.graph.add_node(Atom {protons: 0, isotope: Some(isotope.unwrap_or(0)), charge: 0, aromatic: true}),
                    Some(c) if c.is_ascii_uppercase() => {
                        let start = self.index;
                        let len = self.input[(self.index + 1)..].iter().copied().take_while(u8::is_ascii_lowercase).count();
                        let elem = &self.input[start..(start + len + 1)];
                        self.index += len;
                        let protons = if elem == b"R" {
                            isotope.get_or_insert(0);
                            0
                        } else {
                            ATOM_DATA.iter().enumerate().find(|(_, a)| a.sym.as_bytes() == elem).ok_or(SmilesError::new(start, UnknownElement(bstr::BStr::new(elem).into())))?.0 as _
                        };
                        self.graph.add_node(Atom {protons, isotope, charge: 0, aromatic: false})
                    }
                    _ => Err(SmilesError::new(self.index, ExpectedAtom))?
                };
                self.index += 1;
                if self.input.get(self.index) == Some(&b'H') {
                    self.index += 1;
                    let (mut h, used) = usize::from_radix_10(&self.input[self.index..]);
                    if used == 0 {
                        h = 1;
                    }
                    else {
                        self.index += used;
                    }
                    for _ in 0..h {
                        let hy = self.graph.add_node(Atom {protons: 1, isotope: None, charge: 0, aromatic: false});
                        self.graph.add_edge(atom, hy, Bond::Single);
                    }
                }
                match self.input.get(self.index) {
                    Some(&b'+') => {
                        self.index += 1;
                        let (mut charge, used) = i8::from_radix_10(&self.input[self.index..]);
                        if used == 0 {
                            let count = self.input[self.index..].iter().take_while(|&&c| c == b'+').count();
                            self.index += count;
                            charge = (count + 1) as i8;
                        }
                        else {
                            self.index += used;
                        }
                        self.graph[atom].charge = charge;
                    }
                    Some(&b'-') => {
                        self.index += 1;
                        let (mut charge, used) = i8::from_radix_10(&self.input[self.index..]);
                        if used == 0 {
                            let count = self.input[self.index..].iter().take_while(|&&c| c == b'-').count();
                            self.index += count;
                            charge = (count + 1) as i8;
                        }
                        else {
                            self.index += used;
                        }
                        self.graph[atom].charge = -charge;
                    }
                    _ => {}
                }
                if self.input.get(self.index) == Some(&b']') {
                    self.index += 1;
                    Ok(Some(atom))
                } else {
                    Err(SmilesError::new(self.index, ExpectedClosingBracket))
                }
            }
            _ => Err(SmilesError::new(self.index, ExpectedAtom))
        }
    }
    /// Handle loops. Since this needs to know which bond to use, it also parses a bond.
    fn handle_loops(&mut self, last_atom: NodeIndex) -> Result<(Bond, bool), SmilesError<'a>> {
        let bond_idx = self.index;
        let prev_bond = self.get_bond();
        if self.index >= self.input.len() {
            return Ok(prev_bond);
        }
        let c = self.input[self.index];
        let num_idx = self.index;
        let num = if c.is_ascii_digit() {
            self.index += 1;
            (c - b'0') as usize
        } else if c == b'%' {
            self.index += 1;
            let (num, used) = usize::from_radix_10(&self.input[self.index..]);
            if used == 0 {
                return Err(SmilesError::new(self.index, ExpectedRingNumber));
            }
            self.index += used;
            num
        } else {
            return Ok(prev_bond);
        };
        use std::collections::hash_map::Entry;
        match self.rings.entry(num) {
            Entry::Occupied(e) => {
                let (other, old_bond) = e.remove();
                let (bond, ex) = match (prev_bond, old_bond) {
                    (bond, None) => bond,
                    ((_, false), Some(bond)) => (bond, true),
                    ((b1, true), Some(b2)) => {
                        if b1 == b2 {
                            (b1, true)
                        } else {
                            return Err(SmilesError::new(bond_idx, LoopBondMismatch(b1, b2)))
                        }
                    }
                };
                if self.graph.contains_edge(last_atom, other) {
                    return Err(SmilesError::new(num_idx, DuplicateBond));
                }
                if bond != Bond::Non {
                    self.graph.add_edge(last_atom, other, if !ex && self.graph[last_atom].aromatic && self.graph[other].aromatic {Bond::Aromatic} else {bond});
                }
            }
            Entry::Vacant(e) => {
                e.insert((last_atom, prev_bond.1.then_some(prev_bond.0)));
            }
        }
        Ok(self.get_bond())
    }
    /// Parse a single bond
    fn get_bond(&mut self) -> (Bond, bool) {
        let (bond, incr) = match self.input.get(self.index) {
            Some(&b'.') => (Bond::Non, true),
            Some(&b'-') => (Bond::Single, true),
            Some(&b'=') => (Bond::Double, true),
            Some(&b'#') => (Bond::Triple, true),
            Some(&b'$') => (Bond::Quad, true),
            Some(&b':') => (Bond::Aromatic, true),
            Some(&b'/') => (Bond::Left, true),
            Some(&b'\\') => (Bond::Right, true),
            _ => (Bond::Single, false)
        };
        if incr {
            self.index += 1;
        }
        (bond, incr)
    }
    /// Saturate all atoms with hydrogens
    fn saturate(&mut self) {
        for atom in self.graph.node_indices() {
            let ex_bonds = match self.graph[atom].protons {
                x @ 6..=9 => Some(10 - x),
                x @ 14..=17 => Some(18 - x),
                _ => None
            };
            if let Some(ex_bonds) = ex_bonds {
                let bond_count = self.graph.edges(atom).fold(0f32, |c, b| c + b.weight().bond_count()).floor().clamp(0.0, 255.0) as u8;
                if bond_count < ex_bonds {
                    for _ in 0..(ex_bonds - bond_count) {
                        let hy = self.graph.add_node(Atom {protons: 1, isotope: None, charge: 0, aromatic: false});
                        self.graph.add_edge(atom, hy, Bond::Single);
                    }
                }
            }
        }
    }
    /// Parse the molecule, consuming self. This is taken by value to avoid cleanup.
    pub fn parse(mut self) -> Result<MoleculeGraph, SmilesError<'a>> {
        self.parse_chain(false)?;
        self.saturate();
        Ok(self.graph)
    }
}
