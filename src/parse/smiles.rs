use crate::atom_info::ATOM_DATA;
use crate::molecule::*;
use atoi::FromRadix10;
use petgraph::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::collections::HashSet;
use thiserror::Error;
use SmilesErrorKind::*;

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
    #[error("expected a closing parenthesis")]
    ExpectedClosingParen,
    #[error("bonds in loop don't match: {0} vs {1}")]
    LoopBondMismatch(Bond, Bond),
    #[error("duplicate bonds between atoms")]
    DuplicateBond,
    #[error("multiple bonds on an R-group aren't allowed")]
    MultiBondedR,
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
            ExpectedClosingParen => ExpectedClosingParen,
            LoopBondMismatch(b1, b2) => LoopBondMismatch(b1, b2),
            DuplicateBond => DuplicateBond,
            MultiBondedR => MultiBondedR,
        }
    }
}

/// Something went wrong trying to parse a SMILES string
#[derive(Debug, Clone, Error)]
#[error("an error occured at byte {index} in the SMILES string: {kind}")]
pub struct SmilesError<'a> {
    pub index: usize,
    pub kind: SmilesErrorKind<'a>,
}
impl<'a> SmilesError<'a> {
    /// Convenience method
    pub const fn new(index: usize, kind: SmilesErrorKind<'a>) -> Self {
        Self { index, kind }
    }
    /// Convert self into a static lifetime by heap-allocating the string
    pub fn into_owned(self) -> SmilesError<'static> {
        SmilesError {
            index: self.index,
            kind: self.kind.into_owned(),
        }
    }
}

/// Parser for a SMILES string
pub struct SmilesParser<'a> {
    // use byte slice because SMILES shouldn't have non-ASCII data
    pub input: &'a [u8],
    index: usize,
    rings: HashMap<usize, (NodeIndex, Option<Bond>)>,
    graph: MoleculeGraph,
    pub suppress_hydrogens: bool,
}
impl<'a> SmilesParser<'a> {
    pub fn new<I: AsRef<[u8]> + ?Sized>(input: &'a I) -> Self {
        let input = input.as_ref();
        debug_assert!(input.is_ascii());
        Self {
            input,
            index: 0,
            rings: Default::default(),
            graph: Default::default(),
            suppress_hydrogens: true,
        }
    }
    pub fn new_unsuppressed<I: AsRef<[u8]> + ?Sized>(input: &'a I) -> Self {
        let input = input.as_ref();
        debug_assert!(input.is_ascii());
        Self {
            input,
            index: 0,
            rings: Default::default(),
            graph: Default::default(),
            suppress_hydrogens: false,
        }
    }
    /// Parse a "chain". This can really be anything, though, it just returns the first atom in the
    /// group so it can be bonded to something else
    fn parse_chain(
        &mut self,
        nested: bool,
    ) -> Result<Option<(NodeIndex, Bond, bool)>, SmilesError<'a>> {
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
            let (bond, ex) = self.handle_loops(last_atom)?;
            match self.input.get(self.index) {
                Some(&b'(') => {
                    if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom))?;
                    }
                    self.index += 1;
                    let start_idx = self.index;
                    let (atom, bond, ex) = self
                        .parse_chain(true)?
                        .ok_or(SmilesError::new(self.index, ExpectedAtom))?;
                    if self.graph.contains_edge(last_atom, atom) {
                        Err(SmilesError::new(start_idx, DuplicateBond))?
                    }
                    if self.graph[last_atom].protons == 0
                        && self.graph.edges(last_atom).next().is_some()
                    {
                        Err(SmilesError::new(start_idx, MultiBondedR))?
                    }
                    self.graph.add_edge(
                        last_atom,
                        atom,
                        if !ex && self.graph[last_atom].scratch & self.graph[atom].scratch & 2 != 0
                        {
                            Bond::Aromatic
                        } else {
                            bond
                        },
                    );
                }
                Some(&b')') if nested => {
                    if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom))?;
                    }
                    self.index += 1;
                    return Ok(Some((first_atom, first_bond.0, first_bond.1)));
                }
                Some(_) => {
                    let start_idx = self.index;
                    let atom = self.get_atom()?;
                    if let Some(atom) = atom {
                        if self.graph[last_atom].protons == 0
                            && self.graph.edges(last_atom).next().is_some()
                        {
                            Err(SmilesError::new(start_idx, MultiBondedR))?
                        }
                        if bond != Bond::Non {
                            self.graph.add_edge(
                                last_atom,
                                atom,
                                if !ex
                                    && self.graph[last_atom].scratch & self.graph[atom].scratch & 2
                                        != 0
                                {
                                    Bond::Aromatic
                                } else {
                                    bond
                                },
                            );
                        }
                        last_atom = atom;
                    } else if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom))?
                    } else {
                        break;
                    }
                }
                None => {
                    if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom))?
                    } else {
                        break;
                    }
                }
            }
        }
        if nested {
            Err(SmilesError::new(self.index, ExpectedClosingParen))
        } else {
            Ok(Some((first_atom, first_bond.0, first_bond.1)))
        }
    }
    /// Parse an atom or an atom with hydrogens attached (in brackets)
    fn get_atom(&mut self) -> Result<Option<NodeIndex>, SmilesError<'a>> {
        match self.input.get(self.index) {
            None => Ok(None),
            Some(&b'B') => {
                self.index += 1;
                if self.input.get(self.index) == Some(&b'r') {
                    self.index += 1;
                    Ok(Some(self.graph.add_node(Atom::new(35))))
                } else {
                    Ok(Some(self.graph.add_node(Atom::new(5))))
                }
            }
            Some(&b'C') => {
                self.index += 1;
                if self.input.get(self.index) == Some(&b'l') {
                    self.index += 1;
                    Ok(Some(self.graph.add_node(Atom::new(17))))
                } else {
                    Ok(Some(self.graph.add_node(Atom::new(6))))
                }
            }
            Some(&b'N') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(7))))
            }
            Some(&b'O') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(8))))
            }
            Some(&b'P') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(15))))
            }
            Some(&b'S') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(16))))
            }
            Some(&b'F') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(9))))
            }
            Some(&b'I') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(53))))
            }
            Some(&b'R') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new(0))))
            }
            Some(&b'n') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(7, 1))))
            }
            Some(&b'o') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(8, 1))))
            }
            Some(&b'p') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(15, 1))))
            }
            Some(&b's') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(16, 1))))
            }
            Some(&b'b') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(5, 1))))
            }
            Some(&b'c') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(6, 1))))
            }
            Some(&b'r') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(0, 1))))
            }
            Some(&b'[') => {
                self.index += 1;
                let (isotope, used) = u16::from_radix_10(&self.input[self.index..]);
                self.index += used;
                let atom = match self.input.get(self.index) {
                    None => Err(SmilesError::new(self.index, ExpectedAtom))?,
                    Some(&b'b') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(5, isotope, 1)),
                    Some(&b'c') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(6, isotope, 1)),
                    Some(&b'n') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(7, isotope, 1)),
                    Some(&b'o') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(8, isotope, 1)),
                    Some(&b'p') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(15, isotope, 1)),
                    Some(&b's') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(16, isotope, 1)),
                    Some(&b'r') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(0, isotope, 1)),
                    Some(c) if c.is_ascii_uppercase() => {
                        let start = self.index;
                        let len = self.input[(self.index + 1)..]
                            .iter()
                            .copied()
                            .take_while(u8::is_ascii_lowercase)
                            .count();
                        let elem = &self.input[start..(start + len + 1)];
                        self.index += len;
                        let protons = if elem == b"R" {
                            0
                        } else {
                            ATOM_DATA
                                .iter()
                                .enumerate()
                                .find(|(_, a)| a.sym.as_bytes() == elem)
                                .ok_or(SmilesError::new(
                                    start,
                                    UnknownElement(bstr::BStr::new(elem).into()),
                                ))?
                                .0 as _
                        };
                        self.graph.add_node(Atom::new_isotope(protons, isotope))
                    }
                    _ => Err(SmilesError::new(self.index, ExpectedAtom))?,
                };
                self.index += 1;
                if self.input.get(self.index) == Some(&b'@') {
                    self.index += 1;
                    if self.input.get(self.index) == Some(&b'@') {
                        self.index += 1;
                        self.graph[atom].chirality = Chirality::Cw;
                    } else {
                        self.graph[atom].chirality = Chirality::Ccw;
                    }
                }
                if self.input.get(self.index) == Some(&b'H') {
                    if self.graph[atom].protons == 0 {
                        Err(SmilesError::new(self.index, MultiBondedR))?
                    }
                    self.index += 1;
                    let (mut h, used) = u8::from_radix_10(&self.input[self.index..]);
                    if used == 0 {
                        h = 1;
                    } else {
                        self.index += used as usize;
                    }
                    if self.suppress_hydrogens {
                        self.graph[atom].hydrogens = h;
                    } else {
                        for _ in 0..h {
                            let hy = self.graph.add_node(Atom::new(1));
                            self.graph.add_edge(atom, hy, Bond::Single);
                        }
                    }
                    self.graph[atom].scratch |= 4;
                }
                match self.input.get(self.index) {
                    Some(&b'+') => {
                        self.index += 1;
                        let (mut charge, used) = i8::from_radix_10(&self.input[self.index..]);
                        if used == 0 {
                            let count = self.input[self.index..]
                                .iter()
                                .take_while(|&&c| c == b'+')
                                .count();
                            self.index += count;
                            charge = (count + 1) as i8;
                        } else {
                            self.index += used;
                        }
                        self.graph[atom].charge = charge;
                    }
                    Some(&b'-') => {
                        self.index += 1;
                        let (mut charge, used) = i8::from_radix_10(&self.input[self.index..]);
                        if used == 0 {
                            let count = self.input[self.index..]
                                .iter()
                                .take_while(|&&c| c == b'-')
                                .count();
                            self.index += count;
                            charge = (count + 1) as i8;
                        } else {
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
            _ => Err(SmilesError::new(self.index, ExpectedAtom)),
        }
    }
    /// Handle loops. Since this needs to know which bond to use, it also parses a bond.
    fn handle_loops(&mut self, last_atom: NodeIndex) -> Result<(Bond, bool), SmilesError<'a>> {
        loop {
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
                                return Err(SmilesError::new(bond_idx, LoopBondMismatch(b1, b2)));
                            }
                        }
                    };
                    if self.graph.contains_edge(last_atom, other) {
                        return Err(SmilesError::new(num_idx, DuplicateBond));
                    }
                    if bond != Bond::Non {
                        self.graph.add_edge(
                            last_atom,
                            other,
                            if !ex
                                && self.graph[last_atom].scratch & self.graph[other].scratch & 2
                                    != 0
                            {
                                Bond::Aromatic
                            } else {
                                bond
                            },
                        );
                    }
                }
                Entry::Vacant(e) => {
                    e.insert((last_atom, prev_bond.1.then_some(prev_bond.0)));
                }
            }
        }
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
            _ => (Bond::Single, false),
        };
        if incr {
            self.index += 1;
        }
        (bond, incr)
    }

    /// Saturate all atoms with hydrogens
    fn update_hydrogens(&mut self) {
        for atom in self.graph.node_indices() {
            let ex_bonds = match self.graph[atom].protons {
                x @ 6..=9 => Some(10 - (x as i8) + self.graph[atom].charge),
                x @ 14..=17 => Some(18 - (x as i8) + self.graph[atom].charge),
                _ => None,
            };
            if let Some(ex_bonds) = ex_bonds {
                let bond_count = (self
                    .graph
                    .edges(atom)
                    .fold(0f32, |c, b| c + b.weight().bond_count())
                    .ceil()
                    .clamp(0.0, 127.0) as i8)
                    + (self.graph[atom].hydrogens as i8);
                if self.suppress_hydrogens {
                    let mut walk = self.graph.neighbors(atom).detach();
                    while let Some((e, n)) = walk.next(&self.graph) {
                        if self.graph[n].protons == 1 && self.graph[e] == Bond::Single {
                            self.graph[n].hydrogens += 1;
                            self.graph.remove_node(n);
                        }
                    }
                    if bond_count < ex_bonds && self.graph[atom].scratch & 4 == 0 {
                        self.graph[atom].hydrogens += (ex_bonds - bond_count) as u8;
                    }
                } else if bond_count < ex_bonds && self.graph[atom].scratch & 4 == 0 {
                    for _ in 0..(ex_bonds - bond_count) {
                        let hy = self.graph.add_node(Atom::new(1));
                        self.graph.add_edge(atom, hy, Bond::Single);
                    }
                }
            }
        }
    }

    /// Update R-groups to avoid any duplicates
    fn update_rs(&mut self) {
        let mut found = HashSet::new();
        for n in self.graph.node_weights_mut() {
            if n.protons != 0 {
                continue;
            }
            if found.contains(&n.isotope) {
                n.isotope = (0..).find(|i| !found.contains(i)).unwrap();
            }
            found.insert(n.isotope);
        }
    }

    /// Parse the molecule, consuming self. This is taken by value to avoid cleanup.
    pub fn parse(mut self) -> Result<MoleculeGraph, SmilesError<'a>> {
        self.parse_chain(false)?;
        self.update_hydrogens();
        self.update_rs();
        if let Some(id) = self.rings.into_keys().next() {
            Err(SmilesError::new(self.index, UnclosedLoop(id)))
        } else {
            Ok(self.graph)
        }
    }
}
