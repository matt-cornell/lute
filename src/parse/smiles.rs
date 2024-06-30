use crate::atom_info::ATOM_DATA;
use crate::core::{TooManyBonds, *};
use crate::utils::echar::*;
use atoi::FromRadix10;
use itertools::Itertools;
use petgraph::prelude::*;
use petgraph::visit::*;
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::collections::HashMap;
use thiserror::Error;
use tracing::*;
use SmilesErrorKind::*;

#[macro_export]
macro_rules! smiles {
    ($smiles:literal) => {
        $crate::parse::smiles::SmilesParser::new($smiles)
            .parse()
            .expect(concat!("Failed to parse SMILES ", $smiles))
    };
}

/// Inner enum for `SmilesError`
#[derive(Debug, Clone, Copy, Error)]
pub enum SmilesErrorKind {
    #[error(transparent)]
    TooManyBonds(#[from] TooManyBonds),
    #[error("{0} is not a recognized element")]
    UnknownElement(EChar),
    #[error("a loop was opened without an atom")]
    LoopWithoutAtom,
    #[error("loop {0} was not closed by the end of the formula")]
    UnclosedLoop(usize),
    #[error("expected a ring number")]
    ExpectedRingNumber,
    #[error("expected an atom, found {}", .0.map_or_else(|| "EOF".to_string(), |e| TextByte(e).to_string()))]
    ExpectedAtom(Option<u8>),
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
    #[error("isotopes can't be specified on an R-group")]
    GlobIsotope,
}

/// Something went wrong trying to parse a SMILES string
#[derive(Debug, Clone, Copy, Error)]
#[error("an error occured at {} in the SMILES string: {kind}", IdxPrint(*.index))]
pub struct SmilesError {
    pub index: usize,
    pub kind: SmilesErrorKind,
}
impl SmilesError {
    /// Convenience method
    pub const fn new(index: usize, kind: SmilesErrorKind) -> Self {
        Self { index, kind }
    }
}
impl From<TooManyBonds> for SmilesError {
    fn from(value: TooManyBonds) -> Self {
        Self::new(usize::MAX, value.into())
    }
}

#[derive(Debug, Clone)]
/// Parser for a SMILES string
pub struct SmilesParser<'a> {
    // use byte slice because SMILES shouldn't have non-ASCII data
    pub input: &'a [u8],
    index: usize,
    rings: HashMap<usize, (NodeIndex, Option<Bond>)>,
    graph: MoleculeGraph,
    pub suppress: bool,
    pub validate: bool,
}
impl<'a> SmilesParser<'a> {
    /// Create a new parser from an input string
    pub fn new<I: AsRef<[u8]> + ?Sized>(input: &'a I) -> Self {
        let input = input.as_ref();
        debug_assert!(input.is_ascii());
        Self {
            input,
            index: 0,
            rings: Default::default(),
            graph: Default::default(),
            suppress: true,
            validate: cfg!(debug_assertions),
        }
    }

    /// Create a new parser from an input string, without suppressing hydrogens and r-groups
    pub fn new_unsuppressed<I: AsRef<[u8]> + ?Sized>(input: &'a I) -> Self {
        let input = input.as_ref();
        debug_assert!(input.is_ascii());
        Self {
            input,
            index: 0,
            rings: Default::default(),
            graph: Default::default(),
            suppress: false,
            validate: true,
        }
    }

    pub fn with_suppression(mut self, suppress: bool) -> Self {
        self.suppress = suppress;
        self
    }
    pub fn with_validation(mut self, validate: bool) -> Self {
        self.validate = validate;
        self
    }
    pub fn set_suppression(&mut self, suppress: bool) -> &mut Self {
        self.suppress = suppress;
        self
    }
    pub fn set_validation(&mut self, validate: bool) -> &mut Self {
        self.validate = validate;
        self
    }

    /// Parse a "chain". This can really be anything, though, it just returns the first atom in the
    /// group so it can be bonded to something else
    #[instrument(level = "debug", skip(self), fields(self.input, self.index))]
    fn parse_chain(
        &mut self,
        nested: bool,
    ) -> Result<Option<(NodeIndex, Bond, bool)>, SmilesError> {
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
                        Err(SmilesError::new(self.index, ExpectedAtom(Some(b'('))))?;
                    }
                    self.index += 1;
                    let start_idx = self.index;
                    let (atom, bond, ex) = self.parse_chain(true)?.ok_or(SmilesError::new(
                        self.index,
                        ExpectedAtom(Some(self.input[start_idx])),
                    ))?;
                    if bond == Bond::Non {
                        continue;
                    }
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
                        if !ex
                            && self.graph[last_atom].data.scratch()
                                & self.graph[atom].data.scratch()
                                & 2
                                != 0
                        {
                            Bond::Aromatic
                        } else {
                            bond
                        },
                    );
                    self.graph[atom].with_scratch(|s| *s |= 8);
                }
                Some(&b')') if nested => {
                    if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom(Some(b')'))))?;
                    }
                    self.index += 1;
                    return Ok(Some((first_atom, first_bond.0, first_bond.1)));
                }
                Some(&b'&') => {
                    self.graph[last_atom].add_rs(1)?;
                    self.index += 1;
                }
                Some(_) => {
                    let atom = self.get_atom()?;
                    if let Some(atom) = atom {
                        if bond != Bond::Non {
                            self.graph.add_edge(
                                last_atom,
                                atom,
                                if !ex
                                    && self.graph[last_atom].data.scratch()
                                        & self.graph[atom].data.scratch()
                                        & 2
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
                        Err(SmilesError::new(
                            self.index,
                            ExpectedAtom(self.input.get(self.index).copied()),
                        ))?
                    } else {
                        break;
                    }
                }
                None => {
                    if ex {
                        Err(SmilesError::new(self.index, ExpectedAtom(None)))?
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
    #[instrument(level = "trace", skip_all, fields(self.input, self.index))]
    fn get_atom(&mut self) -> Result<Option<NodeIndex>, SmilesError> {
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
            Some(&b'*') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_isotope(0, 0xFFFF))))
            }
            Some(&b'A') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_isotope(0, 0xFFFE))))
            }
            Some(&b'Q') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_isotope(0, 0xFFFC))))
            }
            Some(&b'X') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_isotope(0, 0x0100))))
            }
            Some(&b'M') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_isotope(0, 0x042B))))
            }
            Some(&b'n') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(7, 2))))
            }
            Some(&b'o') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(8, 2))))
            }
            Some(&b'p') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(15, 2))))
            }
            Some(&b's') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(16, 2))))
            }
            Some(&b'b') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(5, 2))))
            }
            Some(&b'c') => {
                self.index += 1;
                Ok(Some(self.graph.add_node(Atom::new_scratch(6, 2))))
            }
            Some(&b'a') => {
                self.index += 1;
                Ok(Some(
                    self.graph.add_node(Atom::new_isotope_scratch(0, 0xFFFE, 2)),
                ))
            }
            Some(&b'q') => {
                self.index += 1;
                Ok(Some(
                    self.graph.add_node(Atom::new_isotope_scratch(0, 0xFFFC, 2)),
                ))
            }
            Some(&b'x') => {
                self.index += 1;
                Ok(Some(
                    self.graph.add_node(Atom::new_isotope_scratch(0, 0x0100, 2)),
                ))
            }
            Some(&b'm') => {
                self.index += 1;
                Ok(Some(
                    self.graph.add_node(Atom::new_isotope_scratch(0, 0x042B, 2)),
                ))
            }
            Some(&b'[') => {
                let _span = trace_span!("parsing element block", index = self.index).entered();
                self.index += 1;
                let (isotope, used) = u16::from_radix_10(&self.input[self.index..]);
                self.index += used;
                let atom = match self.input.get(self.index) {
                    None => Err(SmilesError::new(self.index, ExpectedAtom(None)))?,
                    Some(&b'b') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(5, isotope, 2)),
                    Some(&b'c') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(6, isotope, 2)),
                    Some(&b'n') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(7, isotope, 2)),
                    Some(&b'o') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(8, isotope, 2)),
                    Some(&b'p') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(15, isotope, 2)),
                    Some(&b's') => self
                        .graph
                        .add_node(Atom::new_isotope_scratch(16, isotope, 2)),
                    Some(&b'a') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope_scratch(0, 0xFFFE, 2))
                    }
                    Some(&b'q') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope_scratch(0, 0xFFFC, 2))
                    }
                    Some(&b'x') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope_scratch(0, 0x0100, 2))
                    }
                    Some(&b'm') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope_scratch(0, 0x042B, 2))
                    }
                    Some(&b'*') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope_scratch(0, 0xFFFF, 2))
                    }
                    Some(&b'A') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope(0, 0xFFFE))
                    }
                    Some(&b'Q') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope(0, 0xFFFC))
                    }
                    Some(&b'X') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope(0, 0x0100))
                    }
                    Some(&b'M') => {
                        if used > 0 {
                            Err(SmilesError::new(self.index, GlobIsotope))?
                        }
                        self.graph.add_node(Atom::new_isotope(0, 0x042B))
                    }
                    Some(c) if c.is_ascii_alphabetic() => {
                        let _span = trace_span!("parsing arbitrary element", start = c).entered();
                        let scratch = if c.is_ascii_lowercase() { 2 } else { 0 };
                        let c = c.to_ascii_uppercase();
                        let start = self.index;
                        let len = self.input[(self.index + 1)..]
                            .iter()
                            .copied()
                            .take_while(u8::is_ascii_lowercase)
                            .take(4)
                            .count();
                        let elem = &self.input[(start + 1)..(start + len + 1)];
                        if event_enabled!(Level::TRACE, elem, "found an element") {
                            trace!(
                                elem = String::from_utf8_lossy(elem).to_string(),
                                "found an element"
                            );
                        }
                        self.index += len;
                        let protons = ATOM_DATA
                            .iter()
                            .position(|a| {
                                a.sym.as_bytes()[0] == c && a.sym.as_bytes()[1..] == *elem
                            })
                            .ok_or(SmilesError::new(
                                start,
                                UnknownElement(
                                    EChar::new(&self.input[start..], (len + 1) as _).unwrap(),
                                ),
                            ))? as _;
                        self.graph
                            .add_node(Atom::new_isotope_scratch(protons, isotope, scratch))
                    }
                    c => Err(SmilesError::new(self.index, ExpectedAtom(c.copied())))?,
                };
                self.index += 1;
                if self.input.get(self.index) == Some(&b'@') {
                    trace!("handling chirality");
                    self.index += 1;
                    if self.input.get(self.index) == Some(&b'@') {
                        self.index += 1;
                        self.graph[atom].data.set_chirality(Chirality::R);
                    } else {
                        self.graph[atom].data.set_chirality(Chirality::S);
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
                        self.index += used;
                    }
                    trace!(count = h, "adding explicit hydrogens");
                    if self.suppress {
                        self.graph[atom].add_hydrogens(h)?;
                    } else {
                        for _ in 0..h {
                            let hy = self.graph.add_node(Atom::new(1));
                            self.graph.add_edge(atom, hy, Bond::Single);
                        }
                    }
                    self.graph[atom].with_scratch(|s| *s |= 1);
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
                        trace!(charge, "adding positive charge");
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
                        trace!(charge, "adding negative charge");
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
            c => Err(SmilesError::new(self.index, ExpectedAtom(c.copied()))),
        }
    }

    /// Handle loops. Since this needs to know which bond to use, it also parses a bond.
    #[instrument(level = "debug", skip_all, fields(self.input, self.index))]
    fn handle_loops(&mut self, last_atom: NodeIndex) -> Result<(Bond, bool), SmilesError> {
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
                                && self.graph[last_atom].data.scratch()
                                    & self.graph[other].data.scratch()
                                    & 2
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
    #[instrument(level = "trace", skip_all, fields(self.input, self.index))]
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

    /// Update bond counts
    #[instrument(level = "debug", skip_all, fields(self.input, self.index))]
    fn update_bonds(&mut self) -> Result<(), SmilesError> {
        for n in self.graph.node_indices() {
            let mut sc = 0;
            let mut bc = 0;
            for b in self.graph.edges(n) {
                if b.weight() == &Bond::Single {
                    sc += 1;
                } else {
                    bc += 1;
                }
            }
            self.graph[n].set_single_bonds(
                sc.try_into()
                    .map_err(|_| TooManyBonds(TooMany::Single, sc))?,
            )?;
            self.graph[n].set_other_bonds(
                bc.try_into()
                    .map_err(|_| TooManyBonds(TooMany::Other, bc))?,
            )?;
        }
        Ok(())
    }

    /// Saturate all atoms with hydrogens
    #[instrument(level = "debug", skip_all, fields(self.input, self.index))]
    fn update_hydrogens(&mut self) -> Result<(), SmilesError> {
        for atom in self.graph.node_indices() {
            let ex_bonds = match self.graph[atom].protons {
                x @ 6..=9 => Some(10 - (x as i8) + self.graph[atom].charge),
                x @ 14..=17 => Some(18 - (x as i8) + self.graph[atom].charge),
                _ => None,
            };
            let a = self.graph[atom];
            trace!(atom = %a, ex_bonds, "checking atom");
            if let Some(ex_bonds) = ex_bonds {
                let bond_count = (self
                    .graph
                    .edges(atom)
                    .fold(0f32, |c, b| c + b.weight().bond_count())
                    .ceil()
                    .clamp(0.0, 127.0) as i8)
                    + (a.data.hydrogen() + a.data.unknown()) as i8;
                if self.suppress {
                    let mut walk = self.graph.neighbors(atom).detach();
                    let mut count = 0;
                    while let Some((e, n)) = walk.next(&self.graph) {
                        if self.graph[n].protons == 1 && self.graph[e] == Bond::Single {
                            self.graph[n].add_hydrogens(1)?;
                            self.graph.remove_node(n);
                            count += 1;
                        }
                    }
                    if count > 0 {
                        debug!(count, id = atom.index(), "pruned hydrogens");
                    }
                    if bond_count < ex_bonds && self.graph[atom].data.scratch() & 1 == 0 {
                        self.graph[atom].add_hydrogens((ex_bonds - bond_count) as u8)?;
                    }
                } else if bond_count < ex_bonds && self.graph[atom].data.scratch() & 1 == 0 {
                    for _ in 0..(ex_bonds - bond_count) {
                        let hy = self.graph.add_node(Atom::new(1));
                        self.graph.add_edge(atom, hy, Bond::Single);
                    }
                }
            }
        }
        Ok(())
    }

    /// Update stereochemistry
    #[instrument(level = "debug", skip_all, fields(self.input, self.index))]
    fn update_stereo(&mut self) {
        use crate::molecule::*;

        // E/Z assignment

        {
            let _span = debug_span!("E/Z assignment").entered();
            for id in (0..self.graph.edge_count()).map(petgraph::graph::EdgeIndex::new) {
                if self.graph[id] != Bond::Double {
                    continue;
                }
                trace!(id = id.index(), "checking double bond");
                let (a, b) = self.graph.edge_endpoints(id).unwrap();
                let (mut left, mut right, mut unbound) = (None, None, None);
                let mut es = SmallVec::<_, 4>::new();
                for e in self.graph.edges(a) {
                    let other = if e.source() == a {
                        e.target()
                    } else {
                        e.source()
                    };
                    let r = match *e.weight() {
                        Bond::Left => &mut left,
                        Bond::Right => &mut right,
                        Bond::Single => &mut unbound,
                        _ => continue,
                    };
                    es.push(e.id());
                    *r = Some(self.graph.cip_priority(other, a));
                }
                let aord = match (left, right, unbound) {
                    (Some(lhs), Some(rhs), _)
                    | (Some(lhs), None, Some(rhs))
                    | (None, Some(rhs), Some(lhs)) => Ord::cmp(&lhs, &rhs),
                    (Some(_), None, None) => Ordering::Greater,
                    (None, Some(_), None) => Ordering::Less,
                    _ => Ordering::Equal,
                };
                trace!(id = id.index(), ord = aord as i8, "checked first end");
                for e in es.drain(..) {
                    self.graph[e] = Bond::Single;
                }
                if aord == Ordering::Equal {
                    continue;
                }
                (left, right, unbound) = (None, None, None);
                for e in self.graph.edges(b) {
                    let other = if e.source() == b {
                        e.target()
                    } else {
                        e.source()
                    };
                    let r = match *e.weight() {
                        Bond::Left => &mut left,
                        Bond::Right => &mut right,
                        Bond::Single => &mut unbound,
                        _ => continue,
                    };
                    es.push(e.id());
                    *r = Some(self.graph.cip_priority(other, b));
                }
                let bord = match (left, right, unbound) {
                    (Some(lhs), Some(rhs), _)
                    | (Some(lhs), None, Some(rhs))
                    | (None, Some(rhs), Some(lhs)) => Ord::cmp(&lhs, &rhs),
                    (Some(_), None, None) => Ordering::Greater,
                    (None, Some(_), None) => Ordering::Less,
                    _ => Ordering::Equal,
                };
                trace!(id = id.index(), ord = bord as i8, "checked other end");
                for e in es.drain(..) {
                    self.graph[e] = Bond::Single;
                }
                if bord == Ordering::Equal {
                    continue;
                }
                self.graph[id] = if aord == bord {
                    Bond::DoubleZ
                } else {
                    Bond::DoubleE
                };
            }
        }

        // R/S assignment
        {
            let _span = debug_span!("R/S assignment").entered();
            for id in (0..self.graph.node_count()).map(petgraph::graph::NodeIndex::new) {
                let ch = self.graph[id].data.chirality();
                if !ch.is_chiral() {
                    continue;
                }
                trace!(id = id.index(), "checking atom");

                let mut ns: SmallVec<_, 4> = self
                    .graph
                    .neighbors(id)
                    .map(|n| self.graph.cip_priority(n, id))
                    .collect();
                if ns.len() > 4 {
                    warn!(bonds = ?ns, "attempting to determine chirality for atom with more than 4 bonds");
                }
                if ch == Chirality::R {
                    ns.reverse();
                }
                if ns.len() == 4 {
                    let idx = ns.iter().position_min().unwrap();
                    ns.remove(idx);
                }
                let ch = if let [a, b, c] = &ns[..] {
                    'blk: {
                        let _span = trace_span!("comparing neighbors").entered();
                        let ab = a.cmp(b);
                        trace!(cmp = ab as i8, "a-b cmp");
                        if ab == Ordering::Equal {
                            break 'blk Chirality::None;
                        }
                        let bc = b.cmp(c);
                        trace!(cmp = bc as i8, "b-c cmp");
                        match (ab, bc) {
                            // a b c
                            (Ordering::Less, Ordering::Less) => Chirality::S,
                            // c b a
                            (Ordering::Greater, Ordering::Greater) => Chirality::R,
                            // b is max
                            (Ordering::Less, Ordering::Greater) => match c.cmp(a) {
                                Ordering::Equal => Chirality::None,
                                // c a b
                                Ordering::Greater => Chirality::S,
                                // a c b
                                Ordering::Less => Chirality::R,
                            },
                            (Ordering::Greater, Ordering::Less) => match c.cmp(a) {
                                Ordering::Equal => Chirality::None,
                                // b a c
                                Ordering::Greater => Chirality::R,
                                // b c a
                                Ordering::Less => Chirality::S,
                            },
                            (Ordering::Equal, _) | (_, Ordering::Equal) => Chirality::None,
                        }
                    }
                } else {
                    Chirality::None
                };
                std::mem::drop(ns);
                self.graph[id].data.set_chirality(ch);
            }
        }
    }

    /// Perform some checks on the molecule. Panics on failure (which should be impossible).
    #[instrument(level = "trace", skip_all, fields(self.input, self.index))]
    fn validate(&self) {
        for node in self.graph.node_references() {
            assert_eq!(
                self.graph.edges(node.id()).count(),
                (node.weight().data.single() + node.weight().data.other()) as usize
            );
        }
        for &edge in self.graph.edge_weights() {
            assert_ne!(edge, Bond::Non);
            assert_ne!(edge, Bond::Left);
            assert_ne!(edge, Bond::Right);
        }
    }

    /// Parse the molecule, consuming self. This is taken by value to avoid cleanup.
    #[instrument(level = "debug", skip_all)]
    pub fn parse(mut self) -> Result<MoleculeGraph, SmilesError> {
        self.parse_chain(false)?;
        self.update_bonds()?;
        self.update_hydrogens()?;
        self.update_stereo();
        if self.validate {
            self.validate();
        }
        if let Some(id) = self.rings.into_keys().next() {
            Err(SmilesError::new(self.index, UnclosedLoop(id)))
        } else {
            Ok(self.graph)
        }
    }
}
