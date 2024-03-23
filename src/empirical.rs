use crate::atom_info::*;
use crate::core::*;
use crate::utils::echar::*;
use bstr::ByteSlice;
use fmtastic::*;
use std::collections::BTreeMap;
use std::fmt::{self, Debug, Display, Formatter};
use std::ops::{Add, AddAssign};
use thiserror::Error;

#[macro_export]
macro_rules! empirical {
    ($form:literal) => {
        EmpiricalFormula::parse_str($form).expect(concat!("Failed to parse formula ", $form))
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Error)]
pub enum EmpiricalErrorKind {
    #[error("unknown atom {0}")]
    UnknownElement(EChar),
    #[error("unexpected character after charge: {0}")]
    AfterCharge(EChar),
    #[error("unexpected character: {0}")]
    UnexpectedChar(EChar),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Error)]
#[error("error at {} in empirical formula: {kind}", IdxPrint(*.idx))]
pub struct EmpiricalError {
    kind: EmpiricalErrorKind,
    idx: usize,
}
impl EmpiricalError {
    pub const fn new(idx: usize, kind: EmpiricalErrorKind) -> Self {
        Self { idx, kind }
    }
}

#[derive(Default, Clone, PartialEq, Eq, Hash)]
pub struct EmpiricalFormula {
    lower: [usize; 18],
    spill: BTreeMap<u8, usize>,
    pub charge: i8,
}
impl EmpiricalFormula {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn parse_str(s: impl AsRef<[u8]>) -> Result<Self, EmpiricalError> {
        use atoi::FromRadix10;
        use EmpiricalErrorKind::*;
        let input = s.as_ref();
        let mut idx = 0;
        let mut this = Self::new();
        loop {
            match input.get(idx) {
                None => return Ok(this),
                Some(&b'+' | &b'-') => {
                    let (mut c, u) = i8::from_radix_10(&input[(idx + 1)..]);
                    if u == 0 {
                        c = 1;
                    }
                    if input[idx] == b'-' {
                        c = -c;
                    }
                    idx += u + 1;
                    this.charge = c;
                    return if idx == input.len() {
                        Ok(this)
                    } else {
                        let rem = &input[idx..];
                        let len = rem.char_indices().next().unwrap().1;
                        Err(EmpiricalError::new(
                            idx,
                            AfterCharge(EChar::new(rem, len as _).unwrap()),
                        ))
                    };
                }
                Some(&c) if c.is_ascii_uppercase() => {
                    let start = idx;
                    let len = input[(idx + 1)..]
                        .iter()
                        .copied()
                        .take_while(u8::is_ascii_lowercase)
                        .take(3)
                        .count();
                    let elem = &input[start..(start + len + 1)];
                    idx += len + 1;
                    let protons = ATOM_DATA
                        .iter()
                        .enumerate()
                        .find(|(_, a)| a.sym.as_bytes() == elem)
                        .ok_or(EmpiricalError::new(
                            start,
                            UnknownElement(EChar::new(elem, (len + 1) as _).unwrap()),
                        ))?
                        .0 as _;
                    let (mut n, u) = usize::from_radix_10(&input[idx..]);
                    if u == 0 {
                        n = 1;
                    }
                    idx += u;
                    this.add_atom(protons, n);
                }
                Some(_) => {
                    let len = input[idx..].char_indices().next().unwrap().1;
                    Err(EmpiricalError::new(
                        idx,
                        UnexpectedChar(EChar::new(&input[idx..], len as _).unwrap()),
                    ))?
                }
            }
        }
    }

    pub fn get_atom(&mut self, atom: u8) -> usize {
        if atom < 18 {
            self.lower[atom as usize]
        } else {
            self.spill.get(&atom).copied().unwrap_or(0)
        }
    }
    pub fn set_atom(&mut self, atom: u8, count: usize) {
        if atom < 18 {
            self.lower[atom as usize] = count;
        } else if count == 0 {
            self.spill.remove(&atom);
        } else {
            self.spill.insert(atom, count);
        }
    }
    pub fn add_atom(&mut self, atom: u8, count: usize) {
        if atom < 18 {
            self.lower[atom as usize] += count;
        } else if count > 0 {
            self.spill
                .entry(atom)
                .and_modify(|c| *c += count)
                .or_insert(count);
        }
    }
}

impl Debug for EmpiricalFormula {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        struct AtomsHelper<'a>(&'a [usize], &'a BTreeMap<u8, usize>);
        impl Debug for AtomsHelper<'_> {
            fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
                f.debug_map()
                    .entries(self.0.iter().enumerate())
                    .entries(self.1)
                    .finish()
            }
        }
        f.debug_struct("EmpiricalFormula")
            .field("atoms", &AtomsHelper(&self.lower, &self.spill))
            .field("charge", &self.charge)
            .finish()
    }
}
impl Display for EmpiricalFormula {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
        struct ElemHelper(&'static str, usize);
        impl Display for ElemHelper {
            fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
                match self.1 {
                    0 => Ok(()),
                    1 => f.write_str(self.0),
                    n => {
                        f.write_str(self.0)?;
                        if f.alternate() {
                            Display::fmt(&Subscript(n), f)
                        } else {
                            Display::fmt(&n, f)
                        }
                    }
                }
            }
        }
        write!(
            f,
            "{}{}{}",
            ElemHelper("R", self.lower[0]),
            ElemHelper("C", self.lower[6]),
            ElemHelper("H", self.lower[1])
        )?;
        let mut segs = self
            .lower
            .iter()
            .enumerate()
            .filter(|(n, &c)| (c > 0 && ![0, 1, 6].contains(n)))
            .map(|(n, &c)| ElemHelper(ATOM_DATA[n].sym, c))
            .chain(
                self.spill
                    .iter()
                    .map(|(&n, &c)| ElemHelper(ATOM_DATA[n as usize].sym, c)),
            )
            .collect::<Vec<_>>();
        segs.sort();
        for s in &segs {
            Display::fmt(s, f)?;
        }
        if f.alternate() {
            match self.charge {
                0 => {}
                1 => f.write_str("⁺")?,
                -1 => f.write_str("⁻")?,
                _ => write!(f, "{:+}", Superscript(self.charge))?,
            }
        } else {
            match self.charge {
                0 => {}
                1 => f.write_str("+")?,
                -1 => f.write_str("-")?,
                _ => write!(f, "{:+}", self.charge)?,
            }
        }
        Ok(())
    }
}

impl std::str::FromStr for EmpiricalFormula {
    type Err = EmpiricalError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse_str(s)
    }
}

impl Extend<Atom> for EmpiricalFormula {
    fn extend<T: IntoIterator<Item = Atom>>(&mut self, iter: T) {
        iter.into_iter().for_each(|a| {
            self.add_atom(a.protons, 1);
            self.add_atom(1, a.data.hydrogen() as _);
            self.charge += a.charge;
        });
    }
}
impl<'a> Extend<&'a Atom> for EmpiricalFormula {
    fn extend<T: IntoIterator<Item = &'a Atom>>(&mut self, iter: T) {
        self.extend(iter.into_iter().copied());
    }
}
impl FromIterator<Atom> for EmpiricalFormula {
    fn from_iter<T: IntoIterator<Item = Atom>>(iter: T) -> Self {
        let mut out = Self::new();
        out.extend(iter);
        out
    }
}
impl<'a> FromIterator<&'a Atom> for EmpiricalFormula {
    fn from_iter<T: IntoIterator<Item = &'a Atom>>(iter: T) -> Self {
        let mut out = Self::new();
        out.extend(iter);
        out
    }
}

impl AddAssign<&Self> for EmpiricalFormula {
    fn add_assign(&mut self, other: &Self) {
        for (l, r) in self.lower.iter_mut().zip(&other.lower) {
            *l += r;
        }
        for (&n, &c) in &other.spill {
            self.spill.entry(n).and_modify(|a| *a += c).or_insert(c);
        }
        self.charge += other.charge;
    }
}
impl AddAssign for EmpiricalFormula {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}
impl Add<&Self> for EmpiricalFormula {
    type Output = Self;
    fn add(mut self, other: &Self) -> Self {
        self += other;
        self
    }
}
impl Add for EmpiricalFormula {
    type Output = Self;
    fn add(mut self, other: Self) -> Self {
        self += other;
        self
    }
}
impl Add for &EmpiricalFormula {
    type Output = EmpiricalFormula;
    fn add(self, other: Self) -> EmpiricalFormula {
        let mut this = self.clone();
        this += other;
        this
    }
}
