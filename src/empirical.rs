use crate::atom_info::*;
use crate::core::*;
use fmtastic::*;
use std::collections::BTreeMap;
use std::fmt::{self, Debug, Display, Formatter};
use std::ops::{Add, AddAssign};

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
    pub fn get(&mut self, atom: u8) -> usize {
        if atom < 18 {
            self.lower[atom as usize]
        } else {
            self.spill.get(&atom).copied().unwrap_or(0)
        }
    }
    pub fn set(&mut self, atom: u8, count: usize) {
        if atom < 18 {
            self.lower[atom as usize] = count;
        } else if count == 0 {
            self.spill.remove(&atom);
        } else {
            self.spill.insert(atom, count);
        }
    }
    pub fn add(&mut self, atom: u8, count: usize) {
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
        if f.alternate() {
            if self.lower[0] > 0 {
                write!(f, "R{}", Subscript(self.lower[0]))?;
            }
            if self.lower[6] > 0 {
                write!(f, "C{}", Subscript(self.lower[6]))?;
            }
            if self.lower[1] > 0 {
                write!(f, "H{}", Subscript(self.lower[1]))?;
            }
            let mut segs = self
                .lower
                .iter()
                .enumerate()
                .filter_map(|(n, &c)| {
                    (c > 0 && ![0, 1, 6].contains(&n))
                        .then(|| format!("{}{}", ATOM_DATA[n].sym, Subscript(c)))
                })
                .chain(
                    self.spill
                        .iter()
                        .map(|(&n, &c)| format!("{}{}", ATOM_DATA[n as usize].sym, Subscript(c))),
                )
                .collect::<Vec<_>>();
            segs.sort();
            f.write_str(&segs.join(""))?;
            match self.charge {
                0 => {}
                1 => f.write_str("⁺")?,
                -1 => f.write_str("⁻")?,
                _ => write!(f, "{:+}", Superscript(self.charge))?,
            }
        } else {
            if self.lower[0] > 0 {
                write!(f, "R{}", self.lower[0])?;
            }
            if self.lower[6] > 0 {
                write!(f, "C{}", self.lower[6])?;
            }
            if self.lower[1] > 0 {
                write!(f, "H{}", self.lower[1])?;
            }
            let mut segs = self
                .lower
                .iter()
                .enumerate()
                .filter_map(|(n, &c)| {
                    (c > 0 && ![0, 1, 6].contains(&n))
                        .then(|| format!("{}{c}", ATOM_DATA[n].sym))
                })
                .chain(
                    self.spill
                        .iter()
                        .map(|(&n, c)| format!("{}{c}", ATOM_DATA[n as usize].sym)),
                )
                .collect::<Vec<_>>();
            segs.sort();
            f.write_str(&segs.join(""))?;
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
impl Extend<Atom> for EmpiricalFormula {
    fn extend<T: IntoIterator<Item = Atom>>(&mut self, iter: T) {
        iter.into_iter().for_each(|a| {
            self.add(a.protons, 1);
            self.add(1, a.data.hydrogen());
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
