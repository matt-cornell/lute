#![allow(clippy::identity_op, clippy::excessive_precision)]
//! Utilities for modeling molecules as graphs
// This top-level file handles definitions and parsing

use crate::atom_info::ATOM_DATA;
use c_enum::*;
use modular_bitfield::prelude::*;
use petgraph::prelude::*;
use std::fmt::{self, Display, Formatter};
use std::hash::{Hash, Hasher};
use thiserror::Error;

const PROTON_MASS: f32 = 1.00727647;
const NEUTRON_MASS: f32 = 1.008665;
const ELECTRON_MASS: f32 = 0.00054858;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TooMany {
    H,
    R,
    Other,
}
impl TooMany {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::H => "hydrogen",
            Self::R => "unknown",
            Self::Other => "other",
        }
    }
}
impl Display for TooMany {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Error)]
#[error("too many {0} bonds: attempted to set {1}, the max is 16")]
pub struct TooManyBonds(pub TooMany, pub usize);

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, BitfieldSpecifier)]
pub enum Chirality {
    #[default]
    None,
    R,
    S,
    /// Shouldn't ever be used
    Error,
}
impl Chirality {
    pub fn is_chiral(self) -> bool {
        self != Self::None
    }
    pub fn is_valid(self) -> bool {
        self != Self::Error
    }
}

/// Bit-packed field for tracking bonds, plus two scratch bits
#[allow(clippy::identity_op)]
#[bitfield]
#[derive(Debug, Default, Clone, Copy)]
pub struct AtomData {
    pub hydrogen: B4,
    pub unknown: B4,
    pub other: B4,
    pub chirality: Chirality,
    pub scratch: B10,
}
impl AtomData {
    /// Count the total number of bonded atoms
    #[inline(always)]
    pub fn total_bonds(self) -> u8 {
        self.hydrogen() + self.unknown() + self.other()
    }
    /// Check if two bonding layouts are compatible
    pub fn is_compatible(self, other: Self) -> bool {
        let h1 = self.hydrogen();
        let o1 = self.other();
        let mut u1 = self.unknown();
        let h2 = other.hydrogen();
        let o2 = other.other();
        let mut u2 = other.unknown();
        if h1 + u1 + o1 != h2 + u2 + o2 {
            return false;
        }
        if h1 > h2 {
            let diff = h1 - h2;
            if u2 >= diff {
                u2 -= diff;
            } else {
                return false;
            }
        } else {
            let diff = h2 - h1;
            if u1 >= diff {
                u1 -= diff;
            } else {
                return false;
            }
        }
        if o1 > o2 {
            let diff = o1 - o2;
            if u2 >= diff {
                u2 -= diff;
            } else {
                return false;
            }
        } else {
            let diff = o2 - o1;
            if u1 >= diff {
                u1 -= diff;
            } else {
                return false;
            }
        }
        debug_assert_eq!(u1, u2, "u1 == u2 should hold because their sums are equal!");
        true
    }
}
impl PartialEq for AtomData {
    fn eq(&self, other: &Self) -> bool {
        self.chirality() == other.chirality() && self.is_compatible(*other)
    }
}
impl Eq for AtomData {}
impl Hash for AtomData {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u8(self.total_bonds());
    }
}

/// An atom in the molecule graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct Atom {
    pub protons: u8,
    pub charge: i8,
    pub isotope: u16,
    pub data: AtomData,
}
impl Atom {
    pub fn new(protons: u8) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            data: AtomData::new(),
        }
    }
    pub fn new_scratch(protons: u8, scratch: u16) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            data: AtomData::new().with_scratch(scratch),
        }
    }
    pub fn new_isotope(protons: u8, isotope: u16) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            data: AtomData::new(),
        }
    }
    pub fn new_isotope_scratch(protons: u8, isotope: u16, scratch: u16) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            data: AtomData::new().with_scratch(scratch),
        }
    }

    pub fn add_hydrogens(&mut self, h: u8) -> Result<(), TooManyBonds> {
        let h = self.data.hydrogen() + h;
        if h < 16 {
            self.data.set_hydrogen(h);
            Ok(())
        } else {
            Err(TooManyBonds(TooMany::H, h as _))
        }
    }
    pub fn add_rs(&mut self, r: u8) -> Result<(), TooManyBonds> {
        let r = self.data.unknown() + r;
        if r < 16 {
            self.data.set_unknown(r);
            Ok(())
        } else {
            Err(TooManyBonds(TooMany::R, r as _))
        }
    }
    pub fn set_other_bonds(&mut self, b: u8) -> Result<(), TooManyBonds> {
        if b < 16 {
            self.data.set_other(b);
            Ok(())
        } else {
            Err(TooManyBonds(TooMany::Other, b as _))
        }
    }

    #[inline(always)]
    pub fn map_scratch<F: FnOnce(u16) -> u16>(&mut self, f: F) {
        self.data.set_scratch(f(self.data.scratch()));
    }
    #[inline(always)]
    pub fn with_scratch<R, F: FnOnce(&mut u16) -> R>(&mut self, f: F) -> R {
        let mut scratch = self.data.scratch();
        let r = f(&mut scratch);
        self.data.set_scratch(scratch);
        r
    }

    pub fn mass(self) -> f32 {
        let base = if self.isotope == 0 {
            ATOM_DATA[self.protons as usize].mass
        } else {
            PROTON_MASS * self.protons as f32
                + NEUTRON_MASS * (self.isotope as f32 - self.protons as f32)
        };
        base - ELECTRON_MASS * self.charge as f32
    }

    /// The first atom, used in a pattern, could refer to the second.
    /// If the first atom has no isotope or chirality and the second does, this can still return true.
    /// Otherwise, acts the same as `eq_or_r`
    pub fn matches(&self, other: &Self) -> bool {
        if !self.data.is_compatible(other.data) {
            return false;
        }
        if self.protons == 0 {
            match other.protons {
                // this is a superset of the other
                0 => !self.isotope & other.isotope == 0,
                // bit 0 is hydrogen match
                1 => self.isotope & 1 != 0,
                // bit 1 is carbon match
                6 => self.isotope & 2 != 0,
                // the next bits are based on the ordering in atom_info.rs
                _ => self.isotope & (1 << (ATOM_DATA[other.protons as usize].group as u8 + 2)) != 0,
            }
        } else {
            self.protons == other.protons
                && self.charge == other.charge
                && (self.isotope == 0 || self.isotope == other.isotope)
                && (self.data.chirality() == Chirality::None
                    || self.data.chirality() == other.data.chirality())
        }
    }
}
impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            write!(f, "{}", ATOM_DATA[self.protons as usize].name)?;
            if self.isotope != 0 || self.protons == 0 {
                write!(f, "-{}", self.isotope)?;
            }
        } else {
            use fmtastic::*;
            if self.protons != 0 && self.isotope != 0 {
                write!(f, "{}", Superscript(self.isotope))?;
            }
            write!(f, "{}", ATOM_DATA[self.protons as usize].sym)?;
            match self.charge {
                0 => {}
                1 => f.write_str("⁺")?,
                -1 => f.write_str("⁻")?,
                _ => write!(f, "{:+}", Superscript(self.charge))?,
            }
        }
        Ok(())
    }
}

c_enum! {
    /// A bond between atoms in the molecule graph
    #[derive(Clone, Copy, PartialEq, Eq)]
    pub enum Bond: u8 {
        /// Non-bond, shouldn't appear in final graph
        Non,
        Single,
        Double,
        Triple,
        Quad,
        Aromatic,
        DoubleE,
        DoubleZ,
        Left,
        Right,
    }
}
impl Bond {
    pub fn bond_count(self) -> f32 {
        match self {
            Self::Single | Self::Left | Self::Right => 1f32,
            Self::Double | Self::DoubleE | Self::DoubleZ => 2f32,
            Self::Triple => 3f32,
            Self::Quad => 4f32,
            Self::Aromatic => 1.5f32,
            Self::Non => 0f32,
            _ => panic!("invalid bond!"),
        }
    }
    pub fn as_static_str(self) -> &'static str {
        match self {
            Self::Non => "non",
            Self::Single => "single",
            Self::Double => "double",
            Self::DoubleE => "trans double",
            Self::DoubleZ => "cis double",
            Self::Triple => "triple",
            Self::Quad => "quad",
            Self::Aromatic => "aromatic",
            Self::Left => "left",
            Self::Right => "right",
            _ => panic!("invalid bond!"),
        }
    }
}
impl Display for Bond {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_static_str())
    }
}
/// A molecule graph is an undirected graph between atoms, connected with bonds
pub type MoleculeGraph = UnGraph<Atom, Bond>;
