//! Utilities for modeling molecules as graphs
// This top-level file handles definitions and parsing

use crate::atom_info::ATOM_DATA;
use c_enum::*;
use petgraph::prelude::*;
use std::fmt::{self, Display, Formatter};
use std::hash::{Hash, Hasher};
use modular_bitfield::prelude::*;
use thiserror::Error;

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

c_enum! {
    #[derive(Clone, Copy, PartialEq, Eq, Hash)]
    #[repr(transparent)]
    pub enum Chirality: u8 {
        None,
        Ccw,
        Cw,
    }
}
impl Chirality {
    pub fn is_chiral(self) -> bool {
        self != Self::None
    }
}
#[bitfield]
#[repr(u16)]
#[derive(Debug, Clone, Copy)]
pub struct AtomData {
    pub hydrogen: B4,
    pub unknown: B4,
    pub other: B4,
    pub scratch: B4,
}

/// An atom in the molecule graph
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct Atom {
    pub protons: u8,
    pub charge: i8,
    pub isotope: u16,
    pub chirality: Chirality,
    pub data: AtomData,
}
impl Atom {
    pub fn new(protons: u8) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            chirality: Chirality::None,
            data: AtomData::new(),
        }
    }
    pub fn new_scratch(protons: u8, scratch: u8) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            chirality: Chirality::None,
            data: AtomData::new().with_scratch(scratch),
        }
    }
    pub fn new_isotope(protons: u8, isotope: u16) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            chirality: Chirality::None,
            data: AtomData::new(),
        }
    }
    pub fn new_isotope_scratch(protons: u8, isotope: u16, scratch: u8) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            chirality: Chirality::None,
            data: AtomData::new().with_scratch(scratch),
        }
    }
    
    pub fn add_hydrogens(&mut self, h: u8) -> Result<(), TooManyBonds> {
        let h = self.data.hydrogen() + h;
        if h < 16 {
            self.data.set_hydrogen(h);
            Ok(())
        }
        else {
            Err(TooManyBonds(TooMany::H, h as _))
        }
    }
    pub fn add_rs(&mut self, r: u8) -> Result<(), TooManyBonds> {
        let h = self.data.unknown() + r;
        if h < 16 {
            self.data.set_unknown(r);
            Ok(())
        }
        else {
            Err(TooManyBonds(TooMany::R, r as _))
        }
    }
    pub fn set_other_bonds(&mut self, b: u8) -> Result<(), TooManyBonds> {
        if b < 16 {
            self.data.set_other(b);
            Ok(())
        }
        else {
            Err(TooManyBonds(TooMany::Other, b as _))
        }
    }
    
    #[inline(always)]
    pub fn map_scratch<F: FnOnce(u8) -> u8>(&mut self, f: F) {
        self.data.set_scratch(f(self.data.scratch()));
    }
    #[inline(always)]
    pub fn with_scratch<R, F: FnOnce(&mut u8) -> R>(&mut self, f: F) -> R {
        let mut scratch = self.data.scratch();
        let r = f(&mut scratch);
        self.data.set_scratch(scratch);
        r
    }

    pub fn mass(self) -> f32 {
        ATOM_DATA[self.protons as usize].mass
    }

    pub fn eq_or_r(&self, other: &Self) -> bool {
        self.protons == 0 || other.protons == 0 || self.eq_match_r(other)
    }
    pub fn eq_match_r(&self, other: &Self) -> bool {
        (self.protons == 0 && other.protons == 0) || self == other
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
            if self.protons == 0 {
                write!(f, "{}", Subscript(self.isotope))?;
            }
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
impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.protons == other.protons
            && self.isotope == other.isotope
            && self.charge == other.charge
            && (self.chirality == Chirality::None
                || other.chirality == Chirality::None
                || self.chirality == other.chirality)
    }
}
impl Eq for Atom {}
impl Hash for Atom {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.protons.hash(state);
        self.isotope.hash(state);
        self.charge.hash(state);
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
        Left,
        Right,
    }
}
impl Bond {
    pub fn bond_count(self) -> f32 {
        match self {
            Self::Single | Self::Left | Self::Right => 1f32,
            Self::Double => 2f32,
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
