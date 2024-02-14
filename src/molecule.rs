//! Utilities for modeling molecules as graphs
// This top-level file handles definitions and parsing

use crate::atom_info::ATOM_DATA;
use c_enum::*;
use petgraph::prelude::*;
use std::fmt::{self, Display, Formatter};
use std::hash::{Hash, Hasher};

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

/// An atom in the molecule graph
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct Atom {
    pub protons: u8,
    pub charge: i8,
    pub isotope: u16,
    pub chirality: Chirality,
    pub hydrogens: u8,
    /// Scratch buffer for use in a parser, ignored by all equality checks. Probably shouldn't be
    /// used after parsing
    pub scratch: u16,
}
impl Atom {
    pub fn new(protons: u8) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            chirality: Chirality::None,
            hydrogens: 0,
            scratch: 0,
        }
    }
    pub fn new_scratch(protons: u8, scratch: u16) -> Self {
        Self {
            protons,
            isotope: 0,
            charge: 0,
            chirality: Chirality::None,
            hydrogens: 0,
            scratch,
        }
    }
    pub fn new_isotope(protons: u8, isotope: u16) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            chirality: Chirality::None,
            hydrogens: 0,
            scratch: 0,
        }
    }
    pub fn new_isotope_scratch(protons: u8, isotope: u16, scratch: u16) -> Self {
        Self {
            protons,
            isotope,
            charge: 0,
            chirality: Chirality::None,
            hydrogens: 0,
            scratch,
        }
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
