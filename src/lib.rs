#![allow(clippy::module_inception)]
#![warn(clippy::infinite_loop)]
#![feature(cmp_minmax, get_many_mut, maybe_uninit_write_slice)]

#[rustfmt::skip]
pub mod atom_info;
pub mod arena;
pub mod core;
pub mod disp;
pub mod empirical;
pub mod graph;
pub mod molecule;
pub mod parse;
pub mod utils;

#[cfg(feature = "rand")]
pub mod rand;

pub mod prelude {
    pub use crate::arena::{Arena, Container, Molecule};
    pub use crate::core::{Atom, Bond, Chirality};
    pub use crate::disp::fmt_as_dot;
    pub use crate::disp::smiles::*;
    #[cfg(feature = "coordgen")]
    pub use crate::disp::SvgFormatter;
    pub use crate::empirical::EmpiricalFormula;
    pub use crate::graph::misc::DataValueMap;
    pub use crate::molecule::Molecule as MolTrait;
    pub use crate::parse::smiles::SmilesParser;
    pub use crate::{arena, empirical, smiles};
}

#[cfg(test)]
mod tests;
