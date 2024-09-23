//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

pub mod access;
pub mod arena;
pub mod molecule;

pub use access::{ArenaAccessible, ArenaAccessibleMut, ArenaAccessor, ArenaAccessorMut};
pub use arena::Arena;
pub use molecule::Molecule;

use crate::core::*;
use petgraph::graph::IndexType;
use tracing::*;
