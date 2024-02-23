//! The purpose of this is to efficiently handle large groups of similar molecules. Because
//! reactions often only involve a small part of the molecule, it would be inefficient to make
//! copies of everything.

pub mod access;
pub mod arena;
pub mod container;
pub mod molecule;

pub use access::{ArenaAccessible, ArenaAccessibleMut, ArenaAccessor, ArenaAccessorMut};
pub use arena::Arena;
pub use container::Container;
pub use molecule::Molecule;

use petgraph::graph::IndexType;
use crate::molecule::*;
