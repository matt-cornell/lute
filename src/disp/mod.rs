#[cfg(feature = "coordgen")]
pub mod cgen;
pub mod dot;
pub mod svg;

mod color;

#[cfg(feature = "coordgen")]
pub use cgen::fmt_with_cg;
pub use dot::fmt_as_dot;
