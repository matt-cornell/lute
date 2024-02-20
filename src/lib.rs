#[rustfmt::skip]
pub mod atom_info;
pub mod disp;
pub mod graph;
pub mod molecule;
pub mod parse;
pub mod utils;

pub mod prelude {
    pub use crate::disp::fmt_as_dot;
    #[cfg(feature = "mol-svg")]
    pub use crate::disp::fmt_as_svg;
    #[cfg(feature = "mol-bmp")]
    pub use crate::disp::{make_img, make_img_vec};
    pub use crate::graph::arena::{Arena, Container, Molecule};
    pub use crate::graph::compact::GraphCompactor;
    pub use crate::molecule::{Atom, Bond};
    pub use crate::parse::smiles::SmilesParser;
}
