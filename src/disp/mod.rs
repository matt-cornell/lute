pub mod dot;
pub mod img;
pub mod svg;

pub use dot::fmt_as_dot;
pub use img::{make_img, make_img_vec};
pub use svg::fmt_as_svg;
