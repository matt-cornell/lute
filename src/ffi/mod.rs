#![allow(non_snake_case)]
//! FFI for this library. Probably shouldn't be called directly, only visible for documentation

use crate::core::MoleculeGraph as Graph;
use crate::prelude::{self as lib, Arena};

pub mod c;
