//! Complex molecular properties, i.e. something you probably don't want to recompute
use frunk::prelude::*;

pub trait Select<T> {
    fn get_prop(self) -> T;
}
