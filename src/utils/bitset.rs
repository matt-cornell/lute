use num_traits::*;
use smallvec::SmallVec;
use std::fmt::{self, Binary, Debug, Formatter};

#[derive(Clone, PartialEq, Eq)]
pub struct BitSet<T, const N: usize>(SmallVec<T, N>);
impl<T: PrimInt + Zero, const N: usize> BitSet<T, N> {
    pub const fn new() -> Self {
        Self(SmallVec::new())
    }
    pub fn with_capacity(cap: usize) -> Self {
        let len = (cap + std::mem::size_of::<T>() - 1) / std::mem::size_of::<T>();
        Self(SmallVec::from_elem(T::zero(), len))
    }

    pub fn as_slice(&self) -> &[T] {
        self.0.as_slice()
    }

    pub fn get(&self, idx: usize) -> bool {
        self.0
            .get(idx / std::mem::size_of::<T>())
            .map_or(false, |&i| {
                i & (T::one() << (idx % std::mem::size_of::<T>())) != T::zero()
            })
    }
    pub fn set(&mut self, idx: usize, bit: bool) {
        let si = idx / std::mem::size_of::<T>();
        let sb = idx % std::mem::size_of::<T>();
        if si >= self.0.len() {
            self.0.resize(si + 1, T::zero());
        }
        let word = self.0[si];
        self.0[si] = if bit {
            word | (T::one() << sb)
        } else {
            word & !(T::one() << sb)
        };
    }
    pub fn clear(&mut self) {
        for i in &mut self.0 {
            *i = T::zero();
        }
    }
}

impl<T: Binary, const N: usize> Debug for BitSet<T, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut l = f.debug_list();
        for i in &self.0 {
            l.entry(&format_args!("{i:0>8b}"));
        }
        l.finish()
    }
}
