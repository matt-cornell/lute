use crate::arena::molecule::NodeIndex;
use num_traits::*;
use petgraph::graph::IndexType;
use petgraph::visit::VisitMap;
use smallvec::SmallVec;
use std::fmt::{self, Binary, Debug, Formatter};

/// A `BitSet` acts as a set of `usize`, but is implemented using bits.
///
/// It wraps a `SmallVec` of an integer, but also can hold an offset to track a compact set of high-valued numbers efficiently.
#[derive(Default, Clone, PartialEq, Eq)]
pub struct BitSet<T, const N: usize> {
    bits: SmallVec<T, N>,
    offset: usize,
    minimum: usize,
}
impl<T: PrimInt, const N: usize> BitSet<T, N> {
    pub const fn new() -> Self {
        Self {
            bits: SmallVec::new(),
            offset: 0,
            minimum: usize::MAX,
        }
    }

    /// Get the bit at a specified index, defaulting to false.
    pub fn get(&self, idx: usize) -> bool {
        let bits = T::zero().leading_zeros() as usize;
        self.bits
            .get(idx / bits - self.offset)
            .map_or(false, |&i| i & (T::one() << (idx % bits)) != T::zero())
    }
    /// Set the bit at an index, returning the previous value.
    pub fn set(&mut self, idx: usize, bit: bool) -> bool {
        let zero = T::zero();
        let bits = zero.leading_zeros() as usize;

        let mut si = idx / bits;
        let sb = idx % bits;
        let mask = T::one() << sb;
        if bit {
            if si < self.offset {
                let l = self.bits.len();
                self.bits.resize(l + si, zero);
                self.bits.copy_within(..l, si);
                self.offset = si;
            } else if si >= self.bits.len() {
                self.bits.resize(si + 1, zero);
            }
            si -= self.offset;
            if idx < self.minimum || self.minimum == usize::MAX {
                self.minimum = idx;
            }
            let blk = self.bits[si];
            let ret = blk & mask != zero;
            self.bits[si] = blk | mask;
            ret
        } else {
            si -= self.offset;
            if idx == self.minimum {
                // start iterating at si since idx is the smallest element, everything below must
                // be zero
                if let Some((n, &blk)) = self.bits[si..]
                    .iter()
                    .enumerate()
                    .find(|(_, &blk)| blk != zero)
                {
                    self.minimum = (n + si + self.offset) * bits + blk.trailing_zeros() as usize;
                } else {
                    let ret = self.bits.get(si).map_or(false, |&b| b & mask != zero);
                    // no other bits were found, this set is empty
                    self.minimum = usize::MAX;
                    self.offset = 0;
                    self.bits.truncate(0);
                    return ret;
                }
            }
            if let Some(blk) = self.bits.get_mut(si) {
                let ret = *blk & mask != zero;
                *blk = *blk & !mask;
                ret
            } else {
                false
            }
        }
    }
    /// Clear all of the bits in this set.
    pub fn clear(&mut self) {
        self.bits.clear();
        self.minimum = usize::MAX;
    }

    /// Check if this set is "empty".
    ///
    /// This isn't necessarily a cheap operation since it doesn't store the size by itself.
    pub fn all_zero(&self) -> bool {
        let zero = T::zero();
        self.bits.iter().all(|&i| i == zero)
    }

    /// Count the number of set bits in this set.
    ///
    /// Just like `all_zero`, this isn't cheap because it doesn't store its size.
    pub fn count_ones(&self) -> usize {
        self.bits.iter().map(|i| i.count_ones() as usize).sum()
    }

    /// Get the highest bit set, or `None` if there are none.
    pub fn max_set(&self) -> Option<usize> {
        let mut i = self.bits.len() - 1;
        let zero = T::zero();
        while i > 0 && self.bits[i] == zero {
            i -= 1;
        }
        let v = self.bits[i];
        if v == zero {
            None
        } else {
            let bits = zero.leading_zeros() as usize;
            Some((i + self.offset) * bits + (bits - v.leading_zeros() as usize))
        }
    }

    /// Get the nth set bit.
    pub fn nth(&self, mut idx: usize) -> Option<usize> {
        use std::cmp::Ordering;
        idx += 1;
        let zero = T::zero();
        let one = T::one();
        let bits = zero.leading_zeros() as usize;
        for (i, &blk) in self.bits.iter().enumerate() {
            let ones = blk.count_ones() as usize;
            match idx.cmp(&ones) {
                Ordering::Greater => idx -= ones,
                Ordering::Equal => {
                    return Some((i + self.offset + 1) * bits - blk.leading_zeros() as usize - 1)
                }
                Ordering::Less => {
                    let mut j = 0;
                    loop {
                        debug_assert!(j < bits);
                        if blk & (one << j) != zero {
                            idx -= 1;
                            if idx == 0 {
                                break;
                            }
                        }
                        j += 1;
                    }
                    return Some((i + self.offset) * bits + j);
                }
            }
        }
        None
    }

    /// Get many bits by there indices.
    ///
    /// This method could be optimized to only loop over the underlying vector once, but as of right now it does it separately.
    pub fn nth_many<const O: usize>(&self, idx: [usize; O]) -> [Option<usize>; O] {
        // TODO: optimize
        idx.map(|i| self.nth(i))
    }

    /// Get the index of a set bit, returning `None` if out of range or unset.
    ///
    /// This is the inverse operation of `nth`.
    pub fn index(&self, bit: usize) -> Option<usize> {
        let zero = T::zero();
        let one = T::one();
        let bits = zero.leading_zeros() as usize;
        let bit = bit.checked_sub(self.offset * bits)?;
        let prev: usize = self
            .bits
            .get(..(bit / bits))?
            .iter()
            .map(|blk| blk.count_ones() as usize)
            .sum();
        let blk = *self.bits.get(bit / bits)?;
        let mask = one << (bit % bits);
        (blk & mask != zero).then_some(())?;
        let curr = (blk & (mask - one)).count_ones() as usize;
        Some(prev + curr)
    }

    /// Iterate over the set bits by reference.
    pub fn iter(&self) -> Iter<'_, T, N> {
        Iter {
            bits: self,
            next: 0,
        }
    }
}

impl<T: Binary, const N: usize> Debug for BitSet<T, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut l = f.debug_list();
        for i in &self.bits {
            l.entry(&format_args!("{i:0>0$b}", std::mem::size_of::<T>() * 8));
        }
        l.finish()
    }
}

impl<T: PrimInt, const N: usize> Extend<usize> for BitSet<T, N> {
    fn extend<I: IntoIterator<Item = usize>>(&mut self, iter: I) {
        for i in iter {
            self.set(i, true);
        }
    }
}
impl<T: PrimInt + Default, const N: usize> FromIterator<usize> for BitSet<T, N> {
    fn from_iter<I: IntoIterator<Item = usize>>(iter: I) -> Self {
        let mut this = Self::default();
        this.extend(iter);
        this
    }
}

impl<T: PrimInt, const N: usize> IntoIterator for BitSet<T, N> {
    type IntoIter = IntoIter<T, N>;
    type Item = usize;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter {
            bits: self,
            next: 0,
        }
    }
}

impl<'a, T: PrimInt, const N: usize> IntoIterator for &'a BitSet<T, N> {
    type IntoIter = Iter<'a, T, N>;
    type Item = usize;

    fn into_iter(self) -> Self::IntoIter {
        Iter {
            bits: self,
            next: 0,
        }
    }
}

impl<Ix: IndexType, T: PrimInt, const N: usize> VisitMap<NodeIndex<Ix>> for BitSet<T, N> {
    fn visit(&mut self, a: NodeIndex<Ix>) -> bool {
        self.set(a.0.index(), true)
    }
    fn is_visited(&self, a: &NodeIndex<Ix>) -> bool {
        self.get(a.0.index())
    }
}

/// Common iterator logic, since `Iter` and `IntoIter` do the same thing.
fn get_next<T: PrimInt, const N: usize>(bits: &BitSet<T, N>, next: &mut usize) -> Option<usize> {
    let zero = T::zero();
    let one = T::one();
    let nbits = zero.leading_zeros() as usize;
    let idx = *next / nbits;
    let word = *bits.bits.get(idx)?;
    let zeros = (word & !((one << (*next % nbits)) - one)).trailing_zeros() as usize;
    if zeros < nbits {
        let res = idx * nbits + zeros;
        *next = res + 1;
        return Some(bits.offset + res);
    }
    let Some((idx, &word)) = bits.bits[(idx + 1)..]
        .iter()
        .enumerate()
        .find(|(_, n)| **n != zero)
    else {
        *next = usize::MAX;
        return None;
    };
    let zeros = (word & !((one << (*next % nbits)) - one)).trailing_zeros() as usize;
    let res = idx * nbits + zeros;
    *next = res + 1;
    Some(bits.offset + res)
}

#[derive(Clone)]
pub struct Iter<'a, T, const N: usize> {
    bits: &'a BitSet<T, N>,
    next: usize,
}
impl<'a, T: PrimInt, const N: usize> Iterator for Iter<'a, T, N> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        get_next(self.bits, &mut self.next)
    }
}

#[derive(Clone)]
pub struct IntoIter<T, const N: usize> {
    bits: BitSet<T, N>,
    next: usize,
}
impl<T: PrimInt, const N: usize> Iterator for IntoIter<T, N> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        get_next(&self.bits, &mut self.next)
    }
}
