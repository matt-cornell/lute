use itertools::{EitherOrBoth::*, Itertools};
use num_traits::*;
use smallvec::SmallVec;
use std::fmt::{self, Binary, Debug, Formatter};
use std::ops::*;

#[derive(Default, Clone, PartialEq, Eq)]
pub struct BitSet<T, const N: usize>(SmallVec<T, N>);
impl<T: PrimInt + Zero, const N: usize> BitSet<T, N> {
    pub const fn new() -> Self {
        Self(SmallVec::new())
    }
    pub fn with_capacity(cap: usize) -> Self {
        let bits = T::zero().leading_zeros() as usize;
        let len = (cap + bits - 1) / bits;
        Self(SmallVec::from_elem(T::zero(), len))
    }
    pub fn from_buf(buf: SmallVec<T, N>) -> Self {
        Self(buf)
    }

    pub fn as_slice(&self) -> &[T] {
        self.0.as_slice()
    }

    pub fn get(&self, idx: usize) -> bool {
        let bits = T::zero().leading_zeros() as usize;
        self.0
            .get(idx / bits)
            .map_or(false, |&i| i & (T::one() << (idx % bits)) != T::zero())
    }
    pub fn set(&mut self, idx: usize, bit: bool) {
        let bits = T::zero().leading_zeros() as usize;
        let si = idx / bits;
        let sb = idx % bits;
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

    pub fn all_zero(&self) -> bool {
        let zero = T::zero();
        self.0.iter().all(|&i| i == zero)
    }
    pub fn all_ones(&self) -> bool {
        let ones = !T::zero();
        self.0.iter().all(|&i| i == ones)
    }
    pub fn count_ones(&self) -> usize {
        self.0.iter().map(|i| i.count_ones() as usize).sum()
    }
    pub fn max_set(&self) -> Option<usize> {
        let mut i = self.0.len() - 1;
        let zero = T::zero();
        while i > 0 && self.0[i] == zero {
            i -= 1;
        }
        let v = self.0[i];
        if v == zero {
            None
        } else {
            let bits = zero.leading_zeros() as usize;
            Some(i * bits + (bits - v.leading_zeros() as usize))
        }
    }
    pub fn nth(&self, mut idx: usize) -> Option<usize> {
        let zero = T::zero();
        let one = T::one();
        let bits = zero.leading_zeros() as usize;
        for (i, &blk) in self.0.iter().enumerate() {
            let ones = blk.count_ones();
            if ones == 0 {
                continue;
            }
            if let Some(i2) = idx.checked_sub(ones as usize) {
                if i2 > 0 {
                    idx = i2;
                    continue;
                }
            }
            let mut j = 0;
            while j < bits && idx > 0 {
                idx -= (blk & (one << j) != zero) as usize;
                j += 1;
            }
            debug_assert_eq!(idx, 0);
            return Some(i * bits + j);
        }
        None
    }
    #[deprecated(note = "doesn't seem to work")]
    pub fn nth_many<const O: usize>(&self, mut idx: [usize; O]) -> [Option<usize>; O] {
        let mut out = [None; O];
        let zero = T::zero();
        let one = T::one();
        let bits = zero.leading_zeros() as usize;
        // this is probably an over-allocation
        let mut buf = SmallVec::<usize, { std::mem::size_of::<usize>() * 8 }>::new();
        for (i, &blk) in self.0.iter().enumerate() {
            buf.clear();
            let ones = blk.count_ones() as usize;
            let mut can_ret = true;
            for (n, idx) in idx.iter_mut().enumerate() {
                if let Some(i2) = idx.checked_sub(ones) {
                    if i2 > 0 {
                        *idx = i2;
                        can_ret = false;
                        continue;
                    }
                }
                if *idx == 0 {
                    out[n] = Some(i * bits);
                } else if let Some(&j) = buf.get(*idx - 1) {
                    *idx = 0;
                    out[n] = Some(i * bits + j);
                } else {
                    let mut j = buf.last().map_or(0, |&j| j);
                    *idx -= buf.len();
                    while j < bits && *idx > 0 {
                        if blk & (one << j) != zero {
                            buf.push(j);
                            *idx -= 1;
                        }
                        j += 1;
                    }
                    debug_assert_eq!(*idx, 0);
                    out[n] = Some(i * bits + j);
                }
            }
            if can_ret {
                return out;
            }
        }
        out
    }
    #[deprecated(note = "doesn't seem to work")]
    pub fn nth_many_short<const O: usize>(&self, mut idx: [usize; O]) -> Option<[usize; O]> {
        let mut out = [0; O];
        let zero = T::zero();
        let one = T::one();
        let bits = zero.leading_zeros() as usize;
        // this is probably an over-allocation
        let mut buf = SmallVec::<usize, { std::mem::size_of::<usize>() * 8 }>::new();
        for (i, &blk) in self.0.iter().enumerate() {
            buf.clear();
            let ones = blk.count_ones() as usize;
            let mut can_ret = true;
            for (n, idx) in idx.iter_mut().enumerate() {
                if let Some(i2) = idx.checked_sub(ones) {
                    if i2 > 0 {
                        *idx = i2;
                        can_ret = false;
                        continue;
                    }
                }
                if *idx == 0 {
                    out[n] = i * bits;
                } else if let Some(&j) = buf.get(*idx - 1) {
                    *idx = 0;
                    out[n] = i * bits + j;
                } else {
                    let mut j = buf.last().map_or(0, |&j| j);
                    *idx -= buf.len();
                    while j < bits && *idx > 0 {
                        if blk & (one << j) != zero {
                            buf.push(j);
                            *idx -= 1;
                        }
                        j += 1;
                    }
                    debug_assert_eq!(*idx, 0);
                    out[n] = i * bits + j;
                }
            }
            if can_ret {
                return Some(out);
            }
        }
        None
    }

    pub fn index(&self, bit: usize) -> Option<usize> {
        let zero = T::zero();
        let one = T::one();
        let bits = T::zero().leading_zeros() as usize;
        let prev: usize = self
            .0
            .get(..(bit / bits))?
            .iter()
            .map(|blk| blk.count_ones() as usize)
            .sum();
        let blk = *self.0.get(bit / bits)?;
        let mask = one << (bit % bits);
        (blk & mask != zero).then_some(())?;
        let curr = (blk & (mask - one)).count_ones() as usize;
        Some(prev + curr)
    }
}

impl<T: Binary, const N: usize> Debug for BitSet<T, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut l = f.debug_list();
        for i in &self.0 {
            l.entry(&format_args!("{i:0>0$b}", std::mem::size_of::<T>() * 8));
        }
        l.finish()
    }
}

impl<T: PrimInt, const N: usize> BitAnd for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(mut self, rhs: Self) -> BitSet<T, N> {
        self &= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitAnd<&Self> for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(mut self, rhs: &Self) -> BitSet<T, N> {
        self &= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitAnd for &BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(self, rhs: Self) -> BitSet<T, N> {
        BitSet(self.0.iter().zip(&rhs.0).map(|(l, r)| *l & *r).collect())
    }
}

impl<T: PrimInt, const N: usize> BitAndAssign for BitSet<T, N> {
    fn bitand_assign(&mut self, rhs: Self) {
        *self &= &rhs;
    }
}

impl<T: PrimInt, const N: usize> BitAndAssign<&Self> for BitSet<T, N> {
    fn bitand_assign(&mut self, rhs: &Self) {
        self.0
            .iter_mut()
            .zip(&rhs.0)
            .for_each(|(l, r)| *l = *l & *r);
    }
}

impl<T: PrimInt, const N: usize> BitOr for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitor(mut self, rhs: Self) -> BitSet<T, N> {
        self |= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitOr<&Self> for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitor(mut self, rhs: &Self) -> BitSet<T, N> {
        self |= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitOr for &BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitor(self, rhs: Self) -> BitSet<T, N> {
        BitSet(
            self.0
                .iter()
                .zip_longest(&rhs.0)
                .map(|e| match e {
                    Left(i) | Right(i) => *i,
                    Both(l, r) => *l | *r,
                })
                .collect(),
        )
    }
}

impl<T: PrimInt, const N: usize> BitOrAssign for BitSet<T, N> {
    fn bitor_assign(&mut self, rhs: Self) {
        *self |= &rhs;
    }
}

impl<T: PrimInt, const N: usize> BitOrAssign<&Self> for BitSet<T, N> {
    fn bitor_assign(&mut self, rhs: &Self) {
        let mut iter = rhs.0.iter();
        self.0
            .iter_mut()
            .zip(iter.by_ref())
            .for_each(|(l, r)| *l = *l | *r);
        self.0.extend(iter.copied());
    }
}

impl<T: PrimInt, const N: usize> BitXor for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitxor(mut self, rhs: Self) -> BitSet<T, N> {
        self ^= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitXor<&Self> for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitxor(mut self, rhs: &Self) -> BitSet<T, N> {
        self ^= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitXor for &BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitxor(self, rhs: Self) -> BitSet<T, N> {
        BitSet(
            self.0
                .iter()
                .zip_longest(&rhs.0)
                .map(|e| match e {
                    Left(i) | Right(i) => *i,
                    Both(l, r) => *l ^ *r,
                })
                .collect(),
        )
    }
}

impl<T: PrimInt, const N: usize> BitXorAssign for BitSet<T, N> {
    fn bitxor_assign(&mut self, rhs: Self) {
        *self ^= &rhs;
    }
}

impl<T: PrimInt, const N: usize> BitXorAssign<&Self> for BitSet<T, N> {
    fn bitxor_assign(&mut self, rhs: &Self) {
        let mut iter = rhs.0.iter();
        self.0
            .iter_mut()
            .zip(iter.by_ref())
            .for_each(|(l, r)| *l = *l ^ *r);
        self.0.extend(iter.copied());
    }
}

impl<'a, T, const N: usize> Not for &'a BitSet<T, N> {
    type Output = InvSet<'a, T>;

    fn not(self) -> InvSet<'a, T> {
        InvSet(&self.0)
    }
}

/// Since `BitSet` is a set, we can't construct its complement in a finite amount of space.
/// Instead, we have an `InvSet` which does basically the same thing so operations can be used.
#[derive(Debug, Clone, Copy)]
pub struct InvSet<'a, T>(&'a [T]);

impl<T: PrimInt, const N: usize> BitAnd<InvSet<'_, T>> for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(mut self, rhs: InvSet<'_, T>) -> BitSet<T, N> {
        self &= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitAnd<&InvSet<'_, T>> for BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(mut self, rhs: &InvSet<'_, T>) -> BitSet<T, N> {
        self &= rhs;
        self
    }
}

impl<T: PrimInt, const N: usize> BitAnd<&InvSet<'_, T>> for &BitSet<T, N> {
    type Output = BitSet<T, N>;

    fn bitand(self, rhs: &InvSet<'_, T>) -> BitSet<T, N> {
        BitSet(self.0.iter().zip(rhs.0).map(|(l, r)| *l & !*r).collect())
    }
}

impl<T: PrimInt, const N: usize> BitAndAssign<InvSet<'_, T>> for BitSet<T, N> {
    fn bitand_assign(&mut self, rhs: InvSet<'_, T>) {
        *self &= &rhs;
    }
}

impl<T: PrimInt, const N: usize> BitAndAssign<&InvSet<'_, T>> for BitSet<T, N> {
    fn bitand_assign(&mut self, rhs: &InvSet<'_, T>) {
        self.0
            .iter_mut()
            .zip(rhs.0)
            .for_each(|(l, r)| *l = *l & !*r);
    }
}
