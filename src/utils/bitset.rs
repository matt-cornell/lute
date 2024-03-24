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
        let len = (cap + std::mem::size_of::<T>() - 1) / std::mem::size_of::<T>();
        Self(SmallVec::from_elem(T::zero(), len))
    }
    pub fn from_buf(buf: SmallVec<T, N>) -> Self {
        Self(buf)
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
