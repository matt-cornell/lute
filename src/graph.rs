//! Utilities for the molecule graph

use slab::Slab;
use std::cell::UnsafeCell;
use std::cmp::Ordering;
use std::fmt::{self, Debug, Formatter};
use std::iter::{ExactSizeIterator, FusedIterator};

#[ouroboros::self_referencing]
struct SortedSlabInner<T: 'static> {
    /// Workaround because ouroboros doesn't let me get a reference to this and I can make sure
    /// that references are never *invalidated*, therefore it's all safe
    backing: UnsafeCell<Slab<T>>,
    #[borrows(backing)]
    #[covariant]
    sorted: copse::BTreeSet<usize, SlabOrdering<'this, T>>,
}

/// This container acts like `Slab` for most purposes, but also can be iterated through in a sorted order.  
pub struct SortedSlab<T: 'static> {
    inner: SortedSlabInner<T>,
}
impl<T: 'static> SortedSlab<T> {
    pub fn new() -> Self {
        Self {
            inner: SortedSlabInner::new(UnsafeCell::new(Slab::new()), |slab| {
                copse::BTreeSet::new(SlabOrdering(Some(unsafe { &*slab.get() })))
            }),
        }
    }
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            inner: SortedSlabInner::new(UnsafeCell::new(Slab::with_capacity(cap)), |slab| {
                copse::BTreeSet::new(SlabOrdering(Some(unsafe { &*slab.get() })))
            }),
        }
    }
    pub fn iter(&self) -> Iter<T> {
        self.inner
            .with(|f| Iter(f.sorted.iter(), unsafe { &*f.backing.get() }))
    }
    /// Get an element by its index.
    pub fn get(&self, index: usize) -> Option<&T> {
        unsafe { (&*self.inner.borrow_backing().get()).get(index) }
    }
}

impl<T: Ord + 'static> SortedSlab<T> {
    /// Modify an element by its index. Because we need to update the b-tree after the mutation,
    /// this can only take a closure.
    pub fn with_mut<R, F: FnOnce(&mut T) -> R>(&mut self, index: usize, f: F) -> Option<R> {
        if self.inner.with_sorted_mut(|s| s.remove(&index)) {
            let ret = self
                .inner
                .with_backing(|b| f(&mut unsafe { &mut *b.get() }[index]));
            self.inner.with_sorted_mut(|s| s.insert(index));
            Some(ret)
        } else {
            None
        }
    }
}

impl<T: Debug + 'static> Debug for SortedSlab<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.debug_set().entries(self.iter()).finish()
    }
}

impl<T: PartialEq + 'static> PartialEq for SortedSlab<T> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().eq(other.iter())
    }
}
impl<T: Eq> Eq for SortedSlab<T> {}
impl<T: PartialOrd + 'static> PartialOrd for SortedSlab<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.iter().partial_cmp(other.iter())
    }
}
impl<T: Ord + 'static> Ord for SortedSlab<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.iter().cmp(other.iter())
    }
}

impl<T: 'static> IntoIterator for SortedSlab<T> {
    type IntoIter = IntoIter<T>;
    type Item = T;
    fn into_iter(mut self) -> IntoIter<T> {
        let iter = self.inner.with_sorted_mut(|s| {
            std::mem::replace(s, copse::BTreeSet::new(SlabOrdering(None))).into_iter()
        });
        let slab = self.inner.into_heads().backing.into_inner();
        IntoIter(iter, slab)
    }
}

impl<'a, T: 'static> IntoIterator for &'a SortedSlab<T> {
    type IntoIter = Iter<'a, T>;
    type Item = &'a T;
    fn into_iter(self) -> Iter<'a, T> {
        self.iter()
    }
}

#[derive(Debug, Default, Clone, Copy)]
struct SlabOrdering<'a, T>(Option<&'a Slab<T>>);
impl<'a, T: Ord> copse::TotalOrder for SlabOrdering<'a, T> {
    type OrderedType = usize;
    fn cmp(&self, lhs: &usize, rhs: &usize) -> Ordering {
        let slab = self.0.unwrap();
        slab[*lhs].cmp(&slab[*rhs])
    }
}
impl<T: Ord> copse::SortableBy<SlabOrdering<'_, T>> for usize {
    fn sort_key(&self) -> &usize {
        self
    }
}

pub struct Iter<'a, T>(copse::btree_set::Iter<'a, usize>, &'a Slab<T>);
impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;
    fn next(&mut self) -> Option<&'a T> {
        self.0.next().map(|i| &self.1[*i])
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}
impl<T> FusedIterator for Iter<'_, T> {}
impl<T> ExactSizeIterator for Iter<'_, T> {
    fn len(&self) -> usize {
        self.0.len()
    }
}
pub struct IntoIter<T>(copse::btree_set::IntoIter<usize>, Slab<T>);
impl<T> Iterator for IntoIter<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        self.0.next().map(|i| self.1.remove(i))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}
impl<T> FusedIterator for IntoIter<T> {}
impl<T> ExactSizeIterator for IntoIter<T> {
    fn len(&self) -> usize {
        self.0.len()
    }
}
