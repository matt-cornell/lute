use crate::utils::bitset::BitSet;

const BITS: usize = std::mem::size_of::<usize>() * 8;

fn make_set(bits: impl IntoIterator<Item = usize>) -> BitSet<usize, 2> {
    let mut out = BitSet::<usize, 2>::new();
    for bit in bits {
        out.set(bit, true);
    }
    out
}

#[cfg(feature = "rand")]
#[test]
fn bits() {
    use rand::prelude::*;
    let mut rng = rand::thread_rng();
    let mut indices = Vec::<u16>::new();
    let mut set = BitSet::<usize, 2>::new();
    for _ in 0..100 {
        let count = rng.gen_range(100..1000);
        indices.resize(count, 0);
        rng.fill(&mut indices[..]);
        indices.sort();
        indices.dedup();

        set.clear();
        for &i in &indices {
            set.set(i as _, true);
        }

        assert_eq!(set.count_ones(), indices.len());

        for i in 0..(indices.as_slice().len() * BITS) {
            let contained = indices.binary_search(&(i as _)).is_ok();
            assert_eq!(set.get(i), contained);
        }
    }
}

#[test]
fn nth() {
    {
        let set = make_set([]);
        assert_eq!(set.nth(0), None);
    }
    {
        let mut set = make_set([65]);
        set.set(65, false);
        assert_eq!(set.nth(0), None);
    }
    {
        let set = make_set([0, 1, 2]);
        assert_eq!(set.nth(0), Some(0));
        assert_eq!(set.nth(1), Some(1));
        assert_eq!(set.nth(2), Some(2));
        assert_eq!(set.nth(3), None);
    }
    {
        let set = make_set([1, 2, 3]);
        assert_eq!(set.nth(0), Some(1));
        assert_eq!(set.nth(1), Some(2));
        assert_eq!(set.nth(2), Some(3));
        assert_eq!(set.nth(3), None);
    }
    {
        let set = make_set([0, 64, 96, 128]);
        assert_eq!(set.nth(0), Some(0));
        assert_eq!(set.nth(1), Some(64));
        assert_eq!(set.nth(2), Some(96));
        assert_eq!(set.nth(3), Some(128));
        assert_eq!(set.nth(4), None);
    }
}
