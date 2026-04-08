use crate::kmer::{BitStorage, Kmer, Mer};
use crate::kmer_collection::MerChunk;
use std::fmt::Debug;

use bitvec::prelude::BitVec;

use std::cmp::max;

const LEFT_OVERLOW_SAFETY: usize = 1000;

pub trait Monoid: Clone + Copy + Debug {
    fn identity() -> Self;
    fn combine(&self, other: &Self) -> Self;
    #[inline(always)]
    fn accumutate(&mut self, other: &Self) {
        *self = self.combine(other);
    }
}

#[derive(Clone, Copy, Debug)]
struct Data<const K: usize, T: BitStorage, V: Monoid> {
    pub hash: u64,
    pub mer: Mer<K, T>,
    pub value: V,
}

impl<const K: usize, T: BitStorage, V: Monoid> Data<K, T, V> {
    fn some() -> Self {
        Self {
            hash: 0,
            mer: Mer::<K, T>::any(),
            value: V::identity(),
        }
    }
}

pub struct MerGroupedMap<const K: usize, T: BitStorage, V: Monoid> {
    data: Vec<Data<K, T, V>>,
    occupied: BitVec,
    size: usize,

    // capacity = 2^log2_capacity
    pub log2_capacity: usize,

    // number of high bits already consumed by upstream clustering
    hash_bit_offset: usize,

    // instrumentation
    collision_count: usize,
    max_run: usize,
}

impl<const K: usize, T: BitStorage, V: Monoid> MerGroupedMap<K, T, V> {
    pub fn new(log2_capacity: usize, hash_bit_offset: usize) -> Self {
        assert!(
            hash_bit_offset + log2_capacity < 64,
            "hash_bit_offset:{hash_bit_offset}, log2_capacity:{log2_capacity}"
        );
        let capacity = (1usize << log2_capacity) + LEFT_OVERLOW_SAFETY;

        Self {
            data: vec![Data::<K, T, V>::some(); capacity],
            occupied: BitVec::repeat(false, capacity),
            size: 0,
            log2_capacity,
            hash_bit_offset,
            collision_count: 0,
            max_run: 0,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.size
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
    #[inline(always)]
    fn capacity(&self) -> usize {
        self.data.capacity()
    }

    /// Extract k bits after skipping the already-consumed MSBs
    #[inline(always)]
    fn extract_from_hash(&self, hash: u64) -> usize {
        let shift = 64 - (self.hash_bit_offset + self.log2_capacity);
        ((hash >> shift) as usize) & ((1 << self.log2_capacity) - 1)
    }

    /// Insert with:
    /// - linear probing (no wraparound)
    /// - ordered-by-hash swapping
    /// - monoid aggregation on equal keys
    #[inline(always)]
    pub fn add(&mut self, mer: Mer<K, T>, v: V) {
        self.add_hash(mer.hash(), mer, v);
    }

    pub fn add_hash(&mut self, mut hash: u64, mer: Mer<K, T>, mut v: V) {
        let mut mer_to_insert = mer;

        let mut offset = self.extract_from_hash(hash);
        let capacity = self.capacity();
        let mut run = 0;

        while offset < capacity {
            if !self.occupied[offset] {
                self.data[offset] = Data {
                    hash,
                    mer: mer_to_insert,
                    value: v,
                };
                self.occupied.set(offset, true);
                self.size += 1;
                self.max_run = max(run, self.max_run);
                return;
            }

            let existing_data = &mut self.data[offset];

            // --- aggregation case ---
            if existing_data.mer == mer_to_insert {
                existing_data.value.accumutate(&v);
                self.max_run = max(run, self.max_run);
                return;
            }

            self.collision_count += 1;
            run += 1;

            // --- ordered insertion (insertion-sort behavior) ---
            if existing_data.hash > hash {
                let tmp = std::mem::replace(
                    existing_data,
                    Data {
                        hash,
                        mer: mer_to_insert,
                        value: v,
                    },
                );
                mer_to_insert = tmp.mer;
                hash = tmp.hash;
                v = tmp.value;
            }
            offset += 1;
        }

        // No wraparound, hard failure if suffix is full
        panic!("MerGroupedMap: insertion overflow (increase capacity or reduce load factor)");
    }

    /// Compact into dense structures.
    /// Output is already sorted by full hash.
    pub fn compactify(&self) -> (Vec<u64>, MerChunk<K, T>, Vec<V>) {
        let mut out_hashes = Vec::with_capacity(self.size);
        let mut out_mers = MerChunk::<K, T>::with_capacity(self.size);
        let mut out_vals = Vec::with_capacity(self.size);

        for i in self.occupied.iter_ones() {
            let data = self.data[i];
            out_hashes.push(data.hash);
            out_mers.push(data.mer);
            out_vals.push(data.value);
        }

        (out_hashes, out_mers, out_vals)
    }
    pub fn collision_count(&self) -> usize {
        self.collision_count
    }

    pub fn max_run(&self) -> usize {
        self.max_run
    }
}

impl Monoid for u16 {
    #[inline(always)]
    fn identity() -> Self {
        0
    }

    #[inline(always)]
    fn combine(&self, other: &Self) -> Self {
        self + other
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NullableMerChunk;
    use crate::kmer::Mer;

    /// Simple additive monoid for testing
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct Count(u32);

    impl Monoid for Count {
        fn identity() -> Self {
            Count(0)
        }
        fn combine(&self, other: &Self) -> Self {
            Count(self.0 + other.0)
        }
    }
    #[test]
    fn u16_identity() {
        let x: u16 = 42;
        assert_eq!(u16::identity().combine(&x), x);
        assert_eq!(x.combine(&u16::identity()), x);
    }

    #[test]
    fn u16_associativity() {
        let a: u16 = 1;
        let b: u16 = 2;
        let c: u16 = 3;

        let left = a.combine(&b).combine(&c);
        let right = a.combine(&b.combine(&c));

        assert_eq!(left, right);
    }

    #[test]
    fn u16_accumulate() {
        let mut x: u16 = 5;
        x.accumutate(&3);
        assert_eq!(x, 8);

        x.accumutate(&u16::identity());
        assert_eq!(x, 8);
    }

    #[test]
    fn count_identity() {
        let x = Count(10);
        assert_eq!(Count::identity().combine(&x), x);
        assert_eq!(x.combine(&Count::identity()), x);
    }

    #[test]
    fn count_associativity() {
        let a = Count(1);
        let b = Count(2);
        let c = Count(3);

        let left = a.combine(&b).combine(&c);
        let right = a.combine(&b.combine(&c));

        assert_eq!(left, right);
    }

    #[test]
    fn count_accumulate() {
        let mut x = Count(7);
        x.accumutate(&Count(5));
        assert_eq!(x, Count(12));

        x.accumutate(&Count::identity());
        assert_eq!(x, Count(12));
    }

    #[test]
    fn mixed_sequence_count() {
        let values = [Count(1), Count(2), Count(3), Count(4)];

        let mut acc = Count::identity();
        for v in &values {
            acc.accumutate(v);
        }

        assert_eq!(acc, Count(10));
    }

    // Helper: create a Mer<K, u8> from a number (for testing)
    fn mer_from_u8(val: u8) -> Mer<4, u8> {
        Mer::<4, u8>(val.into(), val.into())
    }

    fn mer_from_u32(val: u32) -> Mer<31, u32> {
        Mer::<31, u32>(val.into(), val.into())
    }

    #[test]
    fn test_nullable_merchunk_basic() {
        let mut nmc: NullableMerChunk<4, u8> = NullableMerChunk::new(5);

        assert_eq!(nmc.occupied.len(), 5);
        for i in 0..5 {
            assert!(!nmc.is_occupied(i));
            assert!(nmc.get(i).is_none());
        }

        let m = mer_from_u8(42);
        nmc.set(2, m);

        assert!(nmc.is_occupied(2));
        assert_eq!(nmc.get(2), Some(m));
        assert!(nmc.get(1).is_none());
    }

    #[test]
    fn test_mer_grouped_map_insert_and_len() {
        let mut map: MerGroupedMap<4, u8, Count> = MerGroupedMap::new(4, 0);

        let m1 = mer_from_u8(10);
        let m2 = mer_from_u8(20);
        let m3 = mer_from_u8(15);

        map.add(m1, Count(1));
        map.add(m2, Count(2));
        map.add(m3, Count(3));

        assert_eq!(map.len(), 3);
    }

    #[test]
    fn test_mer_grouped_map_duplicate_aggregation() {
        let mut map: MerGroupedMap<4, u8, Count> = MerGroupedMap::new(4, 0);

        let m = mer_from_u8(99);

        map.add(m, Count(1));
        map.add(m, Count(2));

        assert_eq!(map.len(), 1);

        let (_, mers, values) = map.compactify();
        assert_eq!(mers.len(), 1);
        assert_eq!(values[0], Count(3));
    }

    #[test]
    fn test_mer_grouped_map_ordering_after_insert() {
        let mut map: MerGroupedMap<4, u8, Count> = MerGroupedMap::new(4, 0);

        let vals: Vec<u8> = vec![42, 7, 99, 15, 50];
        for &v in &vals {
            map.add(mer_from_u8(v), Count(1));
        }

        let (_, mers, _) = map.compactify();
        let mut prev_hash = 0;
        for mer in mers.iter() {
            let h = mer.hash();
            assert!(h >= prev_hash);
            prev_hash = h;
        }
    }

    #[test]
    #[should_panic]
    fn test_mer_grouped_map_overflow_panics() {
        let mut map: MerGroupedMap<_, _, Count> = MerGroupedMap::new(2, 0);
        // capacity = 4
        for i in 0..1001 {
            map.add(mer_from_u32(i), Count(1));
        }
    }
}
