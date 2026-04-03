use crate::kmer::{BitStorage, Kmer, Mer};
use crate::kmer_collection::MerChunk;
use bitvec::prelude::BitVec;

#[derive(Clone, Debug)]
pub struct NullableMerChunk<const K: usize, T: BitStorage> {
    pub merchunk: MerChunk<K, T>,
    pub occupied: BitVec,
}

impl<const K: usize, T: BitStorage> NullableMerChunk<K, T> {
    /// Creates a new chunk of given size, with all entries inoccupied.
    /// Underlying storage is initialized with a sentinel value.
    pub fn new(size: usize) -> Self {
        let sentinel = Mer::<K, T>::some_kmer();

        let v0 = vec![sentinel.0; size];
        let v1 = vec![sentinel.1; size];
        let occupied = BitVec::repeat(false, size);

        Self {
            merchunk: MerChunk(v0, v1),
            occupied,
        }
    }

    /// Returns the number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.occupied.len()
    }

    /// Returns true if the chunk is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns whether the entry at `offset` is occupied.
    #[inline]
    pub fn is_occupied(&self, offset: usize) -> bool {
        self.occupied[offset]
    }

    /// Safe access: returns `Some` if occupied, otherwise `None`.
    #[inline]
    pub fn get(&self, offset: usize) -> Option<Mer<K, T>> {
        if self.occupied[offset] {
            Some(self.merchunk.get(offset))
        } else {
            None
        }
    }

    /// Returns the value without checking occupiedity.
    /// Caller must ensure `is_occupied(offset)` holds.
    #[inline]
    pub fn get_assume_occupied(&self, offset: usize) -> Mer<K, T> {
        debug_assert!(self.occupied[offset]);
        self.merchunk.get(offset)
    }

    /// Sets a value and marks the slot as occupied.
    #[inline]
    pub fn set(&mut self, offset: usize, mer: Mer<K, T>) {
        self.merchunk.set(offset, mer);
        self.occupied.set(offset, true);
    }

    /// Clears a value (marks as inoccupied).
    /// Underlying data is left unchanged.
    #[inline]
    pub fn clear(&mut self, offset: usize) {
        self.occupied.set(offset, false);
    }

    /// Clears all entries.
    #[inline]
    pub fn clear_all(&mut self) {
        self.occupied.fill(false);
    }

    /// Returns a reference to the underlying chunk.
    /// Restricted to read-only to preserve invariants.
    #[inline]
    pub fn as_chunk(&self) -> &MerChunk<K, T> {
        &self.merchunk
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nullable_new() {
        let size = 5;
        let chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(size);

        assert_eq!(chunk.len(), size);
        for i in 0..size {
            assert!(!chunk.is_occupied(i));
            assert_eq!(chunk.get(i), None);
        }
    }

    #[test]
    fn test_nullable_set_and_get() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(3);

        let mer = Mer::<4, u8>::some_kmer();
        chunk.set(1, mer);

        assert!(!chunk.is_occupied(0));
        assert!(chunk.is_occupied(1));
        assert!(!chunk.is_occupied(2));

        assert_eq!(chunk.get(0), None);
        assert_eq!(chunk.get(1), Some(mer));
        assert_eq!(chunk.get(2), None);
    }

    #[test]
    fn test_nullable_overwrite() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(1);

        let mer1 = Mer::<4, u8>::some_kmer();
        let mer2 = Mer::<4, u8>::some_kmer();

        chunk.set(0, mer1);
        assert_eq!(chunk.get(0), Some(mer1));

        chunk.set(0, mer2);
        assert_eq!(chunk.get(0), Some(mer2));
    }

    #[test]
    fn test_nullable_clear() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(2);

        let mer = Mer::<4, u8>::some_kmer();
        chunk.set(0, mer);

        assert!(chunk.is_occupied(0));
        assert_eq!(chunk.get(0), Some(mer));

        chunk.clear(0);

        assert!(!chunk.is_occupied(0));
        assert_eq!(chunk.get(0), None);
    }

    #[test]
    fn test_nullable_clear_all() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(3);

        let mer = Mer::<4, u8>::some_kmer();
        for i in 0..3 {
            chunk.set(i, mer);
        }

        chunk.clear_all();

        for i in 0..3 {
            assert!(!chunk.is_occupied(i));
            assert_eq!(chunk.get(i), None);
        }
    }

    #[test]
    fn test_nullable_sparse_usage() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(5);

        let mer = Mer::<4, u8>::some_kmer();

        chunk.set(1, mer);
        chunk.set(3, mer);

        for i in 0..5 {
            if i == 1 || i == 3 {
                assert!(chunk.is_occupied(i));
                assert_eq!(chunk.get(i), Some(mer));
            } else {
                assert!(!chunk.is_occupied(i));
                assert_eq!(chunk.get(i), None);
            }
        }
    }

    #[test]
    fn test_nullable_get_assume_occupied() {
        let mut chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(1);

        let mer = Mer::<4, u8>::some_kmer();
        chunk.set(0, mer);

        let retrieved = chunk.get_assume_occupied(0);
        assert_eq!(retrieved, mer);
    }

    #[test]
    #[should_panic]
    fn test_nullable_get_assume_occupied_panics_on_inoccupied_in_debug() {
        let chunk: NullableMerChunk<4, u8> = NullableMerChunk::new(1);

        // Should trigger debug_assert in debug builds
        let _ = chunk.get_assume_occupied(0);
    }
}
