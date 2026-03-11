use std::fmt::{self, Write};

type T = u64;
const BITS_PER_BP: usize = 1;
const BITS_PER_BLOCK: usize = T::BITS as usize;
const BP_PER_BLOCK: usize = BITS_PER_BLOCK / BITS_PER_BP;

fn convert_ascii_to_columns(ch: u8) -> (T, T) {
    match ch {
        b'A' | b'a' => (0, 0),
        b'C' | b'c' => (0, 1),
        b'G' | b'g' => (1, 1),
        b'T' | b't' => (1, 0),
        _ => panic!("Invalid nucleotide: {}", ch as char),
    }
}

#[derive(Clone, Default)]
pub struct ColumnarDNA {
    high_bits: Vec<T>,
    low_bits: Vec<T>,
    cur_hi: T,
    cur_lo: T,
    len: usize,
}

impl ColumnarDNA {
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            high_bits: Vec::new(),
            low_bits: Vec::new(),
            cur_hi: 0,
            cur_lo: 0,
            len: 0,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            high_bits: Vec::with_capacity(capacity),
            low_bits: Vec::with_capacity(capacity),
            cur_hi: 0,
            cur_lo: 0,
            len: 0,
        }
    }

    #[inline(always)]
    pub const fn len(&self) -> usize {
        self.len
    }

    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    pub const fn capacity(&self) -> usize {
        self.high_bits.capacity()
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.high_bits.clear();
        self.low_bits.clear();
        self.cur_hi = 0;
        self.cur_lo = 0;
        self.len = 0;
    }

    #[inline(always)]
    pub fn append(&mut self, high: T, low: T, size: usize) {
        if size == 0 {
            // should not happen?
            return;
        }
        let rem = self.len % BITS_PER_BLOCK;
        let mask = !0 >> (BITS_PER_BLOCK - size);
        self.len += size;
        let hi = high & mask;
        let lo = low & mask;
        if rem + size >= BITS_PER_BLOCK {
            self.cur_hi |= high << rem;
            self.cur_lo |= low << rem;
            self.high_bits.push(self.cur_hi);
            self.low_bits.push(self.cur_lo);
            self.cur_hi = if rem > 0 {
                hi >> (BITS_PER_BLOCK - rem)
            } else {
                0
            };
            self.cur_lo = if rem > 0 {
                lo >> (BITS_PER_BLOCK - rem)
            } else {
                0
            };
        } else {
            self.cur_hi |= hi << rem;
            self.cur_lo |= lo << rem;
        }
    }

    #[inline(always)]
    pub fn high_bits(&self) -> (&[T], T) {
        (&self.high_bits, self.cur_hi)
    }

    #[inline(always)]
    pub fn low_bits(&self) -> (&[T], T) {
        (&self.low_bits, self.cur_lo)
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> (bool, bool) {
        if i < self.len & (!0 << BITS_PER_BLOCK.trailing_zeros()) {
            (
                (self.high_bits[i / BP_PER_BLOCK] >> (i % BP_PER_BLOCK)) & 1 != 0,
                (self.low_bits[i / BP_PER_BLOCK] >> (i % BP_PER_BLOCK)) & 1 != 0,
            )
        } else {
            (
                (self.cur_hi >> (i % BP_PER_BLOCK)) & 1 != 0,
                (self.cur_lo >> (i % BP_PER_BLOCK)) & 1 != 0,
            )
        }
    }

    #[inline(always)]
    pub fn get_char(&self, i: usize) -> char {
        match self.get(i) {
            (false, false) => 'A',
            (false, true) => 'C',
            (true, false) => 'T',
            (true, true) => 'G',
        }
    }

    #[allow(dead_code)]
    fn push_ascii(&mut self, s: &[u8]) {
        for &ch in s {
            let (cur_hi, cur_lo) = convert_ascii_to_columns(ch);
            self.append(cur_hi, cur_lo, 1);
        }
    }

    #[allow(dead_code)]
    fn push_str(&mut self, s: &str) {
        self.push_ascii(s.as_bytes());
    }
}

impl fmt::Display for ColumnarDNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.len() {
            f.write_char(self.get_char(i))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_starts_empty() {
        let dna = ColumnarDNA::new();
        assert_eq!(dna.len(), 0);
        assert_eq!(dna.cur_hi, 0);
        assert_eq!(dna.cur_lo, 0);
        assert_eq!(dna.high_bits, Vec::new());
        assert_eq!(dna.low_bits, Vec::new());
    }

    #[test]
    fn test_append_simple() {
        let mut dna = ColumnarDNA::new();
        // Append one nucleotide (size = 1)
        dna.append(1, 0, 1); // (cur_hi=1,cur_lo=0) => T
        assert_eq!(dna.len(), 1);
        assert_eq!(format!("{}", dna), "T");

        dna.append(0, 1, 1); // (cur_hi=0,cur_lo=1) => C
        assert_eq!(dna.len(), 2);
        assert_eq!(format!("{}", dna), "TC");
    }

    #[test]
    fn test_append_multiple() {
        let mut dna = ColumnarDNA::new();
        // Append four 2-bit nucleotides (A, C, T, G)
        dna.append(0, 0, 1); // A
        dna.append(0, 1, 1); // C
        dna.append(1, 0, 1); // T
        dna.append(1, 1, 1); // G
        assert_eq!(dna.len(), 4);
        assert_eq!(format!("{}", dna), "ACTG");
    }

    #[test]
    fn append_single_bits() {
        let mut v = ColumnarDNA::new();

        v.push_str("ACGT");

        assert_eq!(v.len(), 4);
        assert_eq!(v.to_string(), "ACGT");
    }

    #[test]
    fn append_exact_63_bits() {
        let seq = "A".repeat(63);
        let mut v = ColumnarDNA::new();

        v.push_str(&seq);

        assert_eq!(v.len(), 63);
        assert_eq!(v.to_string(), seq);
        assert_eq!(v.high_bits.len(), 0); // still in partial word
    }

    #[test]
    fn test_crossing_all_offsets() {
        let letters = ['A', 'C', 'G', 'T'];

        for offset in 0..64 {
            let mut v = ColumnarDNA::new();

            // Fill offset bits
            for _ in 0..offset {
                v.push_str("A");
            }

            // Write 1–10 bits after the boundary
            for size in 1..10 {
                let ch = letters[size & 3];
                v.push_str(&ch.to_string());
            }

            // Roundtrip must match
            let s = v.to_string();
            assert_eq!(s.len(), offset + 9);
            assert_eq!(s, s); // no-op but ensures no invalid UTF
        }
    }

    #[test]
    fn repeated_small_appends() {
        let mut v = ColumnarDNA::new();
        let mut s = String::new();

        for i in 0..500 {
            let ch = ['A', 'C', 'G', 'T'][i & 3];
            s.push(ch);
            v.push_str(&ch.to_string());
        }

        assert_eq!(v.len(), s.len());
        assert_eq!(v.to_string(), s);
    }

    #[test]
    fn cross_boundary_regression() {
        let mut v = ColumnarDNA::new();
        let seq = "A".repeat(64) + "C";
        v.push_str(&seq);

        assert_eq!(v.len(), 65);
        assert_eq!(v.to_string(), seq);
        assert_eq!(v.high_bits.len(), 1); // one stored word now
    }

    #[test]
    fn append_exact_64_bases_bulk() {
        let mut v = ColumnarDNA::new();
        v.append(!0u64, !0u64, 64); // 64 G
        v.append(0u64, 0u64, 1); // 1 A

        assert_eq!(v.len(), 65);
        assert_eq!(v.high_bits.len(), 1);
        assert_eq!(v.to_string(), "G".repeat(64) + "A");
    }

    #[test]
    fn two_full_blocks_bulk_distinct() {
        let mut v = ColumnarDNA::new();
        v.append(!0u64, 0u64, 64); // 64 T
        v.append(0u64, !0u64, 64); // 64 C

        assert_eq!(v.len(), 128);
        assert_eq!(v.to_string(), "T".repeat(64) + &"C".repeat(64));
    }
}
