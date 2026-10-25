use alloc::vec::Vec;
use core::fmt::{self, Write};

#[cfg(feature = "packed-seq")]
use packed_seq::{PackedSeq, PackedSeqVec, SeqVec};

type T = u64;
const BITS_PER_BP: usize = 2;
const BITS_PER_BLOCK: usize = T::BITS as usize;
const BP_PER_BLOCK: usize = BITS_PER_BLOCK / BITS_PER_BP;
const WORD_BYTES: usize = size_of::<T>();
// `cur` followed by the 48 bytes of padding of `PackedSeq`
const MARGIN: usize = if cfg!(feature = "packed-seq") {
    1 + 6
} else {
    0
};

/// DNA sequence packed on 2 bits per base: A = 0, C = 1, T = 2, G = 3.
#[derive(Default)]
pub struct PackedDNA {
    /// The main storage.
    /// Never shrinks.
    bits: Vec<u8>,
    /// Trailing characters that have not yet been written to `bits`.
    cur: T,
    num_bits: usize,
    init_bytes: usize,
}

impl Clone for PackedDNA {
    fn clone(&self) -> Self {
        let mut res = Self::with_capacity(self.len());
        res.bits.extend_from_slice(&self.bits);
        res.cur = self.cur;
        res.num_bits = self.num_bits;
        res.init_bytes = res.bits.len();
        res
    }
}

impl PackedDNA {
    /// Creates an empty sequence without allocating.
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            bits: Vec::new(),
            cur: 0,
            num_bits: 0,
            init_bytes: 0,
        }
    }

    /// Creates an empty sequence with room for at least `capacity` bases.
    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: Vec::with_capacity((capacity.div_ceil(BP_PER_BLOCK) + MARGIN) * WORD_BYTES),
            cur: 0,
            num_bits: 0,
            init_bytes: 0,
        }
    }

    /// Returns the number of bases.
    #[inline(always)]
    pub const fn len(&self) -> usize {
        self.num_bits / BITS_PER_BP
    }

    /// Returns `true` if there are no bases.
    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.num_bits == 0
    }

    /// Returns the number of bases that fit without reallocating.
    #[inline(always)]
    pub const fn capacity(&self) -> usize {
        (self.bits.capacity() / WORD_BYTES).saturating_sub(MARGIN) * BP_PER_BLOCK
    }

    #[inline(always)]
    const fn num_words(&self) -> usize {
        self.bits.len() / WORD_BYTES
    }

    #[inline(always)]
    const unsafe fn word_ptr(&mut self, word: usize) -> *mut T {
        unsafe { self.bits.as_mut_ptr().add(word * WORD_BYTES) as *mut T }
    }

    #[inline(always)]
    unsafe fn set_num_words(&mut self, num_words: usize) {
        unsafe { self.bits.set_len(num_words * WORD_BYTES) }
    }

    #[inline(always)]
    fn ensure_spare(&mut self, additional_words: usize) {
        if self.bits.len() + (additional_words + MARGIN) * WORD_BYTES > self.bits.capacity() {
            self.grow(additional_words);
        }
    }

    #[cold]
    #[inline(never)]
    fn grow(&mut self, additional_words: usize) {
        self.bits.reserve((additional_words + MARGIN) * WORD_BYTES);
        // A reallocation only preserves the bytes in `0..len`.
        self.init_bytes = self.bits.len();
    }

    /// Reserves room for at least `additional_bases` more bases.
    #[inline(always)]
    pub fn reserve(&mut self, additional_bases: usize) {
        self.ensure_spare(additional_bases / BP_PER_BLOCK + 2);
    }

    /// Removes all bases, keeping the allocated memory.
    #[inline(always)]
    pub fn clear(&mut self) {
        self.bits.clear();
        self.cur = 0;
        self.num_bits = 0;
    }

    #[inline(always)]
    pub(crate) fn pop_base(&mut self) {
        if self.num_bits >= BITS_PER_BP {
            let rem = self.num_bits % BITS_PER_BLOCK;
            self.num_bits -= BITS_PER_BP;
            if rem == 0 {
                let last_word = self.num_words() - 1;
                unsafe {
                    let last = self.word_ptr(last_word).read_unaligned();
                    self.set_num_words(last_word);
                    self.cur = last & ((1u64 << (BITS_PER_BLOCK - BITS_PER_BP)) - 1);
                }
            } else if rem >= BITS_PER_BP {
                self.cur &= (1u64 << (rem - BITS_PER_BP)) - 1;
            }
        }
    }

    /// Appends 64 bases packed on 2 bits each.
    #[inline(always)]
    pub fn append_full(&mut self, packed: u128) {
        self.ensure_spare(2);
        let lo = packed as u64;
        let hi = (packed >> 64) as u64;
        let rem = self.num_bits % BITS_PER_BLOCK;
        self.num_bits += 128;
        let len = self.num_words();
        unsafe {
            let ptr = self.word_ptr(len);
            if rem == 0 {
                ptr.write_unaligned(lo);
                ptr.add(1).write_unaligned(hi);
            } else {
                let inv = BITS_PER_BLOCK - rem;
                ptr.write_unaligned(self.cur | (lo << rem));
                ptr.add(1).write_unaligned((lo >> inv) | (hi << rem));
                self.cur = hi >> inv;
            }
            self.set_num_words(len + 2);
        }
    }

    #[inline(always)]
    pub(crate) fn pack_blocks<const CONFIG: crate::config::Config>(
        &mut self,
        data: &[u8],
        q: usize,
    ) {
        self.ensure_spare(2 * q);
        let rem = self.num_bits % BITS_PER_BLOCK;
        self.num_bits += 128 * q;
        let start_len = self.num_words();
        unsafe {
            let mut p = self.word_ptr(start_len);
            if rem == 0 {
                for k in 0..q {
                    let m = crate::simd::extract_dna_bitmask::<CONFIG>(&data[64 * k..64 * (k + 1)]);
                    p.write_unaligned(m.two_bits as u64);
                    p.add(1).write_unaligned((m.two_bits >> 64) as u64);
                    p = p.add(2);
                }
            } else {
                let inv = BITS_PER_BLOCK - rem;
                let mut cur = self.cur;
                for k in 0..q {
                    let m = crate::simd::extract_dna_bitmask::<CONFIG>(&data[64 * k..64 * (k + 1)]);
                    let lo = m.two_bits as u64;
                    let hi = (m.two_bits >> 64) as u64;
                    p.write_unaligned(cur | (lo << rem));
                    p.add(1).write_unaligned((lo >> inv) | (hi << rem));
                    cur = hi >> inv;
                    p = p.add(2);
                }
                self.cur = cur;
            }
            self.set_num_words(start_len + 2 * q);
        }
    }

    /// Appends the `num_bases` bases packed in the low bits of `packed`, other bits are ignored.
    /// `num_bases` must be at most 64.
    #[inline(always)]
    pub fn append(&mut self, packed: u128, num_bases: usize) {
        debug_assert!(num_bases <= 2 * BP_PER_BLOCK);
        let num_bits = BITS_PER_BP * num_bases;
        if num_bits == 0 {
            return;
        }
        let lo = packed as u64;
        let hi = (packed >> 64) as u64;
        if num_bits <= 64 {
            self.append_u64(lo, num_bits);
        } else {
            self.append_u64(lo, 64);
            self.append_u64(hi, num_bits - 64);
        }
    }

    #[inline(always)]
    fn append_u64(&mut self, val: u64, num_bits: usize) {
        let mask = if num_bits == 64 {
            !0
        } else {
            (1u64 << num_bits) - 1
        };
        let x = val & mask;
        let rem = self.num_bits % BITS_PER_BLOCK;
        self.num_bits += num_bits;
        if rem + num_bits >= BITS_PER_BLOCK {
            self.ensure_spare(1);
            let len = self.num_words();
            unsafe {
                self.word_ptr(len).write_unaligned(self.cur | (x << rem));
                self.set_num_words(len + 1);
            }
            self.cur = if rem > 0 && num_bits > (BITS_PER_BLOCK - rem) {
                x >> (BITS_PER_BLOCK - rem)
            } else {
                0
            };
        } else {
            self.cur |= x << rem;
        }
    }

    /// The slice of characters, and separately the incomplete last word.
    ///
    /// Use in combination with `.len()`.
    #[inline(always)]
    pub fn bits(&self) -> (&[T], T) {
        let (prefix, words, _) = unsafe { self.bits.align_to::<T>() };
        assert!(prefix.is_empty(), "unaligned PackedDNA buffer");
        (words, self.cur)
    }

    /// The length in characters, and the raw data.
    ///
    /// The final byte might be incomplete.
    #[inline(always)]
    pub fn into_packed_bytes(mut self) -> (usize, Vec<u8>) {
        let len = self.len();
        let word = self.num_words();
        self.bits.resize((word + 1) * WORD_BYTES, 0);
        unsafe { self.word_ptr(word).write_unaligned(self.cur) };
        (len, self.bits)
    }

    /// Returns the 2-bit code of base `i`, which must be lower than [`len`](Self::len).
    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        let block_idx = i / BP_PER_BLOCK;
        let bit_offset = 2 * (i % BP_PER_BLOCK);
        let word = if block_idx < self.num_words() {
            let bytes = &self.bits[block_idx * WORD_BYTES..(block_idx + 1) * WORD_BYTES];
            T::from_ne_bytes(bytes.try_into().unwrap())
        } else {
            self.cur
        };
        ((word >> bit_offset) & 0b11) as u8
    }

    /// Returns base `i` as `A`, `C`, `T` or `G`.
    #[inline(always)]
    pub fn get_char(&self, i: usize) -> char {
        const LUT: [char; 4] = ['A', 'C', 'T', 'G'];
        LUT[self.get(i) as usize]
    }

    /// Write `cur` to the end of `bits`, and add padding bytes for `packed-seq` if the feature is enabled.
    ///
    /// Required before [`as_packed_seq`](Self::as_packed_seq). Parsers call it after each event.
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    pub fn pad(&mut self) {
        self.ensure_spare(0);
        if self.init_bytes < self.bits.len() + MARGIN * WORD_BYTES {
            self.init_padding();
        }
        let len = self.num_words();
        unsafe { self.word_ptr(len).write_unaligned(self.cur) };
    }

    #[cfg(feature = "packed-seq")]
    #[cold]
    #[inline(never)]
    fn init_padding(&mut self) {
        let start = self.init_bytes.max(self.bits.len());
        let end = self.bits.len() + MARGIN * WORD_BYTES;
        unsafe { core::ptr::write_bytes(self.bits.as_mut_ptr().add(start), 0, end - start) };
        self.init_bytes = end;
    }

    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    fn padded_bytes(&self) -> &[u8] {
        static EMPTY: [u8; MARGIN * WORD_BYTES] = [0; MARGIN * WORD_BYTES];
        if self.bits.capacity() == 0 && self.num_bits == 0 {
            return &EMPTY;
        }
        let len = self.bits.len();
        let padded_len = len + MARGIN * WORD_BYTES;
        debug_assert!(
            self.init_bytes >= padded_len,
            "PackedDNA::pad must be called first"
        );
        // Safety: `init_bytes` is reset on reallocation, so these bytes are initialized.
        let padded = unsafe { core::slice::from_raw_parts(self.bits.as_ptr(), padded_len) };
        let stored_cur = T::from_ne_bytes(padded[len..len + WORD_BYTES].try_into().unwrap());
        debug_assert_eq!(stored_cur, self.cur, "PackedDNA::pad must be called first");
        padded
    }

    /// Returns a [`PackedSeq`] view without copying, [`pad`](Self::pad) must be called after modifying it.
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    pub fn as_packed_seq(&self) -> PackedSeq<'_> {
        PackedSeq::from_raw_parts(self.padded_bytes(), 0, self.len())
    }

    /// Converts into a [`PackedSeqVec`] without copying.
    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    pub fn into_packed_seq_vec(mut self) -> PackedSeqVec {
        self.pad();
        let len = self.len();
        let padded_len = self.bits.len() + MARGIN * WORD_BYTES;
        let mut bits = core::mem::take(&mut self.bits);
        // Safety: `pad` initialized `cur` and the padding, `from_raw_parts` then only truncates.
        unsafe { bits.set_len(padded_len) };
        PackedSeqVec::from_raw_parts(bits, len)
    }

    /// Reuses the buffer of a [`PackedSeqVec`], keeping its content.
    #[cfg(feature = "packed-seq")]
    pub fn from_packed_seq_vec(seq: PackedSeqVec) -> Self {
        let len = seq.len();
        let mut bits = seq.into_raw();
        let full_bytes = len / BP_PER_BLOCK * WORD_BYTES;
        let mut last = [0; WORD_BYTES];
        last[..bits.len() - full_bytes].copy_from_slice(&bits[full_bytes..]);
        bits.truncate(full_bytes);
        let rem = 2 * len % BITS_PER_BLOCK;
        Self {
            init_bytes: bits.len(),
            bits,
            cur: T::from_ne_bytes(last) & ((1 << rem) - 1),
            num_bits: 2 * len,
        }
    }
}

impl fmt::Display for PackedDNA {
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

    fn encode(s: &str) -> (u128, usize) {
        let mut packed: u128 = 0;
        for (i, &b) in s.as_bytes().iter().enumerate() {
            let code: u128 = match b {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'T' | b't' => 0b10,
                b'G' | b'g' => 0b11,
                _ => 0,
            };
            packed |= code << (2 * i);
        }
        (packed, s.len())
    }

    fn roundtrip(s: &str) -> String {
        let mut dna = PackedDNA::new();
        let (packed, num_bases) = encode(s);
        dna.append(packed, num_bases);
        format!("{dna}")
    }

    #[test]
    fn append_full_block_at_boundary() {
        let (a_packed, _) = encode(&"A".repeat(BP_PER_BLOCK));
        let (g_packed, _) = encode(&"G".repeat(BP_PER_BLOCK));

        let mut dna = PackedDNA::new();
        dna.append(a_packed, BP_PER_BLOCK);
        dna.append(g_packed, BP_PER_BLOCK);

        assert_eq!(dna.len(), 2 * BP_PER_BLOCK);
        for i in 0..BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b00, "base {i} should be A");
        }
        for i in BP_PER_BLOCK..2 * BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b11, "base {i} should be G");
        }
    }

    #[test]
    fn three_full_blocks_distinct() {
        let (a, _) = encode(&"A".repeat(BP_PER_BLOCK));
        let (c, _) = encode(&"C".repeat(BP_PER_BLOCK));
        let (t, _) = encode(&"T".repeat(BP_PER_BLOCK));

        let mut dna = PackedDNA::new();
        dna.append(a, BP_PER_BLOCK);
        dna.append(c, BP_PER_BLOCK);
        dna.append(t, BP_PER_BLOCK);

        assert_eq!(dna.len(), 3 * BP_PER_BLOCK);
        for i in 0..BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b00, "base {i} A");
        }
        for i in BP_PER_BLOCK..2 * BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b01, "base {i} C");
        }
        for i in 2 * BP_PER_BLOCK..3 * BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b10, "base {i} T");
        }
    }

    #[test]
    fn partial_then_full_block() {
        let partial = BP_PER_BLOCK - 6;
        let (a_partial, a_len) = encode(&"A".repeat(partial));
        let (g_full, _) = encode(&"G".repeat(BP_PER_BLOCK));

        let mut dna = PackedDNA::new();
        dna.append(a_partial, a_len);
        dna.append(g_full, BP_PER_BLOCK);

        assert_eq!(dna.len(), partial + BP_PER_BLOCK);
        for i in 0..partial {
            assert_eq!(dna.get(i), 0b00, "base {i} A");
        }
        for i in partial..partial + BP_PER_BLOCK {
            assert_eq!(dna.get(i), 0b11, "base {i} G");
        }
    }

    #[test]
    fn roundtrip_short() {
        assert_eq!(roundtrip("ACGT"), "ACGT");
    }

    #[test]
    fn roundtrip_exactly_one_block() {
        let s = "ACGT".repeat(BP_PER_BLOCK / 4);
        assert_eq!(s.len(), BP_PER_BLOCK);
        assert_eq!(roundtrip(&s), s);
    }

    #[test]
    fn roundtrip_two_full_blocks() {
        let block = "ACGT".repeat(BP_PER_BLOCK / 4);
        let mut dna = PackedDNA::new();
        let (p, n) = encode(&block);
        dna.append(p, n);
        let (p, n) = encode(&block);
        dna.append(p, n);
        assert_eq!(dna.len(), 2 * BP_PER_BLOCK);
        assert_eq!(format!("{dna}"), block.repeat(2));
    }

    #[test]
    fn append_full_aligned_and_unaligned() {
        let (a_full, _) = encode(&"A".repeat(64));
        let (g_full, _) = encode(&"G".repeat(64));
        let (c_part, c_len) = encode("CCC");

        let mut aligned = PackedDNA::new();
        aligned.append_full(a_full);
        aligned.append_full(g_full);
        assert_eq!(format!("{aligned}"), "A".repeat(64) + &"G".repeat(64));

        let mut unaligned = PackedDNA::new();
        unaligned.append(c_part, c_len);
        unaligned.append_full(a_full);
        unaligned.append_full(g_full);
        assert_eq!(
            format!("{unaligned}"),
            "CCC".to_string() + &"A".repeat(64) + &"G".repeat(64)
        );
    }

    #[test]
    fn test_pop_base() {
        let (a_full, _) = encode(&"A".repeat(64));
        let mut dna = PackedDNA::new();
        dna.append_full(a_full);
        assert_eq!(dna.len(), 64);
        dna.pop_base();
        assert_eq!(dna.len(), 63);
        assert_eq!(format!("{dna}"), "A".repeat(63));

        // Pop across block boundary
        let (g_full, _) = encode(&"G".repeat(64));
        dna.append_full(g_full);
        assert_eq!(dna.len(), 63 + 64);
        dna.pop_base();
        assert_eq!(dna.len(), 63 + 63);
        assert_eq!(format!("{dna}"), "A".repeat(63) + &"G".repeat(63));
    }

    #[cfg(feature = "packed-seq")]
    #[test]
    fn test_packed_seq_padding() {
        use packed_seq::SeqVec;
        let seq = "ACGT".repeat(20) + "TTG";
        let mut dna = PackedDNA::new();
        let (packed, num_bases) = encode(&seq[..60]);
        dna.append(packed, num_bases);
        let (packed, num_bases) = encode(&seq[60..]);
        dna.append(packed, num_bases);
        dna.pad();
        let to_string = |s: PackedSeq| String::from_utf8(s.unpack()).unwrap();
        assert_eq!(to_string(dna.as_packed_seq()), seq);
        let mut copy = dna.clone();
        copy.pad();
        assert_eq!(to_string(copy.as_packed_seq()), seq);
        dna.pop_base();
        dna.pad();
        assert_eq!(to_string(dna.as_packed_seq()), seq[..seq.len() - 1]);
        assert_eq!(
            to_string(dna.into_packed_seq_vec().as_slice()),
            seq[..seq.len() - 1]
        );
        assert_eq!(to_string(PackedDNA::new().as_packed_seq()), "");
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic]
    fn append_rejects_bit_count() {
        PackedDNA::new().append(!0, 128);
    }

    #[cfg(all(feature = "packed-seq", debug_assertions))]
    #[test]
    #[should_panic]
    fn test_packed_seq_requires_padding() {
        let mut dna = PackedDNA::new();
        dna.pad();
        dna.append(0b11, 1);
        dna.as_packed_seq();
    }
}
