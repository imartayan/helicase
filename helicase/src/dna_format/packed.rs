use core::cell::UnsafeCell;
use core::fmt::{self, Write};

#[cfg(feature = "packed-seq")]
use packed_seq::{PackedSeq, PackedSeqVec};

type T = u64;
const BITS_PER_BP: usize = 2;
const BITS_PER_BLOCK: usize = T::BITS as usize;
const BP_PER_BLOCK: usize = BITS_PER_BLOCK / BITS_PER_BP;
#[cfg(feature = "packed-seq")]
const PADDING: usize = 6; // 48 bytes / 8 bytes per u64 = 6 u64s

#[derive(Default)]
pub struct PackedDNA {
    bits: UnsafeCell<Vec<T>>,
    cur: T,
    num_bits: usize,
}

impl Clone for PackedDNA {
    fn clone(&self) -> Self {
        Self {
            bits: UnsafeCell::new(unsafe { &*self.bits.get() }.clone()),
            cur: self.cur,
            num_bits: self.num_bits,
        }
    }
}

impl PackedDNA {
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            bits: UnsafeCell::new(Vec::new()),
            cur: 0,
            num_bits: 0,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: UnsafeCell::new(Vec::with_capacity(capacity)),
            cur: 0,
            num_bits: 0,
        }
    }

    #[inline(always)]
    pub const fn len(&self) -> usize {
        self.num_bits / BITS_PER_BP
    }

    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.num_bits == 0
    }

    #[inline(always)]
    pub fn capacity(&self) -> usize {
        unsafe { (*self.bits.get()).capacity() }
    }

    #[inline(always)]
    pub fn reserve(&mut self, additional_bases: usize) {
        let additional_words = additional_bases / BP_PER_BLOCK + 2;
        self.bits.get_mut().reserve(additional_words);
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.bits.get_mut().clear();
        self.cur = 0;
        self.num_bits = 0;
    }

    #[inline(always)]
    pub fn append_block(&mut self, packed: u128) {
        let lo = packed as u64;
        let hi = (packed >> 64) as u64;
        let rem = self.num_bits % BITS_PER_BLOCK;
        self.num_bits += 128;
        let bits = self.bits.get_mut();
        let len = bits.len();
        if len + 2 > bits.capacity() {
            bits.reserve(128);
        }
        unsafe {
            let ptr = bits.as_mut_ptr().add(len);
            if rem == 0 {
                ptr.write(lo);
                ptr.add(1).write(hi);
            } else {
                let inv = BITS_PER_BLOCK - rem;
                ptr.write(self.cur | (lo << rem));
                ptr.add(1).write((lo >> inv) | (hi << rem));
                self.cur = hi >> inv;
            }
            bits.set_len(len + 2);
        }
    }

    #[inline(always)]
    pub fn bits_mut(&mut self) -> &mut Vec<u64> {
        self.bits.get_mut()
    }

    #[inline(always)]
    pub fn add_bits(&mut self, bits: usize) {
        self.num_bits += bits;
    }

    #[inline(always)]
    pub fn append(&mut self, packed: u128, num_bits: usize) {
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
        let bits = self.bits.get_mut();
        if rem + num_bits >= BITS_PER_BLOCK {
            bits.push(self.cur | (x << rem));
            self.cur = if rem > 0 && num_bits > (BITS_PER_BLOCK - rem) {
                x >> (BITS_PER_BLOCK - rem)
            } else {
                0
            };
        } else {
            self.cur |= x << rem;
        }
    }

    #[inline(always)]
    pub fn bits(&self) -> (&[T], T) {
        let bits = unsafe { &*self.bits.get() };
        (&bits[..self.num_bits / BITS_PER_BLOCK], self.cur)
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        let bits = unsafe { &*self.bits.get() };
        let block_idx = i / BP_PER_BLOCK;
        let bit_offset = 2 * (i % BP_PER_BLOCK);
        if block_idx < bits.len() {
            ((bits[block_idx] >> bit_offset) & 0b11) as u8
        } else {
            ((self.cur >> bit_offset) & 0b11) as u8
        }
    }

    #[inline(always)]
    pub fn get_char(&self, i: usize) -> char {
        const LUT: [char; 4] = ['A', 'C', 'T', 'G'];
        LUT[self.get(i) as usize]
    }

    #[cfg(feature = "packed-seq")]
    #[inline(always)]
    pub(crate) fn ensure_padding(&self) {
        let bits = unsafe { &mut *self.bits.get() };
        let len = self.num_bits / BITS_PER_BLOCK;
        bits.resize(len + 1 + PADDING, 0);
        bits[len] = self.cur;
    }

    #[cfg(feature = "packed-seq")]
    #[allow(clippy::missing_transmute_annotations)]
    #[inline(always)]
    pub fn as_packed_seq(&self) -> PackedSeq<'_> {
        self.ensure_padding();
        let bits = unsafe { &*self.bits.get() };
        let len = self.len();
        let data = bits.as_ptr() as *const u8;
        let data_len = bits.len() * BITS_PER_BLOCK / 8;
        let seq = unsafe { core::slice::from_raw_parts(data, data_len) };
        PackedSeq::from_raw_parts(seq, 0, len)
    }

    #[cfg(feature = "packed-seq")]
    #[allow(
        clippy::missing_transmute_annotations,
        clippy::unsound_collection_transmute
    )]
    #[inline(always)]
    pub fn into_packed_seq_vec(mut self) -> PackedSeqVec {
        self.ensure_padding();
        let len = self.len();
        let mut bits = core::mem::take(self.bits.get_mut());
        let data_len = bits.len() * BITS_PER_BLOCK / 8;
        let data_cap = bits.capacity() * BITS_PER_BLOCK / 8;
        let data = bits.as_mut_ptr() as *mut u8;
        core::mem::forget(bits);
        let seq = unsafe { Vec::from_raw_parts(data, data_len, data_cap) };
        PackedSeqVec::from_raw_parts(seq, len)
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
        (packed, s.len() * 2)
    }

    fn roundtrip(s: &str) -> String {
        let mut dna = PackedDNA::new();
        let (packed, num_bits) = encode(s);
        dna.append(packed, num_bits);
        format!("{dna}")
    }

    #[test]
    fn append_full_block_at_boundary() {
        let (a_packed, _) = encode(&"A".repeat(BP_PER_BLOCK));
        let (g_packed, _) = encode(&"G".repeat(BP_PER_BLOCK));

        let mut dna = PackedDNA::new();
        dna.append(a_packed, BITS_PER_BLOCK);
        dna.append(g_packed, BITS_PER_BLOCK);

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
        dna.append(a, BITS_PER_BLOCK);
        dna.append(c, BITS_PER_BLOCK);
        dna.append(t, BITS_PER_BLOCK);

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
        let (a_partial, a_bits) = encode(&"A".repeat(partial));
        let (g_full, _) = encode(&"G".repeat(BP_PER_BLOCK));

        let mut dna = PackedDNA::new();
        dna.append(a_partial, a_bits);
        dna.append(g_full, BITS_PER_BLOCK);

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
    fn append_block_aligned_and_unaligned() {
        let (a_full, _) = encode(&"A".repeat(64));
        let (g_full, _) = encode(&"G".repeat(64));
        let (c_part, c_bits) = encode("CCC");

        let mut aligned = PackedDNA::new();
        aligned.append_block(a_full);
        aligned.append_block(g_full);
        assert_eq!(format!("{aligned}"), "A".repeat(64) + &"G".repeat(64));

        let mut unaligned = PackedDNA::new();
        unaligned.append(c_part, c_bits);
        unaligned.append_block(a_full);
        unaligned.append_block(g_full);
        assert_eq!(
            format!("{unaligned}"),
            "CCC".to_string() + &"A".repeat(64) + &"G".repeat(64)
        );
    }
}
