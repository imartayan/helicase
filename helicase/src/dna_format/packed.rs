use std::fmt::{self, Write};

type T = u128;
const BITS_PER_BP: usize = 2;
const BITS_PER_BLOCK: usize = T::BITS as usize;
const BP_PER_BLOCK: usize = BITS_PER_BLOCK / BITS_PER_BP;
// const PADDING: usize = 3;

#[derive(Clone, Default)]
pub struct PackedDNA {
    bits: Vec<T>,
    cur: T,
    num_bits: usize,
}

impl PackedDNA {
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            bits: Vec::new(),
            cur: 0,
            num_bits: 0,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: Vec::with_capacity(capacity),
            cur: 0,
            num_bits: 0,
        }
    }

    #[inline(always)]
    pub const fn len(&self) -> usize {
        self.num_bits / 2
    }

    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.num_bits == 0
    }

    #[inline(always)]
    pub const fn capacity(&self) -> usize {
        self.bits.capacity()
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.bits.clear();
        self.num_bits = 0;
    }

    #[inline(always)]
    pub fn append(&mut self, packed: T, num_bits: usize) {
        if num_bits == 0 {
            // should not happen?
            return;
        }
        let rem = self.num_bits % BITS_PER_BLOCK;
        let mask = !0 >> (BITS_PER_BLOCK - num_bits);
        self.num_bits += num_bits;
        let x = packed & mask;
        if rem + num_bits >= BITS_PER_BLOCK {
            self.cur |= packed << rem;
            self.bits.push(self.cur);
            self.cur = x >> (BITS_PER_BLOCK - rem);
        } else {
            self.cur |= x << rem;
        }
    }

    #[inline(always)]
    pub fn bits(&self) -> (&[T], T) {
        (&self.bits, self.cur)
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        if i < self.len() & (!0 << BP_PER_BLOCK.trailing_zeros()) {
            ((self.bits[i / BP_PER_BLOCK] >> (2 * (i % BP_PER_BLOCK))) & 0b11) as u8
        } else {
            ((self.cur >> (2 * (i % BP_PER_BLOCK))) & 0b11) as u8
        }
    }

    #[inline(always)]
    pub fn get_char(&self, i: usize) -> char {
        const LUT: [char; 4] = ['A', 'C', 'T', 'G'];
        LUT[self.get(i) as usize]
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
