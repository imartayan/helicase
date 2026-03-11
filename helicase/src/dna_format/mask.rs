type T = u64;
const BITS_PER_BP: usize = 1;
const BITS_PER_BLOCK: usize = T::BITS as usize;
const BP_PER_BLOCK: usize = BITS_PER_BLOCK / BITS_PER_BP;

#[derive(Clone, Default)]
pub struct BitMask {
    bits: Vec<T>,
    cur: T,
    len: usize,
}

impl BitMask {
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            bits: Vec::new(),
            cur: 0,
            len: 0,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: Vec::with_capacity(capacity),
            cur: 0,
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
        self.bits.capacity()
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.bits.clear();
        self.cur = 0;
        self.len = 0;
    }

    #[inline(always)]
    pub fn append(&mut self, x: T, size: usize) {
        if size == 0 {
            // should not happen?
            return;
        }
        let rem = self.len % BITS_PER_BLOCK;
        let mask = !0 >> (BITS_PER_BLOCK - size);
        self.len += size;
        let y = x & mask;
        if rem + size >= BITS_PER_BLOCK {
            self.cur |= x << rem;
            self.bits.push(self.cur);
            self.cur = if rem > 0 {
                y >> (BITS_PER_BLOCK - rem)
            } else {
                0
            };
        } else {
            self.cur |= y << rem;
        }
    }

    #[inline(always)]
    pub fn bits(&self) -> (&[T], T) {
        (&self.bits, self.cur)
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> bool {
        if i < self.len & (!0 << BITS_PER_BLOCK.trailing_zeros()) {
            (self.bits[i / BP_PER_BLOCK] >> (i % BP_PER_BLOCK)) & 1 != 0
        } else {
            (self.cur >> (i % BP_PER_BLOCK)) & 1 != 0
        }
    }

    #[cfg(test)]
    pub fn get_all(&self) -> Vec<bool> {
        (0..self.len).map(|i| self.get(i)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn append_exact_64_bits_bulk() {
        let mut m = BitMask::new();
        m.append(!0u64, 64); // 64 ones
        m.append(0u64, 1); // 1 zero

        assert_eq!(m.len(), 65);
        assert_eq!(m.bits.len(), 1);
        let got = m.get_all();
        assert!(got[..64].iter().all(|&b| b));
        assert!(!got[64]);
    }

    #[test]
    fn two_full_blocks_bulk_distinct() {
        let mut m = BitMask::new();
        m.append(!0u64, 64); // 64 ones
        m.append(0u64, 64); // 64 zeros

        assert_eq!(m.len(), 128);
        let got = m.get_all();
        assert!(got[..64].iter().all(|&b| b));
        assert!(got[64..].iter().all(|&b| !b));
    }

    #[test]
    fn partial_then_full_block() {
        let mut m = BitMask::new();
        m.append(!0u64, 10); // 10 ones, rem becomes 10
        m.append(0u64, 64); // 64 zeros crossing the boundary

        assert_eq!(m.len(), 74);
        let got = m.get_all();
        assert!(got[..10].iter().all(|&b| b));
        assert!(got[10..].iter().all(|&b| !b));
    }
}
