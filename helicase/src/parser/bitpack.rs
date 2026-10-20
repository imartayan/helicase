use crate::config::{advanced::*, *};
use crate::dna_format::{BitMask, ColumnarDNA, PackedDNA};
use crate::lexer::Chunk;
use crate::simd::extract_dna_bitmask;
use core::mem::swap;

/// Unified container managing on-the-fly 2-bit bitpacking representations
/// across FASTA and FASTQ parsers.
#[derive(Default)]
pub struct DnaBitpackers {
    pub packed: PackedDNA,
    pub columnar: ColumnarDNA,
    pub mask_non_actg: BitMask,
    pub mask_n: BitMask,
}

impl DnaBitpackers {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn clear<const CONFIG: Config>(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.packed.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.columnar.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n.clear();
        }
    }

    #[inline(always)]
    pub fn pop_base<const CONFIG: Config>(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.packed.pop_base();
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.columnar.pop_base();
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg.pop_base();
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n.pop_base();
        }
    }

    #[inline(always)]
    pub fn reserve<const CONFIG: Config>(&mut self, additional: usize) {
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.packed.reserve(additional);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.columnar.reserve(additional);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg.reserve(additional);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n.reserve(additional);
        }
    }

    #[inline(always)]
    pub fn take_packed(&mut self) -> PackedDNA {
        let mut res = PackedDNA::with_capacity(self.packed.capacity());
        swap(&mut res, &mut self.packed);
        res
    }

    #[inline(always)]
    pub fn take_columnar(&mut self) -> ColumnarDNA {
        let mut res = ColumnarDNA::with_capacity(self.columnar.capacity());
        swap(&mut res, &mut self.columnar);
        res
    }

    #[inline(always)]
    pub fn take_mask_non_actg(&mut self) -> BitMask {
        let mut res = BitMask::with_capacity(self.mask_non_actg.capacity());
        swap(&mut res, &mut self.mask_non_actg);
        res
    }

    #[inline(always)]
    pub fn take_mask_n(&mut self) -> BitMask {
        let mut res = BitMask::with_capacity(self.mask_n.capacity());
        swap(&mut res, &mut self.mask_n);
        res
    }

    #[inline(always)]
    pub fn append_block_from<const CONFIG: Config>(&mut self, block: &Chunk) {
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.columnar.append_block(block.high_bit, block.low_bit);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.packed.append_block(block.two_bits);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg.append_block(block.mask_non_actg);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n.append_block(block.mask_n);
        }
    }

    #[inline(always)]
    pub fn append_range_from<const CONFIG: Config>(
        &mut self,
        block: &Chunk,
        start: usize,
        len: usize,
    ) {
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.columnar
                .append(block.high_bit >> start, block.low_bit >> start, len);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.packed.append(block.two_bits >> (2 * start), 2 * len);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.mask_non_actg.append(block.mask_non_actg >> start, len);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.mask_n.append(block.mask_n >> start, len);
        }
    }

    /// Appends a wrapped line of length `line_width` bases.
    /// Returns `false` if `SPLIT_NON_ACTG` is enabled and a non-ACTG character was found.
    #[inline(always)]
    pub fn append_line_stride<const CONFIG: Config>(
        &mut self,
        line_data: &[u8],
        line_width: usize,
    ) -> bool {
        self.pack_contiguous::<CONFIG>(line_data, line_width) == line_width
    }

    /// Bulk packs a contiguous slice of pure DNA bases into the active bitpackers.
    ///
    /// Returns the number of bytes processed. If `SPLIT_NON_ACTG` is enabled and
    /// a non-ACTG character is found, it stops early and returns the bytes processed
    /// up to the non-ACTG boundary.
    #[inline(always)]
    pub fn pack_contiguous<const CONFIG: Config>(&mut self, data: &[u8], rem_len: usize) -> usize {
        let need_packed = flag_is_set(CONFIG, COMPUTE_DNA_PACKED);
        let need_columnar = flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR);
        let need_mask_non_actg = flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG);
        let need_mask_n = flag_is_set(CONFIG, COMPUTE_MASK_N);
        let is_keep = flag_is_not_set(CONFIG, SPLIT_NON_ACTG);

        if !need_packed && !need_columnar && !need_mask_non_actg && !need_mask_n {
            return rem_len;
        }

        self.reserve::<CONFIG>(rem_len);

        let q = rem_len / 64;
        let r = rem_len % 64;
        let mut curr_pos = 0;

        let processed_blocks =
            if need_packed && is_keep && !need_columnar && !need_mask_non_actg && !need_mask_n {
                self.packed.pack_blocks::<CONFIG>(data, q);
                q
            } else {
                let mut count = 0;
                for k in 0..q {
                    let block_start = 64 * k;
                    let block_data = &data[block_start..block_start + 64];
                    let m = extract_dna_bitmask::<CONFIG>(block_data);

                    if !is_keep && m.is_dna != !0 {
                        break;
                    }

                    if need_packed {
                        self.packed.append_block(m.two_bits);
                    }
                    if need_columnar {
                        self.columnar.append_block(m.high_bit, m.low_bit);
                    }
                    if need_mask_non_actg {
                        self.mask_non_actg.append_block(m.mask_non_actg);
                    }
                    if need_mask_n {
                        self.mask_n.append_block(m.mask_n);
                    }
                    count += 1;
                }
                count
            };

        if processed_blocks < q {
            return 64 * processed_blocks;
        }

        curr_pos += 64 * q;

        if r > 0 {
            let rem_data = &data[curr_pos..curr_pos + r];
            let m = if curr_pos + 64 <= data.len() {
                extract_dna_bitmask::<CONFIG>(&data[curr_pos..curr_pos + 64])
            } else {
                let mut tmp = [b'A'; 64];
                tmp[..r].copy_from_slice(rem_data);
                extract_dna_bitmask::<CONFIG>(&tmp)
            };

            let r_mask = (1u64 << r) - 1;
            if !is_keep && (m.is_dna & r_mask) != r_mask {
                return curr_pos;
            }

            if need_packed {
                let r_bits = 2 * r;
                let two_bits_mask = if r_bits == 128 {
                    !0
                } else {
                    (1u128 << r_bits) - 1
                };
                self.packed.append(m.two_bits & two_bits_mask, r_bits);
            }
            if need_columnar {
                self.columnar
                    .append(m.high_bit & r_mask, m.low_bit & r_mask, r);
            }
            if need_mask_non_actg {
                self.mask_non_actg.append(m.mask_non_actg & r_mask, r);
            }
            if need_mask_n {
                self.mask_n.append(m.mask_n & r_mask, r);
            }
            curr_pos += r;
        }

        curr_pos
    }
}
