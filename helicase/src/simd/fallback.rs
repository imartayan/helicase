use crate::config::{advanced::*, *};
use crate::lexer::*;

const UPPERCASE: u8 = 0b11011111;
const LUT_ACTG: [u8; 16] = *b"A_C_T_G_________";

#[inline(always)]
pub fn extract_bitmask<const CONFIG: Config, const _MODE: u8>(buf: &[u8]) -> Bitmask {
    let mut line_feeds = 0;
    let mut is_dna = !0;
    let mut two_bits = 0;
    let mut high_bit = 0;
    let mut low_bit = 0;
    let mut mask_actg = 0;
    let mut mask_non_actg = 0;
    let mut mask_n = 0;

    for (i, &x) in buf.iter().enumerate().take(64) {
        let bit = 1 << i;
        line_feeds |= if x == b'\n' { bit } else { 0 };

        let code = (x >> 1) & 0b11;
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit |= ((code >> 1) as u64) << i;
            low_bit |= ((code & 1) as u64) << i;
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            two_bits |= (code as u128) << (2 * i);
        }

        if flag_is_set(
            CONFIG,
            SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG | COMPUTE_MASK_N,
        ) {
            mask_actg |= if (x & UPPERCASE) == LUT_ACTG[(x & 0b110) as usize] {
                bit
            } else {
                0
            };
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                is_dna = mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                mask_non_actg = !mask_actg;
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                mask_n |= if (x & UPPERCASE) == b'N' { bit } else { 0 };
            }
        }
    }

    Bitmask {
        line_feeds,
        is_dna,
        two_bits,
        high_bit,
        low_bit,
        mask_non_actg,
        mask_n,
    }
}

// 2-bit code, whether it is not ACTG, whether it is N, encoded like `extract_bitmask`.
#[inline(always)]
const fn bitpack_byte(x: u8) -> (u8, bool, bool) {
    let is_actg = (x & UPPERCASE) == LUT_ACTG[(x & 0b110) as usize];
    ((x & 0b110) >> 1, !is_actg, (x & UPPERCASE) == b'N')
}

#[inline(always)]
pub fn extract_dna_bitmask<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    let mut is_dna = !0;
    let mut two_bits = 0;
    let mut high_bit = 0;
    let mut low_bit = 0;
    let mut mask_non_actg = 0;
    let mut mask_n = 0;

    for (i, &b) in buf.iter().enumerate() {
        let (actg, non_actg, n) = bitpack_byte(b);
        two_bits |= (actg as u128) << (2 * i);
        high_bit |= ((actg >> 1) as u64) << i;
        low_bit |= ((actg & 1) as u64) << i;
        is_dna &= !((non_actg as u64) << i);
        mask_non_actg |= (non_actg as u64) << i;
        mask_n |= (n as u64) << i;
    }

    Bitmask {
        line_feeds: 0,
        is_dna,
        two_bits,
        high_bit,
        low_bit,
        mask_non_actg,
        mask_n,
    }
}
