use crate::config::{advanced::*, *};
use crate::lexer::*;

const UPPERCASE: u8 = 0b11011111;
const LUT_ACTG: [u8; 16] = *b"A_C_T_G_________";

#[inline(always)]
pub fn load(buf: &[u8]) -> &[u8] {
    buf
}

#[inline(always)]
pub fn line_feeds(buf: &[u8]) -> u64 {
    (buf.iter().take(64).enumerate()).fold(0, |acc, (i, &x)| acc | ((x == b'\n') as u64) << i)
}

#[inline(always)]
pub fn scan_line_feeds(buf: &[u8]) -> u64 {
    line_feeds(buf)
}

#[inline(always)]
pub fn bitpack_dna<const CONFIG: Config>(buf: &[u8]) -> Bitmask {
    let mut mask = Bitmask::ALL_DNA;
    let mut mask_actg = 0;
    for (i, &x) in buf.iter().enumerate().take(64) {
        let code = (x >> 1) & 0b11;
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            mask.high_bit |= ((code >> 1) as u64) << i;
            mask.low_bit |= ((code & 1) as u64) << i;
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            mask.two_bits |= (code as u128) << (2 * i);
        }
        if (x & UPPERCASE) == LUT_ACTG[(x & 0b110) as usize] {
            mask_actg |= 1 << i;
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) && (x & UPPERCASE) == b'N' {
            mask.mask_n |= 1 << i;
        }
    }
    if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
        mask.is_dna = mask_actg;
    }
    if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
        mask.mask_non_actg = !mask_actg;
    }
    mask
}
