use crate::config::{advanced::*, *};
use crate::lexer::*;

const UPPERCASE: u8 = 0b11011111;
const LUT_ACTG: [u8; 16] = *b"A_C_T_G_________";

#[inline(always)]
pub fn extract_fasta_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastaBitmask {
    let mut open_bracket = 0;
    let mut line_feeds = 0;
    let mut is_dna = !0;
    let mut two_bits = 0;
    let mut high_bit = 0;
    let mut low_bit = 0;

    for (i, &x) in buf.iter().enumerate().take(64) {
        let bit = 1 << i;
        open_bracket |= if x == b'>' { bit } else { 0 };
        line_feeds |= if x == b'\n' { bit } else { 0 };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit |= ((x & 0b100) as u64) << i.wrapping_sub(2);
            low_bit |= ((x & 0b10) as u64) << i.wrapping_sub(1);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            two_bits |= ((x & 0b110) as u128) << (2 * i).wrapping_sub(1);
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            is_dna |= if (x & UPPERCASE) == LUT_ACTG[(x & 0b110) as usize] {
                bit
            } else {
                0
            };
        }
    }

    FastaBitmask {
        open_bracket,
        line_feeds,
        is_dna,
        two_bits,
        high_bit,
        low_bit,
    }
}

#[inline(always)]
pub fn extract_fastq_bitmask<const CONFIG: Config>(buf: &[u8]) -> FastqBitmask {
    let mut line_feeds = 0;
    let mut is_dna = !0;
    let mut two_bits = 0;
    let mut high_bit = 0;
    let mut low_bit = 0;

    for (i, &x) in buf.iter().enumerate().take(64) {
        let bit = 1 << i;
        line_feeds |= if x == b'\n' { bit } else { 0 };

        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            high_bit |= ((x & 0b100) as u64) << i.wrapping_sub(2);
            low_bit |= ((x & 0b10) as u64) << i.wrapping_sub(1);
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            two_bits |= ((x & 0b110) as u128) << (2 * i).wrapping_sub(1);
        }

        if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
            is_dna |= if (x & UPPERCASE) == LUT_ACTG[(x & 0b110) as usize] {
                bit
            } else {
                0
            };
        }
    }

    FastqBitmask {
        line_feeds,
        is_dna,
        two_bits,
        high_bit,
        low_bit,
    }
}
