use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;

use core::mem::swap;
use core::ops::Range;

/// A parser for the [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).
pub struct FastqParser<'a, const CONFIG: Config, I: InputData<'a>> {
    lexer: FastqLexer<'a, CONFIG, I>,
    finished: bool,
    line_count: usize,
    block: FastqChunk,
    block_counter: usize,
    pos_in_block: usize,
    header_range: Range<usize>,
    quality_range: Range<usize>,
    dna_range: Range<usize>,
    cur_header: Vec<u8>,
    cur_quality: Vec<u8>,
    cur_dna_string: Vec<u8>,
    cur_dna_columnar: ColumnarDNA,
    cur_dna_packed: PackedDNA,
    dna_len: usize,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    fn from_lexer(mut lexer: FastqLexer<'a, CONFIG, I>) -> Self {
        let mut finished: bool = false;
        let block = match lexer.next() {
            Some(c) => c,
            None => {
                finished = true;
                FastqChunk::default()
            }
        };
        Self {
            lexer,
            finished,
            line_count: 0,
            block,
            block_counter: 0,
            pos_in_block: 0,
            header_range: 0..0,
            quality_range: 0..0,
            dna_range: 0..0,
            cur_header: Vec::new(),
            cur_quality: Vec::new(),
            cur_dna_string: Vec::new(),
            cur_dna_columnar: ColumnarDNA::new(),
            cur_dna_packed: PackedDNA::new(),
            dna_len: 0,
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastqParser<'a, CONFIG, I>
{
    fn from_input(input: I) -> Self {
        Self::from_lexer(FastqLexer::from_input(input))
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Parser for FastqParser<'a, CONFIG, I> {
    #[inline(always)]
    fn format(&self) -> Format {
        Format::Fastq
    }

    #[inline(always)]
    fn clear_record(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.cur_header.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_QUALITY) {
            self.cur_quality.clear();
        }
        self.clear_chunk();
    }

    #[inline(always)]
    fn clear_chunk(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
            self.cur_dna_string.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.cur_dna_columnar.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.cur_dna_packed.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len = 0;
        }
    }

    #[inline(always)]
    fn get_header(&self) -> &[u8] {
        assert!(flag_is_set(CONFIG, COMPUTE_HEADER));
        if I::RANDOM_ACCESS {
            &self.lexer.input.data()[self.header_range.clone()]
        } else {
            &self.cur_header
        }
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        assert!(flag_is_set(CONFIG, COMPUTE_HEADER));
        if I::RANDOM_ACCESS {
            self.lexer.input.data()[self.header_range.clone()].to_vec()
        } else {
            let mut res = Vec::with_capacity(self.cur_header.capacity());
            swap(&mut res, &mut self.cur_header);
            res
        }
    }

    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        assert!(flag_is_set(CONFIG, COMPUTE_QUALITY));
        if I::RANDOM_ACCESS {
            Some(&self.lexer.input.data()[self.quality_range.clone()])
        } else {
            Some(&self.cur_quality)
        }
    }

    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        assert!(flag_is_set(CONFIG, COMPUTE_QUALITY));
        if I::RANDOM_ACCESS {
            Some(self.lexer.input.data()[self.quality_range.clone()].to_vec())
        } else {
            let mut res = Vec::with_capacity(self.cur_quality.capacity());
            swap(&mut res, &mut self.cur_quality);
            Some(res)
        }
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_STRING));
        if I::RANDOM_ACCESS && flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
            &self.lexer.input.data()[self.dna_range.clone()]
        } else {
            &self.cur_dna_string
        }
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_STRING));
        if I::RANDOM_ACCESS && flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
            self.lexer.input.data()[self.dna_range.clone()].to_vec()
        } else {
            let mut res = Vec::with_capacity(self.cur_dna_string.capacity());
            swap(&mut res, &mut self.cur_dna_string);
            res
        }
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR));
        &self.cur_dna_columnar
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR));
        let mut res = ColumnarDNA::with_capacity(self.cur_dna_columnar.capacity());
        swap(&mut res, &mut self.cur_dna_columnar);
        res
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_PACKED));
        &self.cur_dna_packed
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_PACKED));
        let mut res = PackedDNA::with_capacity(self.cur_dna_packed.capacity());
        swap(&mut res, &mut self.cur_dna_packed);
        res
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        assert!(flag_is_set(
            CONFIG,
            COMPUTE_DNA_LEN | COMPUTE_DNA_STRING | COMPUTE_DNA_COLUMNAR | COMPUTE_DNA_PACKED
        ));
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len
        } else if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
            if I::RANDOM_ACCESS && flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                self.dna_range.len()
            } else {
                self.cur_dna_string.len()
            }
        } else if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.cur_dna_columnar.len()
        } else if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.cur_dna_packed.len()
        } else {
            unreachable!()
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    #[inline(always)]
    const fn global_pos(&self) -> usize {
        64 * self.block_counter + self.pos_in_block
    }

    #[inline(always)]
    fn increment_pos(&mut self) {
        if self.pos_in_block + 1 < self.lexer.input.current_chunk_len() {
            self.pos_in_block += 1;
        } else {
            match self.lexer.next() {
                Some(b) => {
                    self.block = b;
                    self.block_counter += 1;
                    self.pos_in_block = 0;
                }
                None => self.finished = true,
            };
        }
    }

    #[inline(always)]
    fn consume_newline(&mut self) {
        self.block.newline &= self.block.newline.wrapping_sub(1);
        self.increment_pos();
        self.line_count += 1;
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastqParser<'a, CONFIG, I> {
    type Item = Event;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.line_count % 4 {
                0 => {
                    // HEADER
                    self.increment_pos();
                    if self.finished {
                        return None;
                    }
                    if flag_is_not_set(CONFIG, MERGE_RECORDS) {
                        self.clear_record();
                    }
                    if flag_is_set(CONFIG, COMPUTE_HEADER) && I::RANDOM_ACCESS {
                        self.header_range.start = self.global_pos();
                    }
                    let mut first_pos = self.pos_in_block;
                    while self.block.newline == 0 {
                        if flag_is_set(CONFIG, COMPUTE_HEADER) && !I::RANDOM_ACCESS {
                            let header_chunk = &self.lexer.input.current_chunk()[first_pos..];
                            self.cur_header.extend_from_slice(header_chunk);
                        }
                        self.block = match self.lexer.next() {
                            Some(b) => b,
                            None => {
                                self.finished = true;
                                return None;
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                        first_pos = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    if flag_is_set(CONFIG, COMPUTE_HEADER) && I::RANDOM_ACCESS {
                        self.header_range.end = self.global_pos();
                    }
                    if flag_is_set(CONFIG, COMPUTE_HEADER) && !I::RANDOM_ACCESS {
                        let header_chunk =
                            &self.lexer.input.current_chunk()[first_pos..self.pos_in_block];
                        self.cur_header.extend_from_slice(header_chunk);
                    }
                    self.consume_newline();
                }
                1 => {
                    // SEQUENCE
                    if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                        // skip to dna or newline
                        let mask = !0 << self.pos_in_block;
                        let mut position = (self.block.is_dna | self.block.newline) & mask;
                        while position == 0 {
                            self.block = match self.lexer.next() {
                                Some(b) => b,
                                None => {
                                    self.finished = true;
                                    return None;
                                }
                            };
                            self.block_counter += 1;
                            self.pos_in_block = 0;
                            position = self.block.is_dna | self.block.newline;
                        }
                        self.pos_in_block = position.trailing_zeros() as usize;
                        if ((1 << self.pos_in_block) & self.block.newline) != 0 {
                            self.consume_newline();
                            continue;
                        }
                    }

                    // skip to non dna
                    let mask = !0 << self.pos_in_block;
                    let mut position = !self.block.is_dna & mask;
                    if flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS) {
                        self.clear_chunk();
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                        && I::RANDOM_ACCESS
                    {
                        self.dna_range.start = self.global_pos();
                    }
                    let mut first_pos = self.pos_in_block;
                    while position == 0 {
                        if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                            && (flag_is_set(CONFIG, SPLIT_NON_ACTG) || !I::RANDOM_ACCESS)
                        {
                            let dna_chunk = &self.lexer.input.current_chunk()[self.pos_in_block..];
                            self.cur_dna_string.extend_from_slice(dna_chunk);
                        }
                        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                            self.cur_dna_columnar.append(
                                self.block.high_bit >> self.pos_in_block,
                                self.block.low_bit >> self.pos_in_block,
                                64 - self.pos_in_block,
                            );
                        }
                        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                            self.cur_dna_packed.append(
                                self.block.two_bits >> (2 * self.pos_in_block),
                                128 - 2 * self.pos_in_block,
                            );
                        }
                        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                            self.dna_len +=
                                self.lexer.input.current_chunk_len() - self.pos_in_block;
                        }
                        self.block = match self.lexer.next() {
                            Some(b) => b,
                            None => {
                                self.finished = true;
                                return None;
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                        first_pos = 0;
                        position = !self.block.is_dna;
                    }
                    self.pos_in_block = position.trailing_zeros() as usize;
                    if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && (flag_is_set(CONFIG, SPLIT_NON_ACTG) || !I::RANDOM_ACCESS)
                    {
                        let dna_chunk =
                            &self.lexer.input.current_chunk()[first_pos..self.pos_in_block];
                        self.cur_dna_string.extend_from_slice(dna_chunk);
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                        self.cur_dna_columnar.append(
                            self.block.high_bit >> first_pos,
                            self.block.low_bit >> first_pos,
                            self.pos_in_block - first_pos,
                        );
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                        self.cur_dna_packed.append(
                            self.block.two_bits >> (2 * first_pos),
                            2 * (self.pos_in_block - first_pos),
                        );
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                        self.dna_len += self.pos_in_block - first_pos;
                    }
                    let return_pos = if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                        self.global_pos()
                    } else {
                        0
                    };
                    if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                        && I::RANDOM_ACCESS
                    {
                        self.dna_range.end = self.global_pos();
                    }
                    if flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                        || ((1 << self.pos_in_block) & self.block.newline) != 0
                    {
                        self.consume_newline();
                    }
                    if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                        return Some(Event::DnaChunk(return_pos));
                    }
                }
                2 => {
                    // PLUS
                    while self.block.newline == 0 {
                        self.block = match self.lexer.next() {
                            Some(b) => b,
                            None => {
                                self.finished = true;
                                return None;
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    self.consume_newline();
                }
                3 => {
                    // QUALITY
                    if flag_is_set(CONFIG, COMPUTE_QUALITY) && I::RANDOM_ACCESS {
                        self.quality_range.start = self.global_pos();
                    }
                    let mut first_pos = self.pos_in_block;
                    while self.block.newline == 0 {
                        if flag_is_set(CONFIG, COMPUTE_QUALITY) && !I::RANDOM_ACCESS {
                            let quality_chunk =
                                &self.lexer.input.current_chunk()[self.pos_in_block..];
                            self.cur_quality.extend_from_slice(quality_chunk);
                        }
                        self.block = match self.lexer.next() {
                            Some(b) => b,
                            None => {
                                self.finished = true;
                                break; // return record
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                        first_pos = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    if flag_is_set(CONFIG, COMPUTE_QUALITY) {
                        self.pos_in_block =
                            self.pos_in_block.min(self.lexer.input.current_chunk_len());
                    }
                    if flag_is_set(CONFIG, COMPUTE_QUALITY) && I::RANDOM_ACCESS {
                        self.quality_range.end = self.global_pos();
                    }
                    if flag_is_set(CONFIG, COMPUTE_QUALITY)
                        && !I::RANDOM_ACCESS
                        && self.block.newline != 0
                    {
                        let quality_chunk =
                            &self.lexer.input.current_chunk()[first_pos..self.pos_in_block];
                        self.cur_quality.extend_from_slice(quality_chunk);
                    }
                    self.consume_newline();
                    if flag_is_set(CONFIG, RETURN_RECORD) {
                        return Some(Event::Record(self.global_pos()));
                    }
                }
                _ => unreachable!(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONFIG_HEADER: Config = ParserOptions::default().ignore_dna().config();
    const CONFIG_QUALITY: Config = ParserOptions::default()
        .ignore_headers()
        .ignore_dna()
        .compute_quality()
        .config();
    const CONFIG_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .config();
    const CONFIG_STRING_ACTG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .config()
        & !RETURN_RECORD;
    const CONFIG_STRING_ACTG_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .skip_non_actg()
        .config();
    const CONFIG_COLUMNAR: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .config()
        & !RETURN_RECORD;
    const CONFIG_COLUMNAR_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .config()
        & !RETURN_RECORD;
    const CONFIG_PACKED_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .skip_non_actg()
        .config();

    static FASTQ: &[u8] =
        b"@head\nTTTCTtaAAAAAGAAAAACAAN\n+\n123\n@hhh\nCTCTTANNAAACAAAnAGCTTT\n+\nQQ@@++AA\n@A B C \nCCAC\n+\nQUAL"
            .as_slice();

    #[test]
    fn test_header() {
        let mut f = FastqParser::<CONFIG_HEADER, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        let mut c = 0;
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
            c += 1;
            if c > 3 {
                break;
            }
        }
        assert_eq!(res, vec!["head", "hhh", "A B C "]);
    }

    #[test]
    fn test_quality() {
        let mut f = FastqParser::<CONFIG_QUALITY, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        let mut c = 0;
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_quality_owned().unwrap()).unwrap());
            c += 1;
            if c > 3 {
                break;
            }
        }
        assert_eq!(res, vec!["123", "QQ@@++AA", "QUAL"]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = FastqParser::<CONFIG_STRING, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC"]
        );

        let mut f = FastqParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec![
                "TTTCTtaAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC"]
        );
    }

    #[test]
    fn test_dna_columnar() {
        let mut f = FastqParser::<CONFIG_COLUMNAR, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA", // uppercased
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA", // uppercased
                "CTCTTAAAACAAAAGCTTT",
                "CCAC"
            ]
        );
    }

    #[test]
    fn test_dna_packed() {
        let mut f = FastqParser::<CONFIG_PACKED, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA", // uppercased
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC"
            ]
        );

        let mut f = FastqParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTQ);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA", // uppercased
                "CTCTTAAAACAAAAGCTTT",
                "CCAC"
            ]
        );
    }
}
