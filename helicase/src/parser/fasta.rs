use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;

use core::mem::swap;
use core::ops::Range;

enum State {
    Start,
    Restart,
    Header,
    StartDNA,
    InDNABlock,
    EndDNA,
}

/// A parser for the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
pub struct FastaParser<'a, const CONFIG: Config, I: InputData<'a>> {
    lexer: FastaLexer<'a, CONFIG, I>,
    finished: bool,
    state: State,
    block: FastaChunk,
    block_counter: usize,
    pos_in_block: usize,
    header_range: Range<usize>,
    dna_range: Range<usize>,
    contiguous_dna: bool,
    cur_header: Vec<u8>,
    cur_dna_string: Vec<u8>,
    cur_dna_columnar: ColumnarDNA,
    cur_dna_packed: PackedDNA,
    dna_len: usize,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaParser<'a, CONFIG, I> {
    fn from_lexer(mut lexer: FastaLexer<'a, CONFIG, I>) -> Self {
        let mut finished: bool = false;
        let first = match lexer.next() {
            Some(c) => c,
            None => {
                finished = true;
                FastaChunk::default()
            }
        };
        Self {
            lexer,
            finished,
            state: State::Start,
            block: first,
            block_counter: 0,
            pos_in_block: 0,
            header_range: 0..0,
            dna_range: 0..0,
            contiguous_dna: true,
            cur_header: Vec::new(),
            cur_dna_string: Vec::new(),
            cur_dna_columnar: ColumnarDNA::new(),
            cur_dna_packed: PackedDNA::new(),
            dna_len: 0,
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastaParser<'a, CONFIG, I>
{
    fn from_input(input: I) -> Self {
        Self::from_lexer(FastaLexer::from_input(input))
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Parser for FastaParser<'a, CONFIG, I> {
    #[inline(always)]
    fn format(&self) -> Format {
        Format::Fasta
    }

    #[inline(always)]
    fn clear_record(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.cur_header.clear();
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
        if flag_is_not_set(CONFIG, COMPUTE_HEADER) {
            panic!("Parser config error: headers are ignored")
        }
        if I::RANDOM_ACCESS {
            &self.lexer.input.data()[self.header_range.clone()]
        } else {
            let n = self.cur_header.len();
            if n < 2 {
                &self.cur_header
            } else {
                &self.cur_header[1..(n - 1)]
                // TODO don't store '>' in cur_header
            }
        }
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        if flag_is_not_set(CONFIG, COMPUTE_HEADER) {
            panic!("Parser config error: headers are ignored")
        }
        if I::RANDOM_ACCESS {
            self.lexer.input.data()[self.header_range.clone()].to_vec()
        } else {
            let mut res = Vec::with_capacity(self.cur_header.capacity());
            swap(&mut res, &mut self.cur_header);
            res
        }
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        if I::RANDOM_ACCESS && self.contiguous_dna {
            &self.lexer.input.data()[self.dna_range.clone()]
        } else {
            &self.cur_dna_string
        }
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        if I::RANDOM_ACCESS && self.contiguous_dna {
            self.lexer.input.data()[self.dna_range.clone()].to_vec()
        } else {
            let mut res = Vec::with_capacity(self.cur_dna_string.capacity());
            swap(&mut res, &mut self.cur_dna_string);
            res
        }
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            panic!("Parser config error: dna_columnar is not enabled")
        }
        &self.cur_dna_columnar
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            panic!("Parser config error: dna_columnar is not enabled")
        }
        let mut res = ColumnarDNA::with_capacity(self.cur_dna_columnar.capacity());
        swap(&mut res, &mut self.cur_dna_columnar);
        res
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
            panic!("Parser config error: dna_packed is not enabled")
        }
        &self.cur_dna_packed
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
            panic!("Parser config error: dna_packed is not enabled")
        }
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
            if I::RANDOM_ACCESS && self.contiguous_dna {
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

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaParser<'a, CONFIG, I> {
    #[inline(always)]
    const fn global_pos(&self) -> usize {
        64 * self.block_counter + self.pos_in_block
    }

    #[inline(always)]
    fn prepare_return(&mut self) {
        #[cfg(feature = "packed-seq")]
        {
            if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                self.cur_dna_packed.append_padding();
            }
        }
    }

    #[inline(always)]
    fn skip_to_start_header(&mut self) -> bool {
        let mask = !0 << self.pos_in_block;
        let mut position = self.block.header & mask;
        while position == 0 {
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            position = self.block.header;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        false
    }

    #[inline(always)]
    fn skip_to_header_or_dna(&mut self) -> bool {
        let mask = !0 << self.pos_in_block;
        let mut position = (self.block.is_dna | self.block.header) & mask;
        while position == 0 {
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            position = self.block.is_dna | self.block.header;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        false
    }

    #[inline(always)]
    fn skip_to_end_header(&mut self) -> bool {
        let mask = !0 << self.pos_in_block;
        let mut position = !self.block.header & mask;
        let mut first_pos = self.pos_in_block;
        while position == 0 {
            if flag_is_set(CONFIG, COMPUTE_HEADER) && !I::RANDOM_ACCESS {
                let header_chunk = &self.lexer.input().current_block()[self.pos_in_block..];
                self.cur_header.extend_from_slice(header_chunk);
            }
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            first_pos = 0;
            position = !self.block.header;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        if flag_is_set(CONFIG, COMPUTE_HEADER) && !I::RANDOM_ACCESS {
            let header_chunk = &self.lexer.input().current_block()[first_pos..self.pos_in_block];
            self.cur_header.extend_from_slice(header_chunk);
        }
        false
    }

    #[inline(always)]
    fn skip_to_non_dna(&mut self) -> bool {
        let mask = !0 << self.pos_in_block;
        let mut position = !self.block.is_dna & mask;
        let mut first_pos = self.pos_in_block;
        while position == 0 {
            if flag_is_set(CONFIG, COMPUTE_DNA_STRING) && !(I::RANDOM_ACCESS && self.contiguous_dna)
            {
                let dna_chunk = &self.lexer.input().current_block()[self.pos_in_block..];
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
                self.dna_len += self.block.len - self.pos_in_block;
            }
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    self.pos_in_block = self.block.len;
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            first_pos = 0;
            position = !self.block.is_dna;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        if flag_is_set(CONFIG, COMPUTE_DNA_STRING) && !(I::RANDOM_ACCESS && self.contiguous_dna) {
            let dna_chunk = &self.lexer.input().current_block()[first_pos..self.pos_in_block];
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
        false
    }

    #[inline(always)]
    fn skip_to_dna_or_split_or_header(&mut self) -> bool {
        let mask = !0 << self.pos_in_block;
        let mut position = (self.block.is_dna | self.block.split | self.block.header) & mask;
        while position == 0 {
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    self.pos_in_block = self.block.len;
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            position = self.block.is_dna | self.block.split | self.block.header;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        false
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastaParser<'a, CONFIG, I> {
    type Item = Event;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match &self.state {
                State::Start => {
                    if self.finished {
                        return None;
                    }
                    self.finished = self.skip_to_start_header();
                    if self.block.header != 0 {
                        self.state = State::Header;
                    }
                }
                State::Restart => {
                    if self.finished {
                        self.state = State::Start;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            self.prepare_return();
                            return Some(Event::Record(self.global_pos()));
                        }
                        continue;
                    }
                    self.finished = self.skip_to_header_or_dna();
                    if (1u64 << self.pos_in_block & self.block.header) != 0 {
                        self.state = State::Header;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            self.prepare_return();
                            return Some(Event::Record(self.global_pos()));
                        }
                    } else if (1u64 << self.pos_in_block & self.block.is_dna) != 0 {
                        self.state = State::StartDNA;
                    }
                }
                State::Header => {
                    if flag_is_not_set(CONFIG, MERGE_RECORDS) {
                        self.clear_record();
                    }
                    if flag_is_set(CONFIG, COMPUTE_HEADER) && I::RANDOM_ACCESS {
                        self.header_range.start = self.global_pos() + 1;
                    }
                    self.finished = self.skip_to_end_header();
                    if flag_is_set(CONFIG, COMPUTE_HEADER) && I::RANDOM_ACCESS {
                        self.header_range.end = self.global_pos() - 1;
                    }
                    self.contiguous_dna = true;
                    self.state = State::Restart;
                }
                State::StartDNA => {
                    self.state = State::InDNABlock;
                    if flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS) {
                        self.clear_chunk();
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_STRING) && I::RANDOM_ACCESS {
                        self.dna_range.start = self.global_pos();
                    }
                }
                State::InDNABlock => {
                    if self.skip_to_non_dna() {
                        self.finished = true;
                        self.dna_range.end = self.global_pos();
                        self.state = State::EndDNA;
                        continue;
                    }
                    self.dna_range.end = self.global_pos();
                    if self.skip_to_dna_or_split_or_header() {
                        self.finished = true;
                        self.state = State::EndDNA;
                        continue;
                    }
                    if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && I::RANDOM_ACCESS
                        && self.contiguous_dna
                        && ((1 << self.pos_in_block) & self.block.header) == 0
                    {
                        let dna_chunk = &self.lexer.input.data()[self.dna_range.clone()];
                        self.cur_dna_string.extend_from_slice(dna_chunk);
                        self.contiguous_dna = false;
                    }
                    if ((1 << self.pos_in_block) & self.block.is_dna) != 0 {
                        self.state = State::InDNABlock;
                    } else {
                        self.state = State::EndDNA;
                    }
                }
                State::EndDNA => {
                    self.state = State::Restart;
                    if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                        self.prepare_return();
                        return Some(Event::DnaChunk(self.global_pos() - 1));
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONFIG_HEADER: Config = ParserOptions::default().ignore_dna().config();
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

    static FASTA: &[u8] =
        b">head\nTTTCTtaAAAA\nAGAAAA\nACAAN\n\n>hhh\nCTCTTANNAAA\nCAAAnAGCTTT\n>A B C \nCCAC"
            .as_slice();

    #[test]
    fn test_header() {
        let mut f = FastaParser::<CONFIG_HEADER, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
        }
        assert_eq!(res, vec!["head", "hhh", "A B C ",]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = FastaParser::<CONFIG_STRING, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC",]
        );
        println!();

        let mut f = FastaParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTA);
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
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_dna_columnar() {
        let mut f = FastaParser::<CONFIG_COLUMNAR, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_columnar_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_dna_packed() {
        let mut f = FastaParser::<CONFIG_PACKED, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec![
                "TTTCTTAAAAAAGAAAAACAA",
                "CTCTTA",
                "AAACAAA",
                "AGCTTT",
                "CCAC",
            ]
        );

        let mut f = FastaParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTA);
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }
}
