use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;

use core::mem::swap;
use core::ops::Range;
use std::io;

enum State {
    Start,
    Restart,
    Header,
    StartDNA,
    InDNABlock,
    EndDNA,
}

/// A parser for the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
/// Configured via [`ParserOptions`].
///
/// # Example
///
/// ```rust,no_run
/// use helicase::input::*;
/// use helicase::*;
///
/// // set the options of the parser (at compile-time)
/// const CONFIG: Config = ParserOptions::default().config();
///
/// fn main() {
///     let path = "...";
///
///     // create a parser with the desired options
///     let mut parser = FastaParser::<CONFIG, _>::from_file(&path).expect("Cannot open file");
///
///     // iterate over records
///     while let Some(_event) = parser.next() {
///         let header = parser.get_header();
///         let seq = parser.get_dna_string();
///         // ...
///     }
/// }
/// ```
pub struct FastaParser<'a, const CONFIG: Config, I: InputData<'a>> {
    lexer: FastaLexer<'a, CONFIG, I>,
    finished: bool,
    state: State,
    block: FastaChunk,
    block_counter: usize,
    pos_in_block: usize,
    header_range: Range<usize>,
    dna_range: Range<usize>,
    record_start: usize,
    contiguous_dna: bool,
    cur_dna_string: Vec<u8>,
    cur_dna_columnar: ColumnarDNA,
    cur_dna_packed: PackedDNA,
    cur_mask_non_actg: BitMask,
    cur_mask_n: BitMask,
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
            record_start: 0,
            contiguous_dna: true,
            cur_dna_string: Vec::new(),
            cur_dna_columnar: ColumnarDNA::new(),
            cur_dna_packed: PackedDNA::new(),
            cur_mask_non_actg: BitMask::new(),
            cur_mask_n: BitMask::new(),
            dna_len: 0,
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastaParser<'a, CONFIG, I>
{
    fn from_input(input: I) -> io::Result<Self> {
        let lexer = FastaLexer::from_input(input)?;
        if lexer.input.first_byte() != b'>' {
            return Err(io::Error::other("Invalid record start, expected '>'"));
        }
        Ok(Self::from_lexer(lexer))
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> HelicaseParser for FastaParser<'a, CONFIG, I> {
    #[inline(always)]
    fn format(&self) -> Format {
        Format::Fasta
    }

    #[inline(always)]
    fn clear_record(&mut self) {
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
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.cur_mask_non_actg.clear();
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.cur_mask_n.clear();
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
            let off = self.lexer.input.buffer_offset();
            &self.lexer.input.buffer()[self.header_range.start - off..self.header_range.end - off]
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
            let off = self.lexer.input.buffer_offset();
            self.lexer.input.buffer()[self.header_range.start - off..self.header_range.end - off]
                .to_vec()
        }
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        let zero_copy = self.contiguous_dna
            && !(flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS));
        if I::RANDOM_ACCESS && zero_copy {
            &self.lexer.input.data()[self.dna_range.clone()]
        } else if !I::RANDOM_ACCESS && zero_copy {
            let buf_off = self.lexer.input.buffer_offset();
            &self.lexer.input.buffer()[self.dna_range.start - buf_off..self.dna_range.end - buf_off]
        } else {
            &self.cur_dna_string
        }
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        let zero_copy = self.contiguous_dna
            && !(flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS));
        if I::RANDOM_ACCESS && zero_copy {
            self.lexer.input.data()[self.dna_range.clone()].to_vec()
        } else if !I::RANDOM_ACCESS && zero_copy {
            let buf_off = self.lexer.input.buffer_offset();
            self.lexer.input.buffer()[self.dna_range.start - buf_off..self.dna_range.end - buf_off]
                .to_vec()
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
    fn get_mask_non_actg(&self) -> &BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            panic!("Parser config error: mask_non_actg is not enabled")
        }
        &self.cur_mask_non_actg
    }

    #[inline(always)]
    fn get_mask_non_actg_owned(&mut self) -> BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            panic!("Parser config error: mask_non_actg is not enabled")
        }
        let mut res = BitMask::with_capacity(self.cur_mask_non_actg.capacity());
        swap(&mut res, &mut self.cur_mask_non_actg);
        res
    }

    #[inline(always)]
    fn get_mask_n(&self) -> &BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
            panic!("Parser config error: mask_n is not enabled")
        }
        &self.cur_mask_n
    }

    #[inline(always)]
    fn get_mask_n_owned(&mut self) -> BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
            panic!("Parser config error: mask_n is not enabled")
        }
        let mut res = BitMask::with_capacity(self.cur_mask_n.capacity());
        swap(&mut res, &mut self.cur_mask_n);
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
            if self.contiguous_dna
                && !(flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS))
            {
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

    /// Keep the record anchor in the buffer when the buffer is about to be
    /// exhausted inside a section loop.  CONFIG predicates come first so the
    /// call compiles away when not needed.  Uses `contiguous_dna` to decide
    /// whether zero-copy DNA is still active (FASTA can switch mid-record).
    #[inline(always)]
    fn make_room_record(&mut self) {
        if (flag_is_set(CONFIG, COMPUTE_HEADER) || flag_is_set(CONFIG, COMPUTE_DNA_STRING))
            && !I::RANDOM_ACCESS
            && self.lexer.input.is_end_of_buffer()
            && (flag_is_set(CONFIG, COMPUTE_HEADER) || self.contiguous_dna)
        {
            self.lexer
                .input
                .make_room(if flag_is_set(CONFIG, COMPUTE_HEADER) {
                    self.header_range.start
                } else {
                    self.dna_range.start
                });
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
        while position == 0 {
            if flag_is_set(CONFIG, COMPUTE_HEADER)
                && !I::RANDOM_ACCESS
                && self.lexer.input.is_end_of_buffer()
            {
                // Header-only anchor: dna_range is not yet set at this point.
                self.lexer.input.make_room(self.header_range.start);
            }
            self.block = match self.lexer.next() {
                Some(b) => b,
                None => {
                    return true;
                }
            };
            self.block_counter += 1;
            self.pos_in_block = 0;
            position = !self.block.header;
        }
        self.pos_in_block = position.trailing_zeros() as usize;
        false
    }

    #[inline(always)]
    fn skip_to_non_dna(&mut self) -> bool {
        // Skip per-block copies when the zero-copy fast path is active: contiguous_dna.
        let mask = !0 << self.pos_in_block;
        let mut position = !self.block.is_dna & mask;
        let mut first_pos = self.pos_in_block;
        while position == 0 {
            if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                && (!self.contiguous_dna
                    || (flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && flag_is_set(CONFIG, MERGE_DNA_CHUNKS)))
            {
                let dna_chunk = &self.lexer.input().current_block()[self.pos_in_block..];
                self.cur_dna_string.extend_from_slice(dna_chunk);
            }
            if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                self.cur_dna_columnar.append(
                    self.block.high_bit >> self.pos_in_block,
                    self.block.low_bit >> self.pos_in_block,
                    self.block.len - self.pos_in_block,
                );
            }
            if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                self.cur_dna_packed.append(
                    self.block.two_bits >> (2 * self.pos_in_block),
                    2 * (self.block.len - self.pos_in_block),
                );
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                self.cur_mask_non_actg.append(
                    self.block.mask_non_actg >> self.pos_in_block,
                    self.block.len - self.pos_in_block,
                );
            }
            if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                self.cur_mask_n.append(
                    self.block.mask_n >> self.pos_in_block,
                    self.block.len - self.pos_in_block,
                );
            }
            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                self.dna_len += self.block.len - self.pos_in_block;
            }
            self.make_room_record();
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
        if flag_is_set(CONFIG, COMPUTE_DNA_STRING)
            && (!self.contiguous_dna
                || (flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS)))
        {
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
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.cur_mask_non_actg.append(
                self.block.mask_non_actg >> first_pos,
                self.pos_in_block - first_pos,
            );
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.cur_mask_n.append(
                self.block.mask_n >> first_pos,
                self.pos_in_block - first_pos,
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
            // Keep the record anchor in the buffer if the zero-copy fast path is active.
            self.make_room_record();
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
                            return Some(Event::Record(self.record_start..self.global_pos()));
                        }
                        continue;
                    }
                    self.finished = self.skip_to_header_or_dna();
                    if (1u64 << self.pos_in_block & self.block.header) != 0 {
                        self.state = State::Header;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            self.prepare_return();
                            return Some(Event::Record(self.record_start..self.global_pos()));
                        }
                    } else if (1u64 << self.pos_in_block & self.block.is_dna) != 0 {
                        self.state = State::StartDNA;
                    }
                }
                State::Header => {
                    if flag_is_not_set(CONFIG, MERGE_RECORDS) {
                        self.clear_record();
                    }
                    self.record_start = self.global_pos();
                    if flag_is_set(CONFIG, COMPUTE_HEADER) {
                        self.header_range.start = self.global_pos() + 1;
                    }
                    self.finished = self.skip_to_end_header();
                    if flag_is_set(CONFIG, COMPUTE_HEADER) {
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
                    self.dna_range.start = self.global_pos();
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
                        && self.contiguous_dna
                        && !(flag_is_set(CONFIG, SPLIT_NON_ACTG)
                            && flag_is_set(CONFIG, MERGE_DNA_CHUNKS))
                        && ((1 << self.pos_in_block) & self.block.header) == 0
                    // Split DNA (newline or non-ACTG)
                    {
                        if I::RANDOM_ACCESS {
                            let dna_chunk = &self.lexer.input.data()[self.dna_range.clone()];
                            self.cur_dna_string.extend_from_slice(dna_chunk);
                        } else {
                            // Retroactively copy the first segment from the buffer.
                            let buf_off = self.lexer.input.buffer_offset();
                            let start = self.dna_range.start - buf_off;
                            let end = self.dna_range.end - buf_off;
                            self.cur_dna_string
                                .extend_from_slice(&self.lexer.input.buffer()[start..end]);
                        }
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
                        return Some(Event::DnaChunk(self.dna_range.clone()));
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
    const CONFIG_PACKED_KEEP: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .and_dna_packed()
        .keep_non_actg()
        .config();

    static FASTA: &[u8] =
        b">head\nTTTCTtaAAAA\nAGAAAA\nACAAN\n\n>h>h\nCTCTTANNAAA\nCAAAnAGCTTT\n>A B C \nCCAC"
            .as_slice();

    #[test]
    fn test_header() {
        let mut f = FastaParser::<CONFIG_HEADER, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_header_owned()).unwrap());
        }
        assert_eq!(res, vec!["head", "h>h", "A B C ",]);
    }

    #[test]
    fn test_dna_string() {
        let mut f = FastaParser::<CONFIG_STRING, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(String::from_utf8(f.get_dna_string_owned()).unwrap());
        }
        assert_eq!(
            res,
            vec!["TTTCTtaAAAAAGAAAAACAAN", "CTCTTANNAAACAAAnAGCTTT", "CCAC",]
        );
        println!();

        let mut f = FastaParser::<CONFIG_STRING_ACTG, _>::from_slice(FASTA).unwrap();
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

        let mut f = FastaParser::<CONFIG_STRING_ACTG_MERGE, _>::from_slice(FASTA).unwrap();
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
        let mut f = FastaParser::<CONFIG_COLUMNAR, _>::from_slice(FASTA).unwrap();
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

        let mut f = FastaParser::<CONFIG_COLUMNAR_MERGE, _>::from_slice(FASTA).unwrap();
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
        let mut f = FastaParser::<CONFIG_PACKED, _>::from_slice(FASTA).unwrap();
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

        let mut f = FastaParser::<CONFIG_PACKED_MERGE, _>::from_slice(FASTA).unwrap();
        let mut res = Vec::new();
        while let Some(_) = f.next() {
            res.push(format!("{}", f.get_dna_packed_owned()));
        }
        assert_eq!(
            res,
            vec!["TTTCTTAAAAAAGAAAAACAA", "CTCTTAAAACAAAAGCTTT", "CCAC",]
        );
    }

    #[test]
    fn test_packed_matches_string_long_sequence() {
        let seq = b"ACGT".repeat(40);
        let fasta = [b">seq\n".as_ref(), seq.as_ref()].concat();
        let mut f = FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
        while let Some(_) = f.next() {
            let s = f.get_dna_string();
            let p = f.get_dna_packed();
            assert_eq!(p.len(), s.len(), "length mismatch");
            for (i, base) in s.iter().enumerate() {
                assert_eq!(
                    p.get_char(i) as u8,
                    base.to_ascii_uppercase(),
                    "mismatch at base {i}"
                );
            }
        }
    }

    #[test]
    fn test_packed_length_keep_non_actg_last_chunk_all_dna() {
        let seq = b"A".repeat(194);
        let fasta = [b">s\n".as_ref(), seq.as_ref()].concat();
        assert_eq!(fasta.len() % 64, 197 % 64);
        let mut f = FastaParser::<CONFIG_PACKED_KEEP, _>::from_slice(&fasta).unwrap();
        while let Some(_) = f.next() {
            let string_len = f.get_dna_string().len();
            let packed_len = f.get_dna_packed().len();
            assert_eq!(packed_len, string_len, "packed length mismatch");
        }
    }
}
