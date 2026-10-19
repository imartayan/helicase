use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;
use core::mem::swap;
use core::ops::Range;
use std::io;

/// A parser for the [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).
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
///     let mut parser = FastqParser::<CONFIG, _>::from_file(&path).expect("Cannot open file");
///
///     // iterate over records
///     while let Some(_event) = parser.next() {
///         let header = parser.get_header();
///         let seq = parser.get_dna_string();
///         // ...
///     }
/// }
/// ```
pub struct FastqParser<'a, const CONFIG: Config, I: InputData<'a>> {
    lexer: Lexer<'a, CONFIG, I>,
    finished: bool,
    line_count: usize,
    block: Chunk,
    block_counter: usize,
    pos_in_block: usize,
    header_range: Range<usize>,
    quality_range: Range<usize>,
    dna_range: Range<usize>,
    record_start: usize,
    cur_dna_string: Vec<u8>,
    bitpackers: DnaBitpackers,
    dna_len: usize,
    pending_error: Option<io::Error>,
    #[cfg(feature = "paraseq")]
    record_index: u64,
    #[cfg(feature = "paraseq")]
    record_counter: u64,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    fn from_lexer(mut lexer: Lexer<'a, CONFIG, I>) -> Self {
        let mut finished: bool = false;
        let mut pending_error = None;
        let block = match lexer.next() {
            Some(Ok(c)) => c,
            Some(Err(e)) => {
                finished = true;
                pending_error = Some(e);
                Chunk::default()
            }
            None => {
                finished = true;
                Chunk::default()
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
            record_start: 0,
            cur_dna_string: Vec::new(),
            bitpackers: DnaBitpackers::new(),
            dna_len: 0,
            pending_error,
            #[cfg(feature = "paraseq")]
            record_index: 0,
            #[cfg(feature = "paraseq")]
            record_counter: 0,
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastqParser<'a, CONFIG, I>
{
    fn from_input(input: I) -> io::Result<Self> {
        let lexer = Lexer::from_input(input)?;
        if lexer.input.first_byte() != b'@' {
            return Err(io::Error::other("Invalid record start, expected '@'"));
        }
        Ok(Self::from_lexer(lexer))
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> HelicaseParser for FastqParser<'a, CONFIG, I> {
    #[inline(always)]
    fn format(&self) -> Format {
        Format::Fastq
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
        self.bitpackers.clear();
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len = 0;
        }
    }

    #[inline(always)]
    fn get_header(&self) -> &[u8] {
        if flag_is_not_set(CONFIG, COMPUTE_HEADER) {
            panic!("Parser config error: headers are ignored")
        }
        self.lexer.input.slice(self.header_range.clone())
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        self.get_header().to_vec()
    }

    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        if flag_is_not_set(CONFIG, COMPUTE_QUALITY) {
            panic!("Parser config error: quality is ignored")
        }
        Some(self.lexer.input.slice(self.quality_range.clone()))
    }

    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        self.get_quality().map(|q| q.to_vec())
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        if flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS) {
            &self.cur_dna_string
        } else {
            self.lexer.input.slice(self.dna_range.clone())
        }
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        if flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS) {
            let mut res = Vec::with_capacity(self.cur_dna_string.capacity());
            swap(&mut res, &mut self.cur_dna_string);
            res
        } else {
            self.lexer.input.slice(self.dna_range.clone()).to_vec()
        }
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            panic!("Parser config error: dna_columnar is not enabled")
        }
        &self.bitpackers.columnar
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            panic!("Parser config error: dna_columnar is not enabled")
        }
        self.bitpackers.take_columnar()
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
            panic!("Parser config error: dna_packed is not enabled")
        }
        &self.bitpackers.packed
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
            panic!("Parser config error: dna_packed is not enabled")
        }
        self.bitpackers.take_packed()
    }

    #[inline(always)]
    fn get_mask_non_actg(&self) -> &BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            panic!("Parser config error: mask_non_actg is not enabled")
        }
        &self.bitpackers.mask_non_actg
    }

    #[inline(always)]
    fn get_mask_non_actg_owned(&mut self) -> BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            panic!("Parser config error: mask_non_actg is not enabled")
        }
        self.bitpackers.take_mask_non_actg()
    }

    #[inline(always)]
    fn get_mask_n(&self) -> &BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
            panic!("Parser config error: mask_n is not enabled")
        }
        &self.bitpackers.mask_n
    }

    #[inline(always)]
    fn get_mask_n_owned(&mut self) -> BitMask {
        if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
            panic!("Parser config error: mask_n is not enabled")
        }
        self.bitpackers.take_mask_n()
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len
        } else if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
            if flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS) {
                self.cur_dna_string.len()
            } else {
                self.dna_range.len()
            }
        } else if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.bitpackers.columnar.len()
        } else if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.bitpackers.packed.len()
        } else {
            panic!("Parser config error: dna is ignored")
        }
    }

    #[cfg(feature = "paraseq")]
    #[inline(always)]
    fn record_index(&self) -> u64 {
        self.record_index
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    #[inline(always)]
    fn return_error_or_none(&mut self) -> Option<Event> {
        self.pending_error.take().map(Event::IoError)
    }

    #[inline(always)]
    const fn global_pos(&self) -> usize {
        64 * self.block_counter + self.pos_in_block
    }

    #[inline(always)]
    fn increment_pos(&mut self) {
        if self.pos_in_block + 1 < self.block.len {
            self.pos_in_block += 1;
        } else {
            match self.lexer.next() {
                Some(Ok(b)) => {
                    self.block = b;
                    self.block_counter += 1;
                    self.pos_in_block = 0;
                }
                Some(Err(e)) => {
                    self.pending_error = Some(e);
                    self.finished = true;
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

    #[inline(always)]
    const fn anchor(&self) -> usize {
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.header_range.start
        } else {
            self.dna_range.start
        }
    }

    #[inline(always)]
    fn make_room_record(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.anchor());
        }
    }

    #[inline(always)]
    fn make_room_record_end(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.pos_in_block + 1 >= self.block.len
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.anchor());
        }
    }

    #[inline(always)]
    fn make_room_quality(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_QUALITY)
            && flag_is_not_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.quality_range.start);
        }
    }

    #[inline(always)]
    fn make_room_quality_end(&mut self) {
        if flag_is_set(CONFIG, COMPUTE_QUALITY)
            && flag_is_not_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.pos_in_block + 1 >= self.block.len
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.quality_range.start);
        }
    }

    #[inline(always)]
    fn append_dna_range(&mut self, start: usize, len: usize) {
        self.bitpackers
            .append_range_from::<CONFIG>(&self.block, start, len);
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += len;
        }
    }

    #[inline(always)]
    fn append_dna_block_from(&mut self, block: Chunk) {
        self.bitpackers.append_block_from::<CONFIG>(&block);
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += 64;
        }
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
                        return self.return_error_or_none();
                    }
                    #[cfg(feature = "paraseq")]
                    {
                        self.record_index = self.record_counter;
                        self.record_counter += 1;
                    }
                    self.record_start = self.global_pos() - 1;
                    if flag_is_not_set(CONFIG, MERGE_RECORDS) {
                        self.clear_record();
                    }
                    if flag_is_set(CONFIG, COMPUTE_HEADER) {
                        self.header_range.start = self.global_pos();
                    }
                    while self.block.newline == 0 {
                        if flag_is_set(CONFIG, COMPUTE_HEADER)
                            && !I::RANDOM_ACCESS
                            && self.lexer.input.is_end_of_buffer()
                        {
                            self.lexer.input.make_room(self.header_range.start);
                        }
                        self.block = match self.lexer.next_in_header() {
                            Some(Ok(b)) => b,
                            Some(Err(e)) => {
                                self.finished = true;
                                self.pending_error = Some(e);
                                return self.return_error_or_none();
                            }
                            None => {
                                self.finished = true;
                                return self.return_error_or_none();
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    self.header_range.end = self.global_pos();
                    if flag_is_set(CONFIG, COMPUTE_HEADER)
                        && !I::RANDOM_ACCESS
                        && self.pos_in_block + 1 >= self.block.len
                        && self.lexer.input.is_end_of_buffer()
                    {
                        self.lexer.input.make_room(self.header_range.start);
                    }
                    self.consume_newline();
                }
                1 => {
                    // SEQUENCE
                    if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                        let mask = !0 << self.pos_in_block;
                        let mut position = (self.block.is_dna | self.block.newline) & mask;
                        while position == 0 {
                            self.make_room_record();
                            self.block = match self.lexer.next() {
                                Some(Ok(b)) => b,
                                Some(Err(e)) => {
                                    self.finished = true;
                                    self.pending_error = Some(e);
                                    return self.return_error_or_none();
                                }
                                None => {
                                    self.finished = true;
                                    return self.return_error_or_none();
                                }
                            };
                            self.block_counter += 1;
                            self.pos_in_block = 0;
                            position = self.block.is_dna | self.block.newline;
                        }
                        self.pos_in_block = position.trailing_zeros() as usize;
                        if ((1 << self.pos_in_block) & self.block.newline) != 0 {
                            self.dna_range.end = self.global_pos();
                            self.make_room_record_end();
                            self.consume_newline();
                            continue;
                        }
                    }

                    let mask = !0 << self.pos_in_block;
                    let mut position = !self.block.is_dna & mask;
                    if flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS) {
                        self.clear_chunk();
                    }
                    self.dna_range.start = self.global_pos();
                    let merge_chunks = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && flag_is_set(CONFIG, MERGE_DNA_CHUNKS);
                    if position != 0 {
                        let first_pos = self.pos_in_block;
                        self.pos_in_block = position.trailing_zeros() as usize;
                        self.append_dna_range(first_pos, self.pos_in_block - first_pos);
                    } else {
                        if self.pos_in_block > 0 {
                            self.make_room_record();
                            self.append_dna_range(self.pos_in_block, 64 - self.pos_in_block);
                            self.block = match self.lexer.next() {
                                Some(Ok(b)) => b,
                                Some(Err(e)) => {
                                    self.finished = true;
                                    self.pending_error = Some(e);
                                    return self.return_error_or_none();
                                }
                                None => {
                                    self.finished = true;
                                    return self.return_error_or_none();
                                }
                            };
                            self.block_counter += 1;
                            self.pos_in_block = 0;
                            position = !self.block.is_dna;
                        }
                        let mut block = self.block;
                        let mut dna_blocks: usize = 0;
                        while position == 0 {
                            self.make_room_record();
                            self.append_dna_block_from(block);
                            dna_blocks += 1;

                            if dna_blocks == 16 && !merge_chunks {
                                let anchor = self.anchor();
                                self.lexer.skip_long_dna(
                                    &mut self.bitpackers,
                                    &mut self.dna_len,
                                    &mut self.block_counter,
                                    anchor,
                                    flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING),
                                );
                                dna_blocks = 15;
                            }
                            block = match self.lexer.next() {
                                Some(Ok(b)) => b,
                                Some(Err(e)) => {
                                    self.block = block;
                                    self.finished = true;
                                    self.pending_error = Some(e);
                                    return self.return_error_or_none();
                                }
                                None => {
                                    self.block = block;
                                    self.finished = true;
                                    return self.return_error_or_none();
                                }
                            };
                            self.block_counter += 1;
                            position = !block.is_dna;
                        }
                        self.block = block;
                        self.pos_in_block = position.trailing_zeros() as usize;
                        self.append_dna_range(0, self.pos_in_block);
                    }
                    self.dna_range.end = self.global_pos();
                    if merge_chunks {
                        let dna_chunk = self.lexer.input.slice(self.dna_range.clone());
                        self.cur_dna_string.extend_from_slice(dna_chunk);
                    }
                    self.make_room_record_end();
                    if flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                        || ((1 << self.pos_in_block) & self.block.newline) != 0
                    {
                        self.consume_newline();
                    }
                    if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                        return Some(Event::DnaChunk(self.dna_range.clone()));
                    }
                }
                2 => {
                    // PLUS
                    while self.block.newline == 0 {
                        self.make_room_record();
                        self.block = match self.lexer.next_in_non_dna() {
                            Some(Ok(b)) => b,
                            Some(Err(e)) => {
                                self.finished = true;
                                self.pending_error = Some(e);
                                return self.return_error_or_none();
                            }
                            None => {
                                self.finished = true;
                                return self.return_error_or_none();
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    self.make_room_record_end();
                    self.consume_newline();
                }
                3 => {
                    // QUALITY
                    let qual_start = self.global_pos();
                    if flag_is_set(CONFIG, COMPUTE_QUALITY) {
                        self.quality_range.start = qual_start;
                    }
                    let target_pos = qual_start + self.dna_range.end - (self.header_range.end + 1);
                    if target_pos / 64 > self.block_counter && self.block.newline == 0 {
                        let keep_anchor = flag_is_set(
                            CONFIG,
                            COMPUTE_HEADER | COMPUTE_DNA_STRING | COMPUTE_QUALITY,
                        );
                        let anchor = if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING) {
                            self.anchor()
                        } else {
                            qual_start
                        };
                        self.lexer.skip_quality_line(
                            target_pos,
                            &mut self.block_counter,
                            &mut self.pos_in_block,
                            &mut self.block,
                            anchor,
                            keep_anchor,
                        );
                    }
                    while self.block.newline == 0 {
                        self.make_room_record();
                        self.make_room_quality();
                        self.block = match self.lexer.next_in_non_dna() {
                            Some(Ok(b)) => b,
                            Some(Err(e)) => {
                                self.pending_error = Some(e);
                                self.finished = true;
                                break;
                            }
                            None => {
                                self.finished = true;
                                break;
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                    }
                    self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                    if flag_is_set(CONFIG, COMPUTE_QUALITY) {
                        self.pos_in_block = self.pos_in_block.min(self.block.len);
                        self.quality_range.end = self.global_pos();
                    }
                    self.make_room_record_end();
                    self.make_room_quality_end();
                    self.consume_newline();
                    if flag_is_set(CONFIG, RETURN_RECORD) {
                        return Some(Event::Record(self.record_start..self.global_pos()));
                    }
                }
                _ => unreachable!(),
            }
        }
    }
}
