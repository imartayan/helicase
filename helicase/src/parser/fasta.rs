use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;
use crate::simd::{extract_dna_bitmask, extract_fasta_bitmask};

use core::mem::swap;
use core::ops::Range;
use std::io;

enum State {
    Header,
    Sequence,
    EmitRecord,
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
    pending_error: Option<io::Error>,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaParser<'a, CONFIG, I> {
    fn from_lexer(mut lexer: FastaLexer<'a, CONFIG, I>) -> Self {
        let mut finished: bool = false;
        let mut pending_error = None;
        let first = match lexer.next() {
            Some(Ok(c)) => c,
            Some(Err(e)) => {
                finished = true;
                pending_error = Some(e);
                FastaChunk::default()
            }
            None => {
                finished = true;
                FastaChunk::default()
            }
        };
        Self {
            lexer,
            finished,
            state: State::Header,
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
            pending_error,
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
        self.contiguous_dna = true;
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
    fn return_error_or_none(&mut self) -> Option<Event> {
        self.pending_error.take().map(Event::IoError)
    }

    #[inline(always)]
    const fn global_pos(&self) -> usize {
        64 * self.block_counter + self.pos_in_block
    }

    #[inline(always)]
    fn anchor(&self) -> usize {
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.header_range.start
        } else {
            self.dna_range.start
        }
    }

    #[inline(always)]
    fn make_room_record(&mut self) {
        if (flag_is_set(CONFIG, COMPUTE_HEADER) || flag_is_set(CONFIG, COMPUTE_DNA_STRING))
            && !I::RANDOM_ACCESS
            && self.lexer.input.is_end_of_buffer()
            && (flag_is_set(CONFIG, COMPUTE_HEADER) || self.contiguous_dna)
        {
            self.lexer.input.make_room(self.anchor());
        }
    }

    #[inline(always)]
    fn make_room_record_end(&mut self) {
        if (flag_is_set(CONFIG, COMPUTE_HEADER) || flag_is_set(CONFIG, COMPUTE_DNA_STRING))
            && !I::RANDOM_ACCESS
            && self.pos_in_block + 1 >= self.block.len
            && self.lexer.input.is_end_of_buffer()
            && (flag_is_set(CONFIG, COMPUTE_HEADER) || self.contiguous_dna)
        {
            self.lexer.input.make_room(self.anchor());
        }
    }

    #[inline(always)]
    fn current_byte(&self) -> u8 {
        if self.pos_in_block < self.block.len {
            self.lexer.input.current_block()[self.pos_in_block]
        } else {
            0
        }
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
    }

    #[inline(always)]
    fn append_dna_range(&mut self, start: usize, len: usize) {
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.cur_dna_columnar.append(
                self.block.high_bit >> start,
                self.block.low_bit >> start,
                len,
            );
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.cur_dna_packed
                .append(self.block.two_bits >> (2 * start), 2 * len);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.cur_mask_non_actg
                .append(self.block.mask_non_actg >> start, len);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.cur_mask_n.append(self.block.mask_n >> start, len);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += len;
        }
    }

    #[inline(always)]
    fn append_dna_block_from(&mut self, block: FastaChunk) {
        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
            self.cur_dna_columnar
                .append_block(block.high_bit, block.low_bit);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.cur_dna_packed.append_block(block.two_bits);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
            self.cur_mask_non_actg.append_block(block.mask_non_actg);
        }
        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
            self.cur_mask_n.append_block(block.mask_n);
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += 64;
        }
    }

    #[inline(always)]
    fn finish_record_eof(&mut self) -> Option<Event> {
        self.state = State::Header;
        if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
            if flag_is_set(CONFIG, RETURN_RECORD) {
                self.state = State::EmitRecord;
            }
            return Some(Event::DnaChunk(self.dna_range.clone()));
        }
        if flag_is_set(CONFIG, RETURN_RECORD) {
            return Some(Event::Record(self.record_start..self.global_pos()));
        }
        self.return_error_or_none()
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastaParser<'a, CONFIG, I> {
    type Item = Event;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        'parser_loop: loop {
            match self.state {
                State::EmitRecord => {
                    self.state = State::Header;
                    if flag_is_set(CONFIG, RETURN_RECORD) {
                        return Some(Event::Record(self.record_start..self.global_pos()));
                    }
                }
                State::Header => {
                    if self.finished {
                        return self.return_error_or_none();
                    }
                    if flag_is_not_set(CONFIG, MERGE_RECORDS) {
                        self.clear_record();
                    }
                    self.record_start = self.global_pos();
                    if flag_is_set(CONFIG, COMPUTE_HEADER) {
                        self.header_range.start = self.record_start + 1;
                    }

                    // Scan for newline (current byte is '>', so bit pos_in_block of newline is 0)
                    let mask = !0 << self.pos_in_block;
                    let mut position = self.block.newline & mask;
                    while position == 0 {
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
                                self.pos_in_block = self.block.len;
                                self.finished = true;
                                if flag_is_set(CONFIG, COMPUTE_HEADER) {
                                    self.header_range.end = self.global_pos();
                                }
                                self.state = State::Sequence;
                                break;
                            }
                        };
                        self.block_counter += 1;
                        self.pos_in_block = 0;
                        position = self.block.newline;
                    }

                    if !self.finished {
                        self.pos_in_block = position.trailing_zeros() as usize;
                        if flag_is_set(CONFIG, COMPUTE_HEADER) {
                            self.header_range.end = self.global_pos();
                        }
                        if flag_is_set(CONFIG, COMPUTE_HEADER)
                            && !I::RANDOM_ACCESS
                            && self.pos_in_block + 1 >= self.block.len
                            && self.lexer.input.is_end_of_buffer()
                        {
                            self.lexer.input.make_room(self.header_range.start);
                        }
                        self.consume_newline();
                    }

                    self.contiguous_dna = true;
                    self.dna_range = self.global_pos()..self.global_pos();
                    self.state = State::Sequence;
                }
                State::Sequence => {
                    // Skip blank lines
                    while !self.finished
                        && (self.current_byte() == b'\n' || self.current_byte() == b'\r')
                    {
                        self.consume_newline();
                    }

                    if self.finished {
                        self.state = State::Header;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            return Some(Event::Record(self.record_start..self.global_pos()));
                        }
                        return self.return_error_or_none();
                    }

                    if self.current_byte() == b'>' {
                        self.state = State::Header;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            return Some(Event::Record(self.record_start..self.global_pos()));
                        }
                        continue;
                    }

                    // If SPLIT_NON_ACTG and on non-ACTG: skip non-ACTG characters until DNA or \n
                    if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && ((1 << self.pos_in_block) & self.block.is_dna) == 0
                    {
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
                                    self.pos_in_block = self.block.len;
                                    self.finished = true;
                                    break;
                                }
                            };
                            self.block_counter += 1;
                            self.pos_in_block = 0;
                            position = self.block.is_dna | self.block.newline;
                        }
                        if self.finished {
                            self.state = State::Header;
                            if flag_is_set(CONFIG, RETURN_RECORD) {
                                return Some(Event::Record(self.record_start..self.global_pos()));
                            }
                            return self.return_error_or_none();
                        }
                        self.pos_in_block = position.trailing_zeros() as usize;
                        if ((1 << self.pos_in_block) & self.block.newline) != 0 {
                            self.consume_newline();
                            continue;
                        }

                        // Next chunk started after non-ACTG bases
                        if flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS) {
                            self.clear_chunk();
                        }
                        self.dna_range.start = self.global_pos();
                    }

                    let merge_chunks = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && flag_is_set(CONFIG, MERGE_DNA_CHUNKS);

                    'dna_line: loop {
                        // Scan DNA until !is_dna
                        let line_start = self.global_pos();
                        let copy_str = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                            && (!self.contiguous_dna || merge_chunks);
                        let mask = !0 << self.pos_in_block;
                        let mut position = !self.block.is_dna & mask;
                        if position != 0 {
                            let first_pos = self.pos_in_block;
                            self.pos_in_block = position.trailing_zeros() as usize;
                            if copy_str && !I::RANDOM_ACCESS {
                                let dna_chunk =
                                    &self.lexer.input.current_block()[first_pos..self.pos_in_block];
                                self.cur_dna_string.extend_from_slice(dna_chunk);
                            }
                            self.append_dna_range(first_pos, self.pos_in_block - first_pos);
                        } else {
                            if self.pos_in_block > 0 {
                                if copy_str && !I::RANDOM_ACCESS {
                                    let dna_chunk =
                                        &self.lexer.input.current_block()[self.pos_in_block..];
                                    self.cur_dna_string.extend_from_slice(dna_chunk);
                                }
                                self.make_room_record();
                                self.append_dna_range(self.pos_in_block, 64 - self.pos_in_block);
                                match self.lexer.next() {
                                    Some(Ok(b)) => {
                                        self.block = b;
                                        self.block_counter += 1;
                                        self.pos_in_block = 0;
                                        position = !self.block.is_dna;
                                    }
                                    Some(Err(e)) => {
                                        self.finished = true;
                                        self.pending_error = Some(e);
                                        return self.return_error_or_none();
                                    }
                                    None => {
                                        self.pos_in_block = 64;
                                        self.finished = true;
                                        position = 1;
                                    }
                                }
                            }
                            let mut block = self.block;
                            let mut dna_blocks: usize = 0;
                            while position == 0 {
                                if copy_str && !I::RANDOM_ACCESS {
                                    let dna_chunk = self.lexer.input.current_block();
                                    self.cur_dna_string.extend_from_slice(dna_chunk);
                                }
                                self.make_room_record();
                                self.append_dna_block_from(block);
                                dna_blocks += 1;

                                if dna_blocks == 16
                                    && I::RANDOM_ACCESS
                                    && self.contiguous_dna
                                    && !merge_chunks
                                {
                                    let data = self.lexer.input.data();
                                    let mut curr_pos = (self.block_counter + 1) * 64;
                                    if curr_pos < data.len() {
                                        let newline_offset =
                                            memchr::memchr(b'\n', &data[curr_pos..]);
                                        let raw_len =
                                            newline_offset.unwrap_or(data.len() - curr_pos);
                                        let mut rem_len = raw_len;
                                        if newline_offset.is_some()
                                            && rem_len > 0
                                            && data[curr_pos + rem_len - 1] == b'\r'
                                        {
                                            rem_len -= 1;
                                        }

                                        let mut next_pos = curr_pos + raw_len;
                                        if newline_offset.is_some() {
                                            next_pos += 1;
                                        }
                                        while next_pos < data.len()
                                            && (data[next_pos] == b'\n' || data[next_pos] == b'\r')
                                        {
                                            next_pos += 1;
                                        }
                                        let is_complete_unwrapped_record =
                                            next_pos >= data.len() || data[next_pos] == b'>';

                                        if is_complete_unwrapped_record {
                                            let need_packed =
                                                flag_is_set(CONFIG, COMPUTE_DNA_PACKED);
                                            let need_columnar =
                                                flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR);
                                            let need_mask_non_actg =
                                                flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG);
                                            let need_mask_n = flag_is_set(CONFIG, COMPUTE_MASK_N);
                                            let need_len = flag_is_set(CONFIG, COMPUTE_DNA_LEN);
                                            let is_keep = flag_is_not_set(CONFIG, SPLIT_NON_ACTG);

                                            let any_packing = need_packed
                                                || need_columnar
                                                || need_mask_non_actg
                                                || need_mask_n;

                                            if !any_packing && is_keep {
                                                if need_len {
                                                    self.dna_len += rem_len;
                                                }
                                                self.dna_range.end = curr_pos + rem_len;
                                                self.block_counter = next_pos / 64;
                                                self.pos_in_block = next_pos % 64;
                                                self.lexer.input.seek_block(self.block_counter);
                                                if next_pos >= data.len() {
                                                    self.finished = true;
                                                    return self.finish_record_eof();
                                                }
                                                let block_len = 64.min(
                                                    data.len()
                                                        .saturating_sub(self.block_counter * 64),
                                                );
                                                self.block =
                                                    self.lexer.scan_current_block(block_len);

                                                self.state = State::Header;
                                                if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                                    if flag_is_set(CONFIG, RETURN_RECORD) {
                                                        self.state = State::EmitRecord;
                                                    }
                                                    return Some(Event::DnaChunk(
                                                        self.dna_range.clone(),
                                                    ));
                                                }
                                                if flag_is_set(CONFIG, RETURN_RECORD) {
                                                    return Some(Event::Record(
                                                        self.record_start..self.global_pos(),
                                                    ));
                                                }
                                                continue 'parser_loop;
                                            }

                                            // Align bitpacker to 32 bases (64 bits)
                                            let rem_bases = self.cur_dna_packed.len() % 32;
                                            let needed_bases =
                                                if rem_bases == 0 { 0 } else { 32 - rem_bases };

                                            let mut can_bulk = true;
                                            if needed_bases > 0 && rem_len >= needed_bases {
                                                let rem_data =
                                                    &data[curr_pos..curr_pos + needed_bases];
                                                let m = if curr_pos + 64 <= data.len() {
                                                    extract_dna_bitmask::<CONFIG>(
                                                        &data[curr_pos..curr_pos + 64],
                                                    )
                                                } else {
                                                    let mut tmp = [b'A'; 64];
                                                    tmp[..needed_bases].copy_from_slice(rem_data);
                                                    extract_dna_bitmask::<CONFIG>(&tmp)
                                                };

                                                let n_mask = (1u64 << needed_bases) - 1;
                                                if !is_keep && (m.is_dna & n_mask) != n_mask {
                                                    can_bulk = false;
                                                } else {
                                                    if need_packed {
                                                        let n_bits = 2 * needed_bases;
                                                        let two_bits_mask = (1u128 << n_bits) - 1;
                                                        self.cur_dna_packed.append(
                                                            m.two_bits & two_bits_mask,
                                                            n_bits,
                                                        );
                                                    }
                                                    if need_columnar {
                                                        self.cur_dna_columnar.append(
                                                            m.high_bit & n_mask,
                                                            m.low_bit & n_mask,
                                                            needed_bases,
                                                        );
                                                    }
                                                    if need_mask_non_actg {
                                                        self.cur_mask_non_actg.append(
                                                            m.mask_non_actg & n_mask,
                                                            needed_bases,
                                                        );
                                                    }
                                                    if need_mask_n {
                                                        self.cur_mask_n.append(
                                                            m.mask_n & n_mask,
                                                            needed_bases,
                                                        );
                                                    }
                                                    if need_len {
                                                        self.dna_len += needed_bases;
                                                    }
                                                    curr_pos += needed_bases;
                                                    rem_len -= needed_bases;
                                                }
                                            }

                                            if can_bulk {
                                                let q = rem_len / 64;
                                                let r = rem_len % 64;

                                                if need_packed {
                                                    self.cur_dna_packed.reserve_aligned_blocks(q);
                                                }
                                                if need_columnar {
                                                    self.cur_dna_columnar.reserve_aligned_blocks(q);
                                                }
                                                if need_mask_non_actg {
                                                    self.cur_mask_non_actg
                                                        .reserve_aligned_blocks(q);
                                                }
                                                if need_mask_n {
                                                    self.cur_mask_n.reserve_aligned_blocks(q);
                                                }

                                                let mut stopped_early = false;
                                                let mut processed_blocks = 0;

                                                if need_packed
                                                    && !need_columnar
                                                    && !need_mask_non_actg
                                                    && !need_mask_n
                                                    && is_keep
                                                {
                                                    let out_ptr = unsafe {
                                                        self.cur_dna_packed
                                                            .append_aligned_blocks_ptr(q)
                                                    };
                                                    unsafe {
                                                        bulk_pack_dna_blocks::<CONFIG>(
                                                            &data[curr_pos..],
                                                            q,
                                                            out_ptr,
                                                        );
                                                    }
                                                    processed_blocks = q;
                                                } else {
                                                    for k in 0..q {
                                                        let block_start = curr_pos + 64 * k;
                                                        let block_data =
                                                            &data[block_start..block_start + 64];
                                                        let m = extract_dna_bitmask::<CONFIG>(
                                                            block_data,
                                                        );

                                                        if !is_keep && m.is_dna != !0 {
                                                            stopped_early = true;
                                                            break;
                                                        }

                                                        if need_packed {
                                                            self.cur_dna_packed
                                                                .append_block_aligned_unchecked(
                                                                    m.two_bits,
                                                                );
                                                        }
                                                        if need_columnar {
                                                            self.cur_dna_columnar
                                                                .append_block_aligned_unchecked(
                                                                    m.high_bit, m.low_bit,
                                                                );
                                                        }
                                                        if need_mask_non_actg {
                                                            self.cur_mask_non_actg
                                                                .append_block_aligned_unchecked(
                                                                    m.mask_non_actg,
                                                                );
                                                        }
                                                        if need_mask_n {
                                                            self.cur_mask_n
                                                                .append_block_aligned_unchecked(
                                                                    m.mask_n,
                                                                );
                                                        }
                                                        processed_blocks += 1;
                                                    }
                                                }

                                                if stopped_early {
                                                    let processed_bytes = 64 * processed_blocks;
                                                    if need_len {
                                                        self.dna_len += processed_bytes;
                                                    }
                                                    let resume_pos = curr_pos + processed_bytes;
                                                    self.dna_range.end = resume_pos;
                                                    self.block_counter = resume_pos / 64;
                                                    self.pos_in_block = resume_pos % 64;
                                                    self.lexer.input.seek_block(self.block_counter);
                                                    let block_len =
                                                        64.min(data.len().saturating_sub(
                                                            self.block_counter * 64,
                                                        ));
                                                    self.block =
                                                        self.lexer.scan_current_block(block_len);
                                                    break;
                                                } else {
                                                    let mut remainder_ok = true;
                                                    if r > 0 {
                                                        let rem_start = curr_pos + 64 * q;
                                                        let rem_data =
                                                            &data[rem_start..curr_pos + rem_len];
                                                        let m = if rem_start + 64 <= data.len() {
                                                            extract_dna_bitmask::<CONFIG>(
                                                                &data[rem_start..rem_start + 64],
                                                            )
                                                        } else {
                                                            let mut tmp = [b'A'; 64];
                                                            tmp[..r].copy_from_slice(rem_data);
                                                            extract_dna_bitmask::<CONFIG>(&tmp)
                                                        };

                                                        let r_mask = (1u64 << r) - 1;
                                                        if !is_keep && (m.is_dna & r_mask) != r_mask
                                                        {
                                                            remainder_ok = false;
                                                        } else {
                                                            if need_packed {
                                                                let r_bits = 2 * r;
                                                                let two_bits_mask = if r_bits == 128
                                                                {
                                                                    !0
                                                                } else {
                                                                    (1u128 << r_bits) - 1
                                                                };
                                                                self.cur_dna_packed.append(
                                                                    m.two_bits & two_bits_mask,
                                                                    r_bits,
                                                                );
                                                            }
                                                            if need_columnar {
                                                                self.cur_dna_columnar.append(
                                                                    m.high_bit & r_mask,
                                                                    m.low_bit & r_mask,
                                                                    r,
                                                                );
                                                            }
                                                            if need_mask_non_actg {
                                                                self.cur_mask_non_actg.append(
                                                                    m.mask_non_actg & r_mask,
                                                                    r,
                                                                );
                                                            }
                                                            if need_mask_n {
                                                                self.cur_mask_n
                                                                    .append(m.mask_n & r_mask, r);
                                                            }
                                                        }
                                                    }

                                                    if !remainder_ok {
                                                        let processed_bytes = 64 * q;
                                                        if need_len {
                                                            self.dna_len += processed_bytes;
                                                        }
                                                        let resume_pos = curr_pos + processed_bytes;
                                                        self.dna_range.end = resume_pos;
                                                        self.block_counter = resume_pos / 64;
                                                        self.pos_in_block = resume_pos % 64;
                                                        self.lexer
                                                            .input
                                                            .seek_block(self.block_counter);
                                                        let block_len =
                                                            64.min(data.len().saturating_sub(
                                                                self.block_counter * 64,
                                                            ));
                                                        self.block = self
                                                            .lexer
                                                            .scan_current_block(block_len);
                                                        break;
                                                    } else {
                                                        if need_len {
                                                            self.dna_len += rem_len;
                                                        }
                                                        self.dna_range.end = curr_pos + rem_len;
                                                        self.block_counter = next_pos / 64;
                                                        self.pos_in_block = next_pos % 64;
                                                        self.lexer
                                                            .input
                                                            .seek_block(self.block_counter);
                                                        if next_pos >= data.len() {
                                                            self.finished = true;
                                                            return self.finish_record_eof();
                                                        }
                                                        let block_len =
                                                            64.min(data.len().saturating_sub(
                                                                self.block_counter * 64,
                                                            ));
                                                        self.block = self
                                                            .lexer
                                                            .scan_current_block(block_len);

                                                        self.state = State::Header;
                                                        if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                                            if flag_is_set(CONFIG, RETURN_RECORD) {
                                                                self.state = State::EmitRecord;
                                                            }
                                                            return Some(Event::DnaChunk(
                                                                self.dna_range.clone(),
                                                            ));
                                                        }
                                                        if flag_is_set(CONFIG, RETURN_RECORD) {
                                                            return Some(Event::Record(
                                                                self.record_start
                                                                    ..self.global_pos(),
                                                            ));
                                                        }
                                                        continue 'parser_loop;
                                                    }
                                                }
                                            }
                                        }
                                    }
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
                                        self.pos_in_block = 64;
                                        self.finished = true;
                                        break;
                                    }
                                };
                                self.block_counter += 1;
                                position = !block.is_dna;
                            }
                            self.block = block;
                            if !self.finished {
                                self.pos_in_block = position.trailing_zeros() as usize;
                                if copy_str && !I::RANDOM_ACCESS {
                                    let dna_chunk =
                                        &self.lexer.input.current_block()[..self.pos_in_block];
                                    self.cur_dna_string.extend_from_slice(dna_chunk);
                                }
                                self.append_dna_range(0, self.pos_in_block);
                            }
                        }
                        self.dna_range.end = self.global_pos();
                        if copy_str && I::RANDOM_ACCESS {
                            let dna_chunk =
                                &self.lexer.input.data()[line_start..self.dna_range.end];
                            self.cur_dna_string.extend_from_slice(dna_chunk);
                        }

                        if self.finished {
                            return self.finish_record_eof();
                        }

                        if self.pos_in_block >= self.block.len {
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
                                None => {
                                    self.finished = true;
                                }
                            }
                            if self.finished {
                                return self.finish_record_eof();
                            }
                        }

                        if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                            && ((1 << self.pos_in_block) & self.block.newline) == 0
                        {
                            // Stopped at non-ACTG within line
                            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                return Some(Event::DnaChunk(self.dna_range.clone()));
                            }
                            break 'dna_line;
                        }

                        // Stopped at newline
                        self.make_room_record_end();
                        self.consume_newline();

                        // Skip any subsequent blank lines
                        while !self.finished
                            && (self.current_byte() == b'\n' || self.current_byte() == b'\r')
                        {
                            self.consume_newline();
                        }

                        if self.finished {
                            return self.finish_record_eof();
                        }

                        if self.current_byte() == b'>' {
                            // Next record starts!
                            self.state = State::Header;
                            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                if flag_is_set(CONFIG, RETURN_RECORD) {
                                    self.state = State::EmitRecord;
                                }
                                return Some(Event::DnaChunk(self.dna_range.clone()));
                            }
                            if flag_is_set(CONFIG, RETURN_RECORD) {
                                return Some(Event::Record(self.record_start..self.global_pos()));
                            }
                            break 'dna_line;
                        }

                        // If non-ACTG base starts next line in SPLIT_NON_ACTG mode, previous chunk ends
                        if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                            && ((1 << self.pos_in_block) & self.block.is_dna) == 0
                        {
                            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                return Some(Event::DnaChunk(self.dna_range.clone()));
                            }
                            break 'dna_line;
                        }

                        // Wrapped line continuation!
                        let should_stride = (flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                            || flag_is_set(CONFIG, COMPUTE_DNA_PACKED)
                            || flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR)
                            || flag_is_set(CONFIG, COMPUTE_DNA_LEN))
                            && self.contiguous_dna
                            && !merge_chunks;

                        if should_stride {
                            self.contiguous_dna = false;
                            if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                if I::RANDOM_ACCESS {
                                    let dna_chunk =
                                        &self.lexer.input.data()[self.dna_range.clone()];
                                    self.cur_dna_string.extend_from_slice(dna_chunk);
                                } else {
                                    let buf_off = self.lexer.input.buffer_offset();
                                    let start = self.dna_range.start - buf_off;
                                    let end = self.dna_range.end - buf_off;
                                    self.cur_dna_string
                                        .extend_from_slice(&self.lexer.input.buffer()[start..end]);
                                }
                            }

                            if I::RANDOM_ACCESS {
                                let data = self.lexer.input.data();
                                let line_width = self.dna_range.len();
                                let next_line_start = self.global_pos();
                                let newline_len =
                                    next_line_start.saturating_sub(self.dna_range.end);
                                if line_width > 0 && newline_len > 0 {
                                    let line_stride = line_width + newline_len;
                                    let expected_newline_offset = line_width + newline_len - 1;
                                    let q = line_width / 64;
                                    let r = line_width % 64;
                                    let min_safe_bytes = 64 * q + if r > 0 { 64 } else { 0 };
                                    let stride_margin = line_stride.max(min_safe_bytes);

                                    if next_line_start + stride_margin <= data.len()
                                        && data[next_line_start] != b'>'
                                        && (newline_len == 1
                                            || data[next_line_start + line_width] == b'\r')
                                        && memchr::memchr(
                                            b'\n',
                                            &data[next_line_start..next_line_start + line_stride],
                                        ) == Some(expected_newline_offset)
                                    {
                                        let mut curr_pos = next_line_start;
                                        let estimated_dna_len =
                                            (data.len() - curr_pos) / line_stride * line_width;
                                        if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                            self.cur_dna_string.reserve(estimated_dna_len);
                                        }
                                        if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                                            self.cur_dna_packed.reserve(estimated_dna_len);
                                        }
                                        if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                                            self.cur_dna_columnar.reserve(estimated_dna_len);
                                        }
                                        if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                                            self.cur_mask_non_actg.reserve(estimated_dna_len);
                                        }
                                        if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                                            self.cur_mask_n.reserve(estimated_dna_len);
                                        }

                                        let need_simd = flag_is_set(
                                            CONFIG,
                                            COMPUTE_DNA_PACKED
                                                | COMPUTE_DNA_COLUMNAR
                                                | COMPUTE_MASK_NON_ACTG
                                                | COMPUTE_MASK_N
                                                | SPLIT_NON_ACTG,
                                        );

                                        let mut lines_read = 0;
                                        while curr_pos + stride_margin <= data.len()
                                            && data[curr_pos] != b'>'
                                            && (newline_len == 1
                                                || data[curr_pos + line_width] == b'\r')
                                            && memchr::memchr(
                                                b'\n',
                                                &data[curr_pos..curr_pos + line_stride],
                                            ) == Some(expected_newline_offset)
                                        {
                                            if !need_simd {
                                                if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                                    self.cur_dna_string.extend_from_slice(
                                                        &data[curr_pos..curr_pos + line_width],
                                                    );
                                                }
                                            } else if q == 1 {
                                                let m0 = extract_fasta_bitmask::<CONFIG>(
                                                    &data[curr_pos..curr_pos + 64],
                                                );
                                                let m1 = extract_fasta_bitmask::<CONFIG>(
                                                    &data[curr_pos + 64..curr_pos + 128],
                                                );
                                                let r_mask = (1u64 << r) - 1;
                                                if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                                                    && (m0.is_dna != !0
                                                        || (m1.is_dna & r_mask) != r_mask)
                                                {
                                                    break;
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                                                    self.cur_dna_packed.append_block(m0.two_bits);
                                                    let r_bits = 2 * r;
                                                    let two_bits_mask = if r_bits == 128 {
                                                        !0
                                                    } else {
                                                        (1u128 << r_bits) - 1
                                                    };
                                                    self.cur_dna_packed.append(
                                                        m1.two_bits & two_bits_mask,
                                                        r_bits,
                                                    );
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                                                    self.cur_dna_columnar
                                                        .append_block(m0.high_bit, m0.low_bit);
                                                    self.cur_dna_columnar.append(
                                                        m1.high_bit & r_mask,
                                                        m1.low_bit & r_mask,
                                                        r,
                                                    );
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                                                    self.cur_mask_non_actg
                                                        .append_block(m0.mask_non_actg);
                                                    self.cur_mask_non_actg
                                                        .append(m1.mask_non_actg & r_mask, r);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                                                    self.cur_mask_n.append_block(m0.mask_n);
                                                    self.cur_mask_n.append(m1.mask_n & r_mask, r);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                                    self.cur_dna_string.extend_from_slice(
                                                        &data[curr_pos..curr_pos + line_width],
                                                    );
                                                }
                                            } else if q == 0 && r > 0 {
                                                let m = extract_fasta_bitmask::<CONFIG>(
                                                    &data[curr_pos..curr_pos + 64],
                                                );
                                                let r_mask = (1u64 << r) - 1;
                                                if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                                                    && (m.is_dna & r_mask) != r_mask
                                                {
                                                    break;
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                                                    let r_bits = 2 * r;
                                                    let two_bits_mask = if r_bits == 128 {
                                                        !0
                                                    } else {
                                                        (1u128 << r_bits) - 1
                                                    };
                                                    self.cur_dna_packed
                                                        .append(m.two_bits & two_bits_mask, r_bits);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                                                    self.cur_dna_columnar
                                                        .append(m.high_bit, m.low_bit, r);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                                                    self.cur_mask_non_actg
                                                        .append(m.mask_non_actg, r);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                                                    self.cur_mask_n.append(m.mask_n, r);
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                                    self.cur_dna_string.extend_from_slice(
                                                        &data[curr_pos..curr_pos + line_width],
                                                    );
                                                }
                                            } else {
                                                let mut has_non_actg = false;
                                                if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                                                    for k in 0..q {
                                                        let m = extract_fasta_bitmask::<CONFIG>(
                                                            &data[curr_pos + 64 * k
                                                                ..curr_pos + 64 * (k + 1)],
                                                        );
                                                        if m.is_dna != !0 {
                                                            has_non_actg = true;
                                                            break;
                                                        }
                                                    }
                                                    if !has_non_actg && r > 0 {
                                                        let m = extract_fasta_bitmask::<CONFIG>(
                                                            &data[curr_pos + 64 * q
                                                                ..curr_pos + 64 * q + 64],
                                                        );
                                                        let r_mask = (1u64 << r) - 1;
                                                        if (m.is_dna & r_mask) != r_mask {
                                                            has_non_actg = true;
                                                        }
                                                    }
                                                    if has_non_actg {
                                                        break;
                                                    }
                                                }
                                                for k in 0..q {
                                                    let m = extract_fasta_bitmask::<CONFIG>(
                                                        &data[curr_pos + 64 * k
                                                            ..curr_pos + 64 * (k + 1)],
                                                    );
                                                    if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                                                        self.cur_dna_packed
                                                            .append_block(m.two_bits);
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                                                        self.cur_dna_columnar
                                                            .append_block(m.high_bit, m.low_bit);
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                                                        self.cur_mask_non_actg
                                                            .append_block(m.mask_non_actg);
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                                                        self.cur_mask_n.append_block(m.mask_n);
                                                    }
                                                }
                                                if r > 0 {
                                                    let m = extract_fasta_bitmask::<CONFIG>(
                                                        &data[curr_pos + 64 * q
                                                            ..curr_pos + 64 * q + 64],
                                                    );
                                                    if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                                                        let r_bits = 2 * r;
                                                        let two_bits_mask = if r_bits == 128 {
                                                            !0
                                                        } else {
                                                            (1u128 << r_bits) - 1
                                                        };
                                                        self.cur_dna_packed.append(
                                                            m.two_bits & two_bits_mask,
                                                            r_bits,
                                                        );
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                                                        self.cur_dna_columnar
                                                            .append(m.high_bit, m.low_bit, r);
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                                                        self.cur_mask_non_actg
                                                            .append(m.mask_non_actg, r);
                                                    }
                                                    if flag_is_set(CONFIG, COMPUTE_MASK_N) {
                                                        self.cur_mask_n.append(m.mask_n, r);
                                                    }
                                                }
                                                if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                                    self.cur_dna_string.extend_from_slice(
                                                        &data[curr_pos..curr_pos + line_width],
                                                    );
                                                }
                                            }
                                            curr_pos += line_stride;
                                            lines_read += 1;
                                        }

                                        if lines_read > 0 {
                                            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                                                self.dna_len += lines_read * line_width;
                                            }
                                            while curr_pos < data.len()
                                                && (data[curr_pos] == b'\n'
                                                    || data[curr_pos] == b'\r')
                                            {
                                                curr_pos += 1;
                                            }
                                            self.dna_range.start = curr_pos;
                                            self.dna_range.end = curr_pos;
                                            self.block_counter = curr_pos / 64;
                                            self.pos_in_block = curr_pos % 64;
                                            self.lexer.input.seek_block(self.block_counter);
                                            if curr_pos >= data.len() {
                                                self.finished = true;
                                                return self.finish_record_eof();
                                            }
                                            let block_len = 64.min(
                                                data.len().saturating_sub(self.block_counter * 64),
                                            );
                                            self.block = self.lexer.scan_current_block(block_len);

                                            if self.current_byte() == b'>' {
                                                self.state = State::Header;
                                                if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                                    if flag_is_set(CONFIG, RETURN_RECORD) {
                                                        self.state = State::EmitRecord;
                                                    }
                                                    return Some(Event::DnaChunk(
                                                        self.dna_range.clone(),
                                                    ));
                                                }
                                                if flag_is_set(CONFIG, RETURN_RECORD) {
                                                    return Some(Event::Record(
                                                        self.record_start..self.global_pos(),
                                                    ));
                                                }
                                                break 'dna_line;
                                            }
                                        }
                                    }
                                }
                            }
                        }
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
        .config();
    const CONFIG_STRING_ACTG_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .skip_non_actg()
        .config();
    const CONFIG_COLUMNAR: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .config();
    const CONFIG_COLUMNAR_MERGE: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .skip_non_actg()
        .config();
    const CONFIG_PACKED: Config = ParserOptions::default()
        .ignore_headers()
        .dna_packed()
        .config();
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
