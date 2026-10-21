use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;
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
/// use helicase::parser::FastaParser;
/// use helicase::*;
///
/// // set the options of the parser (at compile-time)
/// const CONFIG: Config = ParserOptions::default().config();
///
/// fn main() -> std::io::Result<()> {
///     let path = "...";
///
///     // create a parser with the desired options
///     let mut parser = FastaParser::<CONFIG, _>::from_file(&path)?;
///
///     // iterate over records
///     while let Some(event) = parser.next() {
///         event?;
///         let header = parser.get_header();
///         let seq = parser.get_dna_string();
///         // ...
///     }
///     Ok(())
/// }
/// ```
pub struct FastaParser<'a, const CONFIG: Config, I: InputData<'a>> {
    lexer: Lexer<'a, CONFIG, I>,
    finished: bool,
    state: State,
    block: Chunk,
    block_counter: usize,
    pos_in_block: usize,
    header_range: Range<usize>,
    dna_range: Range<usize>,
    record_start: usize,
    contiguous_dna: bool,
    cur_dna_string: Vec<u8>,
    bitpackers: DnaBitpackers,
    dna_len: usize,
    pending_error: Option<io::Error>,
    crlf_len: usize,
    crlf_detected: bool,
    #[cfg(feature = "paraseq")]
    record_index: u64,
    #[cfg(feature = "paraseq")]
    record_counter: u64,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaParser<'a, CONFIG, I> {
    fn from_lexer(mut lexer: Lexer<'a, CONFIG, I>) -> Self {
        let mut finished: bool = false;
        let mut pending_error = None;
        let first = match lexer.next() {
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
            state: State::Header,
            block: first,
            block_counter: 0,
            pos_in_block: 0,
            header_range: 0..0,
            dna_range: 0..0,
            record_start: 0,
            contiguous_dna: true,
            cur_dna_string: Vec::new(),
            bitpackers: DnaBitpackers::new(),
            dna_len: 0,
            pending_error,
            crlf_len: 0,
            crlf_detected: false,
            #[cfg(feature = "paraseq")]
            record_index: 0,
            #[cfg(feature = "paraseq")]
            record_counter: 0,
        }
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FromInputData<'a, I>
    for FastaParser<'a, CONFIG, I>
{
    fn from_input(input: I) -> Result<Self> {
        let lexer = Lexer::from_input(input)?;
        if lexer.input.first_byte() != b'>' {
            return Err(Error::invalid_record_start("'>'", lexer.input.first_byte()));
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
        self.bitpackers.clear::<CONFIG>();
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
    fn get_dna_string(&self) -> &[u8] {
        if flag_is_not_set(CONFIG, COMPUTE_DNA_STRING) {
            panic!("Parser config error: dna_string is not enabled")
        }
        let zero_copy = self.contiguous_dna
            && !(flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS));
        if zero_copy {
            self.lexer.input.slice(self.dna_range.clone())
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
        if zero_copy {
            self.lexer.input.slice(self.dna_range.clone()).to_vec()
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
            self.bitpackers.columnar.len()
        } else if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
            self.bitpackers.packed.len()
        } else {
            unreachable!()
        }
    }

    #[cfg(feature = "paraseq")]
    #[inline(always)]
    fn record_index(&self) -> u64 {
        self.record_index
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastaParser<'a, CONFIG, I> {
    #[inline(always)]
    fn return_error_or_none(&mut self) -> Option<Result<Event>> {
        self.pending_error.take().map(Error::from).map(Err)
    }

    #[inline(always)]
    const fn global_pos(&self) -> usize {
        64 * self.block_counter + self.pos_in_block
    }

    #[inline(always)]
    const fn anchor(&self, line_start: usize) -> usize {
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.header_range.start
        } else if self.contiguous_dna {
            self.dna_range.start
        } else {
            line_start
        }
    }

    #[inline(always)]
    fn make_room_record(&mut self, line_start: usize) {
        if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.anchor(line_start));
        }
    }

    #[inline(always)]
    fn make_room_record_end(&mut self, line_start: usize) {
        if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
            && !I::RANDOM_ACCESS
            && self.pos_in_block + 1 >= self.block.len
            && self.lexer.input.is_end_of_buffer()
        {
            self.lexer.input.make_room(self.anchor(line_start));
        }
    }

    #[inline(always)]
    fn current_byte(&self) -> u8 {
        if self.pos_in_block < self.block.len {
            self.lexer.input.current_block()[self.pos_in_block & 63]
        } else {
            0
        }
    }

    #[inline(always)]
    fn skip_blank_lines(&mut self, line_start: usize) {
        while !self.finished {
            if ((1 << self.pos_in_block) & self.block.newline) != 0 {
                self.make_room_record_end(line_start);
                self.consume_newline();
            } else if self.crlf_len != 0 && self.current_byte() == b'\r' {
                self.increment_pos();
            } else {
                break;
            }
        }
    }

    #[inline(always)]
    fn increment_pos(&mut self) {
        if self.pos_in_block + 1 < self.block.len {
            self.pos_in_block += 1;
        } else {
            match self.lexer.next() {
                Some(Ok(b)) => {
                    self.block.set::<CONFIG>(b);
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
        self.bitpackers
            .append_range_from::<CONFIG>(&self.block, start, len);
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += len;
        }
    }

    #[inline(always)]
    fn append_dna_block(&mut self) {
        self.bitpackers.append_block_from::<CONFIG>(&self.block);
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += 64;
        }
    }

    #[inline(always)]
    fn emit_record_events(&mut self) -> Option<Result<Event>> {
        self.state = State::Header;
        if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
            if flag_is_set(CONFIG, RETURN_RECORD) {
                self.state = State::EmitRecord;
            }
            return Some(Ok(Event::DnaChunk(self.dna_range.clone())));
        }
        if flag_is_set(CONFIG, RETURN_RECORD) {
            return Some(Ok(Event::Record(self.record_start..self.global_pos())));
        }
        None
    }

    #[inline(always)]
    fn finish_record_eof(&mut self) -> Option<Result<Event>> {
        self.emit_record_events()
            .or_else(|| self.return_error_or_none())
    }

    #[inline(always)]
    fn advance_header(&mut self) -> Option<Option<Result<Event>>> {
        if self.finished {
            return Some(self.return_error_or_none());
        }
        #[cfg(feature = "paraseq")]
        {
            self.record_index = self.record_counter;
            self.record_counter += 1;
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
            match self.lexer.next_in_header() {
                Some(Ok(b)) => self.block.set::<CONFIG>(b),
                Some(Err(e)) => {
                    self.finished = true;
                    self.pending_error = Some(e);
                    return Some(self.return_error_or_none());
                }
                None => {
                    self.pos_in_block = self.block.len;
                    self.finished = true;
                    if flag_is_set(CONFIG, COMPUTE_HEADER) {
                        let end_pos = self.global_pos();
                        let crlf = if !self.crlf_detected
                            && end_pos > 0
                            && self.lexer.input.slice(end_pos - 1..end_pos) == b"\r"
                        {
                            1
                        } else {
                            self.crlf_len
                        };
                        self.header_range.end = end_pos.saturating_sub(crlf);
                    }
                    self.state = State::Sequence;
                    break;
                }
            }
            self.block_counter += 1;
            self.pos_in_block = 0;
            position = self.block.newline;
        }

        if !self.finished {
            self.pos_in_block = position.trailing_zeros() as usize;
            let nl_pos = self.global_pos();
            if !self.crlf_detected {
                if nl_pos > 0 && self.lexer.input.slice(nl_pos - 1..nl_pos) == b"\r" {
                    self.crlf_len = 1;
                }
                self.crlf_detected = true;
            }
            if flag_is_set(CONFIG, COMPUTE_HEADER) {
                self.header_range.end = nl_pos.saturating_sub(self.crlf_len);
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
        None
    }

    #[inline(always)]
    fn skip_non_actg(&mut self) -> Option<Option<Result<Event>>> {
        if flag_is_set(CONFIG, SPLIT_NON_ACTG)
            && ((1 << self.pos_in_block) & self.block.is_dna) == 0
        {
            let mask = !0 << self.pos_in_block;
            let mut position = (self.block.is_dna | self.block.newline) & mask;
            while position == 0 {
                self.make_room_record(self.global_pos());
                match self.lexer.next() {
                    Some(Ok(b)) => self.block.set::<CONFIG>(b),
                    Some(Err(e)) => {
                        self.finished = true;
                        self.pending_error = Some(e);
                        return Some(self.return_error_or_none());
                    }
                    None => {
                        self.pos_in_block = self.block.len;
                        self.finished = true;
                        break;
                    }
                }
                self.block_counter += 1;
                self.pos_in_block = 0;
                position = self.block.is_dna | self.block.newline;
            }
            if self.finished {
                return Some(self.finish_record_eof());
            }
            self.pos_in_block = position.trailing_zeros() as usize;
            if ((1 << self.pos_in_block) & self.block.newline) != 0 {
                self.make_room_record_end(self.global_pos());
                self.consume_newline();
                return Some(None);
            }

            if flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS) {
                self.clear_chunk();
            }
            self.dna_range.start = self.global_pos();
        }
        None
    }

    #[inline(always)]
    fn try_wrapped_stride(&mut self, next_line_start: usize) -> Option<Option<Result<Event>>> {
        let line_width = self.dna_range.len();
        let newline_len = next_line_start.saturating_sub(self.dna_range.end);
        if line_width == 0 || newline_len == 0 {
            return None;
        }
        let line_stride = line_width + newline_len;
        let expected_newline_offset = line_width + newline_len - 1;
        let min_safe_bytes = line_width.div_ceil(64) * 64;
        let stride_margin = line_stride.max(min_safe_bytes);

        let need_simd = flag_is_set(
            CONFIG,
            COMPUTE_DNA_PACKED
                | COMPUTE_DNA_COLUMNAR
                | COMPUTE_MASK_NON_ACTG
                | COMPUTE_MASK_N
                | SPLIT_NON_ACTG,
        );

        if !I::RANDOM_ACCESS {
            let rem = (self.lexer.input.buffer_offset() + self.lexer.input.buffer().len())
                .saturating_sub(next_line_start);
            if rem < stride_margin.max(65536) {
                let keep = if flag_is_set(CONFIG, COMPUTE_HEADER) {
                    self.header_range.start
                } else {
                    next_line_start
                };
                self.lexer.input.make_room(keep);
            }
        }

        let mut global_curr = next_line_start;
        let mut lines_read = 0;

        loop {
            let buf_off = self.lexer.input.buffer_offset();
            let data = self.lexer.input.buffer();
            let mut curr_pos = global_curr - buf_off;

            if lines_read == 0
                && (curr_pos + stride_margin > data.len()
                    || data[curr_pos] == b'>'
                    || !(newline_len == 1 || data[curr_pos + line_width] == b'\r')
                    || memchr::memchr(b'\n', &data[curr_pos..curr_pos + line_stride])
                        != Some(expected_newline_offset))
            {
                return None;
            }

            let estimated_dna_len = data.len().saturating_sub(curr_pos) / line_stride * line_width;
            if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                self.cur_dna_string.reserve(estimated_dna_len);
            }
            self.bitpackers.reserve::<CONFIG>(estimated_dna_len);

            let mut stopped_early = false;
            while curr_pos + stride_margin <= data.len()
                && data[curr_pos] != b'>'
                && (newline_len == 1 || data[curr_pos + line_width] == b'\r')
                && memchr::memchr(b'\n', &data[curr_pos..curr_pos + line_stride])
                    == Some(expected_newline_offset)
            {
                if need_simd
                    && !self
                        .bitpackers
                        .append_line_stride::<CONFIG>(&data[curr_pos..], line_width)
                {
                    stopped_early = true;
                    break;
                }
                if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                    self.cur_dna_string
                        .extend_from_slice(&data[curr_pos..curr_pos + line_width]);
                }
                curr_pos += line_stride;
                lines_read += 1;
            }

            global_curr = buf_off + curr_pos;
            if I::RANDOM_ACCESS || stopped_early || curr_pos + stride_margin <= data.len() {
                break;
            }
            let prev_end = buf_off + data.len();
            self.lexer.input.seek_block(global_curr / 64);
            let keep = if flag_is_set(CONFIG, COMPUTE_HEADER) {
                self.header_range.start
            } else {
                global_curr
            };
            self.lexer.input.make_room(keep);
            if self.lexer.input.buffer_offset() + self.lexer.input.buffer().len() <= prev_end {
                break;
            }
        }

        if lines_read == 0 {
            return None;
        }

        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len += lines_read * line_width;
        }
        let buf_off = self.lexer.input.buffer_offset();
        let data = self.lexer.input.buffer();
        let mut curr_pos = global_curr - buf_off;
        while (I::RANDOM_ACCESS && curr_pos < data.len() || curr_pos + 1 < data.len())
            && (data[curr_pos] == b'\n' || data[curr_pos] == b'\r')
        {
            curr_pos += 1;
        }
        let global_curr = buf_off + curr_pos;
        let data_len = data.len();
        self.dna_range = global_curr..global_curr;
        self.block_counter = global_curr / 64;
        self.pos_in_block = global_curr % 64;
        self.lexer.input.seek_block(self.block_counter);
        if I::RANDOM_ACCESS && curr_pos >= data_len {
            self.finished = true;
            return Some(self.finish_record_eof());
        }
        let block_len = 64.min((buf_off + data_len).saturating_sub(self.block_counter * 64));
        self.block
            .set::<CONFIG>(self.lexer.scan_current_block(block_len));

        if self.current_byte() == b'>' {
            return Some(self.emit_record_events());
        }

        Some(None)
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for FastaParser<'a, CONFIG, I> {
    type Item = Result<Event>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        'parser_loop: loop {
            match self.state {
                State::EmitRecord => {
                    self.state = State::Header;
                    if flag_is_set(CONFIG, RETURN_RECORD) {
                        return Some(Ok(Event::Record(self.record_start..self.global_pos())));
                    }
                }
                State::Header => {
                    if let Some(event_opt) = self.advance_header() {
                        return event_opt;
                    }
                }
                State::Sequence => {
                    // Skip blank lines
                    self.skip_blank_lines(self.global_pos());

                    if self.finished {
                        return self.finish_record_eof();
                    }

                    if self.current_byte() == b'>' {
                        self.state = State::Header;
                        if flag_is_set(CONFIG, RETURN_RECORD) {
                            return Some(Ok(Event::Record(self.record_start..self.global_pos())));
                        }
                        continue 'parser_loop;
                    }

                    if let Some(event_opt) = self.skip_non_actg() {
                        if let Some(event) = event_opt {
                            return Some(event);
                        }
                        continue 'parser_loop;
                    }

                    let merge_chunks = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                        && flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && flag_is_set(CONFIG, MERGE_DNA_CHUNKS);

                    'dna_line: loop {
                        let line_start = self.global_pos();
                        let copy_str = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
                            && (!self.contiguous_dna || merge_chunks);
                        let mask = !0 << self.pos_in_block;
                        let mut position = !self.block.is_dna & mask;
                        if position != 0 {
                            let first_pos = self.pos_in_block;
                            self.pos_in_block = position.trailing_zeros() as usize;
                            let len = if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                                (self.pos_in_block - first_pos).saturating_sub(self.crlf_len)
                            } else {
                                self.pos_in_block - first_pos
                            };
                            self.append_dna_range(first_pos, len);
                        } else {
                            if self.pos_in_block > 0 {
                                self.make_room_record(line_start);
                                self.append_dna_range(self.pos_in_block, 64 - self.pos_in_block);
                                match self.lexer.next() {
                                    Some(Ok(b)) => {
                                        self.block.set::<CONFIG>(b);
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
                            let mut dna_blocks: usize = 0;
                            while position == 0 {
                                self.make_room_record(line_start);
                                self.append_dna_block();
                                dna_blocks += 1;

                                if dna_blocks == 16 && !merge_chunks {
                                    let anchor = self.anchor(line_start);
                                    self.lexer.skip_long_dna(
                                        &mut self.bitpackers,
                                        &mut self.dna_len,
                                        &mut self.block_counter,
                                        anchor,
                                        flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING),
                                    );
                                    dna_blocks = 15;
                                }

                                match self.lexer.next() {
                                    Some(Ok(b)) => self.block.set::<CONFIG>(b),
                                    Some(Err(e)) => {
                                        self.finished = true;
                                        self.pending_error = Some(e);
                                        return self.return_error_or_none();
                                    }
                                    None => {
                                        self.pos_in_block = 64;
                                        self.finished = true;
                                        break;
                                    }
                                }
                                self.block_counter += 1;
                                position = !self.block.is_dna;
                            }
                            if !self.finished {
                                self.pos_in_block = position.trailing_zeros() as usize;
                                if self.pos_in_block > 0 {
                                    let len = if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                                        self.pos_in_block.saturating_sub(self.crlf_len)
                                    } else {
                                        self.pos_in_block
                                    };
                                    self.append_dna_range(0, len);
                                } else if self.crlf_len == 1
                                    && flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                                {
                                    self.bitpackers.pop_base::<CONFIG>();
                                    if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                                        self.dna_len -= 1;
                                    }
                                }
                            }
                        }
                        if !self.finished {
                            self.dna_range.end = self.global_pos().saturating_sub(
                                if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                                    self.crlf_len
                                } else {
                                    0
                                },
                            );
                        } else {
                            self.dna_range.end = self.global_pos();
                        }
                        if copy_str {
                            let dna_chunk = self.lexer.input.slice(line_start..self.dna_range.end);
                            self.cur_dna_string.extend_from_slice(dna_chunk);
                        }

                        if self.finished {
                            return self.finish_record_eof();
                        }

                        if self.pos_in_block >= self.block.len {
                            self.make_room_record(line_start);
                            match self.lexer.next() {
                                Some(Ok(b)) => {
                                    self.block.set::<CONFIG>(b);
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
                            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                return Some(Ok(Event::DnaChunk(self.dna_range.clone())));
                            }
                            break 'dna_line;
                        }

                        self.make_room_record_end(line_start);
                        self.consume_newline();

                        self.skip_blank_lines(line_start);

                        if self.finished {
                            return self.finish_record_eof();
                        }

                        if self.current_byte() == b'>' {
                            if let Some(event) = self.emit_record_events() {
                                return Some(event);
                            }
                            continue 'parser_loop;
                        }

                        if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                            && ((1 << self.pos_in_block) & self.block.is_dna) == 0
                        {
                            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                                return Some(Ok(Event::DnaChunk(self.dna_range.clone())));
                            }
                            break 'dna_line;
                        }

                        let should_stride = flag_is_set(
                            CONFIG,
                            COMPUTE_DNA_STRING
                                | COMPUTE_DNA_PACKED
                                | COMPUTE_DNA_COLUMNAR
                                | COMPUTE_DNA_LEN,
                        ) && self.contiguous_dna
                            && !merge_chunks;

                        if should_stride {
                            self.contiguous_dna = false;
                            if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                                let dna_chunk = self.lexer.input.slice(self.dna_range.clone());
                                self.cur_dna_string.extend_from_slice(dna_chunk);
                            }

                            if let Some(event_opt) = self.try_wrapped_stride(self.global_pos()) {
                                if let Some(event) = event_opt {
                                    return Some(event);
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
