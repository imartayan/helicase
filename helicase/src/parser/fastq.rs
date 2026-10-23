use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;
use crate::lexer::*;
use core::mem::replace;
use core::ops::Range;
use std::io;

/// A parser for the [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).
/// Configured via [`ParserOptions`].
///
/// # Example
///
/// ```rust,no_run
/// use helicase::input::*;
/// use helicase::parser::FastqParser;
/// use helicase::*;
///
/// // set the options of the parser (at compile-time)
/// const CONFIG: Config = ParserOptions::default().config();
///
/// fn main() -> std::io::Result<()> {
///     let path = "...";
///
///     // create a parser with the desired options
///     let mut parser = FastqParser::<CONFIG, _>::from_file(&path)?;
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
    crlf_len: usize,
    crlf_detected: bool,
    seq_start: usize,
    line_index: LineIndex,
    #[cfg(feature = "paraseq")]
    record_index: u64,
    #[cfg(feature = "paraseq")]
    record_counter: u64,
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    // Copy the string instead of slicing the input when chunks or records are merged.
    const MERGE_DNA_STRING: bool = flag_is_set(CONFIG, COMPUTE_DNA_STRING)
        && (flag_is_set(CONFIG, SPLIT_NON_ACTG) && flag_is_set(CONFIG, MERGE_DNA_CHUNKS)
            || flag_is_set(CONFIG, MERGE_RECORDS));

    const INDEXED: bool = flag_is_set(CONFIG, RETURN_RECORD)
        && flag_is_not_set(CONFIG, RETURN_DNA_CHUNK | MERGE_RECORDS)
        && (flag_is_not_set(CONFIG, SPLIT_NON_ACTG) || flag_is_set(CONFIG, REJECT_NON_ACTG));

    fn from_lexer(mut lexer: Lexer<'a, CONFIG, I>) -> Self {
        let first_block = lexer.next_block::<SCAN_MODE_DNA>().transpose();
        let finished = !matches!(first_block, Ok(Some(_)));
        let (block, pending_error) = match first_block {
            Ok(block) => (block.unwrap_or_default(), None),
            Err(e) => (Chunk::default(), Some(e)),
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
            crlf_len: 0,
            crlf_detected: false,
            seq_start: 0,
            line_index: LineIndex::default(),
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
    fn from_input(input: I) -> Result<Self> {
        let lexer = Lexer::from_input(input)?;
        if !lexer.input.buffer().is_empty() && lexer.input.first_byte() != b'@' {
            return Err(Error::invalid_record_start("'@'", lexer.input.first_byte()));
        }
        Ok(Self::from_lexer(lexer))
    }
}

impl<'a, const CONFIG: Config, I: InputData<'a>> crate::sealed::Sealed
    for FastqParser<'a, CONFIG, I>
{
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
        self.bitpackers.clear::<CONFIG>();
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len = 0;
        }
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
        if Self::MERGE_DNA_STRING {
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
        if Self::MERGE_DNA_STRING {
            let capacity = self.cur_dna_string.capacity();
            replace(&mut self.cur_dna_string, Vec::with_capacity(capacity))
        } else {
            self.lexer.input.slice(self.dna_range.clone()).to_vec()
        }
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len
        } else if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
            if Self::MERGE_DNA_STRING {
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

    shared_getters!();
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    shared_block_helpers!();

    // Skips blank lines before a record, which must then start with `@`.
    // Other invalid record starts are not checked.
    #[cold]
    #[inline(never)]
    fn skip_blank_lines(&mut self) -> Option<Option<Result<Event>>> {
        loop {
            match self.current_byte() {
                b'@' => return None,
                b'\n' => self.consume_newline(),
                b'\r' => self.increment_pos(),
                byte => {
                    self.finished = true;
                    return Some(Some(Err(Error::invalid_record_start("'@'", byte))));
                }
            }
            if self.finished {
                return Some(self.return_error_or_none());
            }
        }
    }

    // Going back to the header ends the iteration on the next call.
    #[cold]
    fn truncated_or_error(&mut self) -> Option<Result<Event>> {
        self.line_count = 0;
        let error = self.pending_error.take();
        Some(Err(error.map_or_else(Error::truncated_record, Error::from)))
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
}

impl<'a, const CONFIG: Config, I: InputData<'a>> FastqParser<'a, CONFIG, I> {
    // Parses a record from 4 indexed line feeds, without touching the block state.
    // Anything unusual (long or malformed record, end of input) goes back to `next_event`.
    #[inline(always)]
    fn next_indexed(&mut self) -> Option<Result<Event>> {
        let (data, offset) = (self.lexer.input.buffer(), self.lexer.input.buffer_offset());
        let [header_end, seq_end, plus_end, qual_end] = self.line_index.peek::<4>(data, offset)?;
        // The block parser jumps over quality lines using the sequence length,
        // and skips blank lines.
        if qual_end - plus_end != seq_end - header_end
            || header_end - self.line_index.next_start <= self.crlf_len
        {
            return None;
        }
        let start = self.line_index.next_start;
        self.line_index.consume::<4>(qual_end + 1);
        #[cfg(feature = "paraseq")]
        {
            self.record_index = self.record_counter;
            self.record_counter += 1;
        }
        self.record_start = start;
        self.clear_record();
        let crlf_len = self.crlf_len;
        if flag_is_set(CONFIG, COMPUTE_HEADER) {
            self.header_range = start + 1..header_end.saturating_sub(crlf_len);
        }
        let seq_len = (seq_end - header_end - 1).saturating_sub(crlf_len);
        self.dna_range = header_end + 1..seq_end.saturating_sub(crlf_len);
        let data = self.lexer.input.buffer();
        let packed = self
            .bitpackers
            .pack_contiguous::<CONFIG>(&data[header_end + 1 - offset..], seq_len);
        if flag_is_set(CONFIG, REJECT_NON_ACTG) && packed < seq_len {
            return self.reject_first_non_actg(self.dna_range.clone());
        }
        if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
            self.dna_len = seq_len;
        }
        if flag_is_set(CONFIG, COMPUTE_QUALITY) {
            self.quality_range = plus_end + 1..qual_end.saturating_sub(crlf_len);
        }
        Some(Ok(Event::Record(start..qual_end + 1)))
    }

    #[inline(always)]
    fn next_event(&mut self) -> Option<Result<Event>> {
        if flag_is_set(CONFIG, REJECT_NON_ACTG) && self.finished {
            return self.return_error_or_none();
        }
        if Self::INDEXED && self.line_index.active {
            if let Some(event) = self.next_indexed() {
                return Some(event);
            }
            self.leave_indexed();
        }
        loop {
            if !flag_is_set(CONFIG, RETURN_DNA_CHUNK) || self.line_count == 0 {
                // HEADER
                // A line feed at the record start, after the `\r` with CRLF, is a blank line.
                if !self.finished
                    && (self.block.newline >> self.pos_in_block) >> self.crlf_len & 1 != 0
                    && let Some(event) = self.skip_blank_lines()
                {
                    return event;
                }
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
                    if !self.next_block::<SCAN_MODE_HEADER>() {
                        return self.truncated_or_error();
                    }
                    self.pos_in_block = 0;
                }
                self.pos_in_block = self.block.newline.trailing_zeros() as usize;
                let nl_pos = self.global_pos();
                if !self.crlf_detected {
                    if nl_pos > 0 && self.lexer.input.slice(nl_pos - 1..nl_pos) == b"\r" {
                        self.crlf_len = 1;
                    }
                    self.crlf_detected = true;
                }
                self.header_range.end = nl_pos.saturating_sub(self.crlf_len);
                if flag_is_set(CONFIG, COMPUTE_HEADER)
                    && !I::RANDOM_ACCESS
                    && self.pos_in_block + 1 >= self.block.len
                    && self.lexer.input.is_end_of_buffer()
                {
                    self.lexer.input.make_room(self.header_range.start);
                }
                self.consume_newline();
                self.seq_start = self.global_pos();
                if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                    self.line_count = 1;
                }
            }

            if !flag_is_set(CONFIG, RETURN_DNA_CHUNK) || self.line_count == 1 {
                // SEQUENCE
                'seq_line: loop {
                    if flag_is_set(CONFIG, SPLIT_NON_ACTG) {
                        let skip_start = self.global_pos();
                        let first_byte = if flag_is_set(CONFIG, REJECT_NON_ACTG) {
                            self.current_byte()
                        } else {
                            0
                        };
                        let mask = !0 << self.pos_in_block;
                        let mut position = (self.block.is_dna | self.block.newline) & mask;
                        while position == 0 {
                            self.make_room_record();
                            if !self.next_block::<SCAN_MODE_DNA>() {
                                return self.truncated_or_error();
                            }
                            self.pos_in_block = 0;
                            position = self.block.is_dna | self.block.newline;
                        }
                        self.pos_in_block = position.trailing_zeros() as usize;
                        let at_newline = ((1 << self.pos_in_block) & self.block.newline) != 0;
                        if flag_is_set(CONFIG, REJECT_NON_ACTG) {
                            // Only the `\r` of a CRLF line ending can be skipped.
                            let skipped = self.global_pos() - skip_start;
                            let crlf = if at_newline { self.crlf_len } else { 0 };
                            if skipped > 0 && (first_byte != b'\r' || skipped > crlf) {
                                return self.reject_base(first_byte, skip_start);
                            }
                        }
                        if at_newline {
                            self.dna_range.end = self.global_pos().saturating_sub(self.crlf_len);
                            // Without any ACTG run, the start is still in the previous record.
                            if self.dna_range.start < self.seq_start {
                                self.dna_range.start = self.dna_range.end;
                            }
                            self.make_room_record_end();
                            self.consume_newline();
                            break 'seq_line;
                        }
                    }

                    let mask = !0 << self.pos_in_block;
                    let mut position = !self.block.is_dna & mask;
                    if flag_is_set(CONFIG, SPLIT_NON_ACTG)
                        && flag_is_not_set(CONFIG, MERGE_DNA_CHUNKS)
                    {
                        self.clear_chunk();
                    }
                    self.dna_range.start = self.global_pos();
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
                            self.make_room_record();
                            self.append_dna_range(self.pos_in_block, 64 - self.pos_in_block);
                            if !self.next_block::<SCAN_MODE_DNA>() {
                                return self.truncated_or_error();
                            }
                            self.pos_in_block = 0;
                            position = !self.block.is_dna;
                        }
                        let mut dna_blocks: usize = 0;
                        while position == 0 {
                            self.make_room_record();
                            self.append_dna_block();
                            dna_blocks += 1;

                            if dna_blocks == 16 && !Self::MERGE_DNA_STRING {
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
                            if !self.next_block::<SCAN_MODE_DNA>() {
                                return self.truncated_or_error();
                            }
                            position = !self.block.is_dna;
                        }
                        self.pos_in_block = position.trailing_zeros() as usize;
                        if self.pos_in_block > 0 {
                            let len = if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                                self.pos_in_block.saturating_sub(self.crlf_len)
                            } else {
                                self.pos_in_block
                            };
                            self.append_dna_range(0, len);
                        } else if self.crlf_len == 1 && flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                            self.bitpackers.pop_base::<CONFIG>();
                            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                                self.dna_len -= 1;
                            }
                        }
                    }
                    self.dna_range.end = self.global_pos().saturating_sub(
                        if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                            self.crlf_len
                        } else {
                            0
                        },
                    );
                    if Self::MERGE_DNA_STRING {
                        let dna_chunk = self.lexer.input.slice(self.dna_range.clone());
                        self.cur_dna_string.extend_from_slice(dna_chunk);
                    }
                    self.make_room_record_end();
                    let at_newline = flag_is_not_set(CONFIG, SPLIT_NON_ACTG)
                        || ((1 << self.pos_in_block) & self.block.newline) != 0;
                    if at_newline {
                        self.consume_newline();
                    }
                    if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                        if at_newline {
                            self.line_count = 2;
                        }
                        return Some(Ok(Event::DnaChunk(self.dna_range.clone())));
                    }
                    if at_newline {
                        break 'seq_line;
                    }
                }
            }

            // PLUS
            while self.block.newline == 0 {
                self.make_room_record();
                if !self.next_block::<SCAN_MODE_QUALITY>() {
                    return self.truncated_or_error();
                }
                self.pos_in_block = 0;
            }
            self.pos_in_block = self.block.newline.trailing_zeros() as usize;
            self.make_room_record_end();
            self.consume_newline();

            // QUALITY
            let qual_start = self.global_pos();
            if flag_is_set(CONFIG, COMPUTE_QUALITY) {
                self.quality_range.start = qual_start;
            }
            let seq_len = self.dna_range.end.saturating_sub(self.seq_start);
            let target_pos = qual_start + seq_len + self.crlf_len;
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
            let mut qual_hit_eof = false;
            while self.block.newline == 0 {
                self.make_room_record();
                self.make_room_quality();
                if !self.next_block::<SCAN_MODE_QUALITY>() {
                    qual_hit_eof = self.pending_error.is_none();
                    break;
                }
                self.pos_in_block = 0;
            }
            self.pos_in_block = self.block.newline.trailing_zeros() as usize;
            // Without a final line feed, the quality line ends with the input.
            let crlf_len = if qual_hit_eof {
                self.pos_in_block = self.pos_in_block.min(self.block.len);
                0
            } else {
                self.crlf_len
            };
            if flag_is_set(CONFIG, COMPUTE_QUALITY) {
                self.quality_range.end = self.global_pos().saturating_sub(crlf_len);
            }
            if qual_hit_eof && self.global_pos() < qual_start + seq_len {
                return self.truncated_or_error();
            }
            self.make_room_record_end();
            self.make_room_quality_end();
            self.consume_newline();
            if flag_is_set(CONFIG, RETURN_DNA_CHUNK) {
                self.line_count = 0;
            }
            if flag_is_set(CONFIG, RETURN_RECORD) {
                if Self::INDEXED && !self.finished {
                    let (data, offset) =
                        (self.lexer.input.buffer(), self.lexer.input.buffer_offset());
                    self.line_index.record_parsed(
                        data,
                        offset,
                        self.record_start,
                        self.global_pos(),
                    );
                }
                return Some(Ok(Event::Record(self.record_start..self.global_pos())));
            }
        }
    }
}

shared_iterator!(FastqParser);
