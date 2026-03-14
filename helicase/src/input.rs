//! Input formats and helpers.

use core::marker::PhantomData;
use deko::read::AnyDecoder;
use memchr::memchr;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, Read, Stdin, stdin};
use std::path::Path;

const BUFFER_DEFAULT_SIZE: usize = 1 << 17;
const BUFFER_DOUBLE_UNTIL: usize = 1 << 24;

pub trait InputData<'a>: Iterator<Item = &'a [u8]> {
    const RANDOM_ACCESS: bool;

    /// Get a reference to the complete slice of data.
    ///
    /// This is not available for reader-based implementations.
    #[inline(always)]
    fn data(&self) -> &[u8] {
        assert!(Self::RANDOM_ACCESS);
        unimplemented!()
    }

    /// Get a reference to the current block of up to 64 bytes.
    ///
    /// If the length is smaller than 64, the following bytes are guaranteed to be zeros.
    fn current_block(&self) -> &[u8];

    /// Get a reference to the internal buffer.
    ///
    /// This is only relevant for reader-based implementations.
    fn buffer(&self) -> &[u8];

    /// Returns the offset of the buffer.
    ///
    /// This is only relevant for reader-based implementations.
    fn buffer_offset(&self) -> usize {
        0
    }

    /// Returns `true` if we are reaching the end of the buffer.
    ///
    /// This is only relevant for reader-based implementations.
    fn is_end_of_buffer(&self) -> bool;

    /// Shift buffered data starting at `keep_from` (absolute file offset) to
    /// the beginning of the buffer, then refill the rest from the reader.
    ///
    /// After this call `buffer_offset() == keep_from` and `buffer()[0..]`
    /// holds everything from `keep_from` onwards.
    ///
    /// This is only relevant for reader-based implementations.
    #[inline(always)]
    fn make_room(&mut self, _keep_from: usize) {}

    /// Grow buffer and load `additional` new bytes.
    ///
    /// This is only relevant for reader-based implementations.
    #[inline(always)]
    fn grow_buffer(&mut self, _additional: usize) {}

    /// Scan the remaining buffer for the next `\n` byte without running the
    /// SIMD lexer.  Advances the internal buffer position past the block that
    /// contains the newline and returns `(delta_blocks, pos_in_block,
    /// block_len)`:
    ///
    /// - `delta_blocks`: number of full 64-byte blocks skipped *before* the
    ///   newline block (i.e. `block_counter` should increase by `delta_blocks + 1`).
    /// - `pos_in_block`: byte offset of `\n` within the newline-containing block.
    /// - `block_len`: number of valid bytes in the newline-containing block
    ///   (`64` for a full block, less for the final partial block).
    ///
    /// Returns `None` if no `\n` is found in the remaining buffered data
    /// (caller should fall back to the block-by-block lexer loop).
    ///
    /// The caller is responsible for re-running the SIMD scan on `current_block()`
    /// to obtain fresh `is_dna`, `newline`, and other bitmasks after this call.
    ///
    /// This is only relevant for reader-based implementations; the default
    /// returns `None` (random-access inputs use the existing block loop).
    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        None
    }

    /// Returns the first byte of the (uncompressed when possible) input.
    fn first_byte(&self) -> u8;

    /// Returns the type of compression format detected.
    ///
    /// This is only available for reader-based implementations.
    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        Ok(None)
    }

    /// Returns `true` if compression has been detected.
    #[inline(always)]
    fn is_compressed(&mut self) -> io::Result<bool> {
        Ok(self.compression_format()?.is_some())
    }
}

pub trait FromInputData<'a, I: InputData<'a>>: Sized {
    /// Build the struct from a type implementing [`InputData`].
    fn from_input(input: I) -> io::Result<Self>;
}

/// Slice input.
/// It supports parallel processing, but not transparent decompression.
pub struct SliceInput<'a> {
    data: &'a [u8],
    pos: usize,
    first_byte: u8,
    last_block: [u8; 64],
}

impl<'a> SliceInput<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        assert!(!data.is_empty());
        let mut last_block = [0; 64];
        last_block[..data.len() % 64].copy_from_slice(&data[(data.len() / 64) * 64..]);
        Self {
            data,
            pos: 0,
            first_byte: data[0],
            last_block,
        }
    }
}

impl<'a> Iterator for SliceInput<'a> {
    type Item = &'a [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos;
        self.pos += 64;
        if pos + 64 <= self.data.len() {
            unsafe { Some(std::slice::from_raw_parts(self.data.as_ptr().add(pos), 64)) }
        } else if pos < self.data.len() {
            unsafe {
                Some(std::slice::from_raw_parts(
                    self.last_block.as_ptr(),
                    self.data.len() % 64,
                ))
            }
        } else {
            None
        }
    }
}

impl<'a> InputData<'a> for SliceInput<'a> {
    const RANDOM_ACCESS: bool = true;

    #[inline(always)]
    fn data(&self) -> &[u8] {
        self.data
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        if 64 <= self.pos && self.pos <= self.data.len() {
            unsafe { std::slice::from_raw_parts(self.data.as_ptr().add(self.pos - 64), 64) }
        } else {
            unsafe { std::slice::from_raw_parts(self.last_block.as_ptr(), self.data.len() % 64) }
        }
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.data
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.pos >= self.data.len()
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        if self.pos >= self.data.len() {
            return None;
        }
        let p = memchr(b'\n', &self.data[self.pos..])?;
        let delta_blocks = p / 64;
        let pos_in_block = p % 64;
        let block_start = self.pos + delta_blocks * 64;
        let block_end = (block_start + 64).min(self.data.len());
        let block_len = block_end - block_start;
        self.pos = block_start + 64;
        Some((delta_blocks, pos_in_block, block_len))
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.first_byte
    }
}

pub trait FromSlice<'a>: FromInputData<'a, SliceInput<'a>> {
    /// Build the struct from a slice.
    /// It supports parallel processing, but not transparent decompression.
    #[inline(always)]
    fn from_slice(data: &'a [u8]) -> io::Result<Self> {
        Self::from_input(SliceInput::new(data))
    }
}

impl<'a, F: FromInputData<'a, SliceInput<'a>>> FromSlice<'a> for F {}

/// Memory mapped file.
/// It supports parallel processing, but not transparent decompression.
pub struct MmapInput<'a> {
    slice: SliceInput<'a>,
    _mmap: Mmap,
}

impl<'a> MmapInput<'a> {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        // Unsafe: mmap are intrisically unsafe.
        // Here we add on the top of that a slice of [u8] that will live as long as file is not dropped.
        // To enforce this, we keep the file in the struct.
        let file = File::open(path)?;
        let _mmap = unsafe { Mmap::map(&file)? };
        let data = unsafe { std::slice::from_raw_parts(_mmap.as_ptr(), _mmap.len()) };
        Ok(Self {
            slice: SliceInput::new(data),
            _mmap,
        })
    }
}

impl<'a> Iterator for MmapInput<'a> {
    type Item = &'a [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.slice.next()
    }
}

impl<'a> InputData<'a> for MmapInput<'a> {
    const RANDOM_ACCESS: bool = true;

    #[inline(always)]
    fn data(&self) -> &[u8] {
        self.slice.data()
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        self.slice.current_block()
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.slice.buffer()
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.slice.is_end_of_buffer()
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        self.slice.skip_to_newline()
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.slice.first_byte()
    }
}

pub trait FromMmap<'a>: FromInputData<'a, MmapInput<'a>> {
    /// Build the struct from a memory mapped file.
    /// It supports parallel processing, but not transparent decompression.
    #[inline(always)]
    fn from_file_mmap<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_input(MmapInput::new(path)?)
    }
}

impl<'a, F: FromInputData<'a, MmapInput<'a>>> FromMmap<'a> for F {}

/// File entirely loaded in RAM, only recommended for small files.
/// It supports parallel processing, but not transparent decompression.
pub struct RamFileInput {
    slice: SliceInput<'static>,
    _vec: Vec<u8>,
}

impl RamFileInput {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let _vec = std::fs::read(path)?;
        let data = unsafe { std::slice::from_raw_parts(_vec.as_ptr(), _vec.len()) };
        Ok(Self {
            slice: SliceInput::new(data),
            _vec,
        })
    }
}

impl Iterator for RamFileInput {
    type Item = &'static [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.slice.next()
    }
}

impl InputData<'static> for RamFileInput {
    const RANDOM_ACCESS: bool = true;

    #[inline(always)]
    fn data(&self) -> &[u8] {
        self.slice.data()
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        self.slice.current_block()
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.slice.buffer()
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.slice.is_end_of_buffer()
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        self.slice.skip_to_newline()
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.slice.first_byte()
    }
}

pub trait FromRamFile: FromInputData<'static, RamFileInput> {
    /// Build the struct from a file entirely loaded in RAM, this is only recommended for small files.
    /// It supports parallel processing, but not transparent decompression.
    #[inline(always)]
    fn from_file_in_ram<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_input(RamFileInput::new(path)?)
    }
}

impl<F: FromInputData<'static, RamFileInput>> FromRamFile for F {}

/// Reader input.
/// It supports transparent decompression, but not parallel processing.
pub struct ReaderInput<'a, R: Read + Send + 'a> {
    data: Vec<u8>,
    len: usize,
    pos: usize,
    offset: usize,
    first_byte: u8,
    decoder: AnyDecoder<R>,
    _phantom: PhantomData<&'a ()>,
}

impl<'a, R: Read + Send + 'a> ReaderInput<'a, R> {
    pub fn new(reader: R) -> Self {
        let mut decoder = AnyDecoder::new(reader);
        let mut data = vec![0; BUFFER_DEFAULT_SIZE];
        let len = decoder
            .read(&mut data[..64])
            .expect("Error while reading data");
        let first_byte = data[0];
        Self {
            data,
            len,
            pos: 0,
            offset: 0,
            first_byte,
            decoder,
            _phantom: PhantomData,
        }
    }
}

impl<'a, R: Read + Send + 'a> Iterator for ReaderInput<'a, R> {
    type Item = &'a [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.len {
            // Save state so we can restore it on EOF, keeping the buffer valid for
            // zero-copy callers (e.g. get_dna_string) that run after next() returns None.
            let saved_offset = self.offset;
            let saved_pos = self.pos;
            let saved_len = self.len;
            self.offset += self.len;
            self.pos = 0;
            self.len = 0;
            while self.len < self.data.len() {
                let n = self
                    .decoder
                    .read(&mut self.data[self.len..])
                    .expect("Error while reading data");
                if n == 0 {
                    break;
                }
                self.len += n;
            }
            if self.len == 0 {
                // EOF: roll back so buffer() still returns the last valid window.
                self.offset = saved_offset;
                self.pos = saved_pos;
                self.len = saved_len;
                return None;
            }
            let padded_len = self.len.next_multiple_of(64);
            self.data[self.len..padded_len].fill(0);
        }
        let pos = self.pos;
        self.pos += 64;
        if pos + 64 <= self.len {
            unsafe { Some(std::slice::from_raw_parts(self.data.as_ptr().add(pos), 64)) }
        } else {
            unsafe {
                Some(std::slice::from_raw_parts(
                    self.data.as_ptr().add(pos),
                    self.len % 64,
                ))
            }
        }
    }
}

impl<'a, R: Read + Send + 'a> InputData<'a> for ReaderInput<'a, R> {
    const RANDOM_ACCESS: bool = false;

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        if 64 <= self.pos && self.pos <= self.len {
            unsafe { std::slice::from_raw_parts(self.data.as_ptr().add(self.pos - 64), 64) }
        } else {
            unsafe {
                std::slice::from_raw_parts(
                    self.data.as_ptr().add((self.len / 64) * 64),
                    self.len % 64,
                )
            }
        }
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        &self.data[..self.len]
    }

    #[inline(always)]
    fn buffer_offset(&self) -> usize {
        self.offset
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.pos >= self.len
    }

    #[inline(always)]
    fn make_room(&mut self, keep_from: usize) {
        // Round the shift down to a 64-byte multiple so that `pos` (which is
        // always a multiple of 64) stays aligned after the copy.  A non-aligned
        // shift would leave `pos` non-aligned, causing `next()` to emit 0-length
        // sentinel blocks at buffer boundaries and desynchronise the parser.
        let shift = keep_from.saturating_sub(self.offset);
        let aligned_shift = (shift / 64) * 64;
        if aligned_shift > 0 {
            self.data.copy_within(aligned_shift..self.len, 0);
            self.len -= aligned_shift;
            self.pos -= aligned_shift;
            self.offset += aligned_shift;
        } else {
            // The shift is smaller than 64 bytes (or keep_from is already behind
            // the buffer start).  Grow the buffer instead so there is room to
            // refill without overwriting data we still need.
            // Policy mirrors needletail: double until 8 MiB, then +8 MiB per step.
            let new_size = if self.data.len() < BUFFER_DOUBLE_UNTIL {
                self.data.len() * 2
            } else {
                self.data.len() + BUFFER_DOUBLE_UNTIL
            };
            self.data.resize(new_size, 0);
        }
        while self.len < self.data.len() {
            let n = self
                .decoder
                .read(&mut self.data[self.len..])
                .expect("Error while reading data");
            if n == 0 {
                break;
            }
            self.len += n;
        }
        let padded = self.len.next_multiple_of(64);
        self.data[self.len..padded].fill(0);
    }

    #[inline(always)]
    fn grow_buffer(&mut self, additional: usize) {
        self.data.resize(self.len + additional, 0);
        let n = self
            .decoder
            .read(&mut self.data[self.len..])
            .expect("Error while reading data");
        self.len += n;
        let padded_len = self.len.next_multiple_of(64);
        self.data[self.len..padded_len].fill(0);
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        // self.pos is always a multiple of 64 and points to the start of the
        // next block that would be returned by Iterator::next().  Scan from
        // there through the end of valid data for '\n'.
        if self.pos >= self.len {
            return None;
        }
        let p = memchr(b'\n', &self.data[self.pos..self.len])?;
        // p is a byte offset from self.pos.
        let delta_blocks = p / 64; // full 64-byte blocks before the \n block
        let pos_in_block = p % 64; // position of \n within its block
        let block_start = self.pos + delta_blocks * 64;
        let block_end = (block_start + 64).min(self.len);
        let block_len = block_end - block_start;
        // Advance past the block that contains the newline so the next
        // Iterator::next() call (from consume_newline → increment_pos) will
        // correctly load the block that follows the quality line.
        self.pos = block_start + 64;
        Some((delta_blocks, pos_in_block, block_len))
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.first_byte
    }

    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        let format = self.decoder.kind()?;
        if format == deko::Format::Verbatim {
            Ok(None)
        } else {
            Ok(Some(format))
        }
    }
}

pub trait FromReader<'a, R: Read + Send + 'a>: FromInputData<'a, ReaderInput<'a, R>> {
    /// Build the struct from a reader.
    /// It supports transparent decompression, but not parallel processing.
    #[inline(always)]
    fn from_reader(reader: R) -> io::Result<Self> {
        Self::from_input(ReaderInput::new(reader))
    }
}

impl<'a, R: Read + Send + 'a, F: FromInputData<'a, ReaderInput<'a, R>>> FromReader<'a, R> for F {}

/// File input.
/// It supports transparent decompression, but not parallel processing.
pub struct FileInput {
    reader: ReaderInput<'static, File>,
}

impl FileInput {
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            reader: ReaderInput::new(File::open(path)?),
        })
    }
}

impl Iterator for FileInput {
    type Item = &'static [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next()
    }
}

impl InputData<'static> for FileInput {
    const RANDOM_ACCESS: bool = false;

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        self.reader.current_block()
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.reader.buffer()
    }

    #[inline(always)]
    fn buffer_offset(&self) -> usize {
        self.reader.buffer_offset()
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.reader.is_end_of_buffer()
    }

    #[inline(always)]
    fn make_room(&mut self, keep_from: usize) {
        self.reader.make_room(keep_from)
    }

    #[inline(always)]
    fn grow_buffer(&mut self, additional: usize) {
        self.reader.grow_buffer(additional)
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        self.reader.skip_to_newline()
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.reader.first_byte()
    }

    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        self.reader.compression_format()
    }
}

pub trait FromFile: FromInputData<'static, FileInput> {
    /// Build the struct from a file.
    /// It supports transparent decompression, but not parallel processing.
    #[inline(always)]
    fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_input(FileInput::new(path)?)
    }
}

impl<F: FromInputData<'static, FileInput>> FromFile for F {}

/// Stdin input.
/// It supports transparent decompression, but not parallel processing.
pub struct StdinInput {
    reader: ReaderInput<'static, Stdin>,
}

impl StdinInput {
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            reader: ReaderInput::new(stdin()),
        }
    }
}

impl Iterator for StdinInput {
    type Item = &'static [u8];

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next()
    }
}

impl InputData<'static> for StdinInput {
    const RANDOM_ACCESS: bool = false;

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        self.reader.current_block()
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.reader.buffer()
    }

    #[inline(always)]
    fn buffer_offset(&self) -> usize {
        self.reader.buffer_offset()
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.reader.is_end_of_buffer()
    }

    #[inline(always)]
    fn make_room(&mut self, keep_from: usize) {
        self.reader.make_room(keep_from)
    }

    #[inline(always)]
    fn grow_buffer(&mut self, additional: usize) {
        self.reader.grow_buffer(additional)
    }

    #[inline(always)]
    fn skip_to_newline(&mut self) -> Option<(usize, usize, usize)> {
        self.reader.skip_to_newline()
    }

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.reader.first_byte()
    }

    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        self.reader.compression_format()
    }
}

pub trait FromStdin: FromInputData<'static, StdinInput> {
    /// Build the struct from stdin.
    /// It supports transparent decompression, but not parallel processing.
    #[inline(always)]
    fn from_stdin() -> io::Result<Self> {
        Self::from_input(StdinInput::new())
    }
}

impl<F: FromInputData<'static, StdinInput>> FromStdin for F {}
