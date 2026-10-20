//! Input types and helpers.

use core::marker::PhantomData;
use core::ops::Range;
use deko::read::AnyDecoder;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

const BUFFER_DEFAULT_SIZE: usize = 1 << 17;
const BUFFER_DOUBLE_UNTIL: usize = 1 << 24;

pub trait InputData<'a>: Iterator<Item = io::Result<&'a [u8]>> {
    const RANDOM_ACCESS: bool;

    /// Returns the first byte of the (uncompressed when possible) input.
    fn first_byte(&self) -> u8;

    /// Get a reference to the current block of up to 64 bytes.
    ///
    /// If the length is smaller than 64, the following bytes are guaranteed to be zeros.
    fn current_block(&self) -> &[u8];

    /// Returns the type of compression format detected, `None` if uncompressed.
    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        Ok(None)
    }

    /// Returns `true` if compression has been detected.
    #[inline(always)]
    fn is_compressed(&mut self) -> io::Result<bool> {
        Ok(self.compression_format()?.is_some())
    }

    /// Get a reference to the complete slice of data.
    ///
    /// This is not available for reader-based implementations.
    #[inline(always)]
    fn data(&self) -> &'a [u8] {
        assert!(Self::RANDOM_ACCESS);
        unimplemented!()
    }

    /// Get a reference to the internal buffer.
    ///
    /// This is only relevant for reader-based implementations.
    fn buffer(&self) -> &[u8];

    /// Returns the offset of the buffer.
    ///
    /// This is only relevant for reader-based implementations.
    #[inline(always)]
    fn buffer_offset(&self) -> usize {
        0
    }

    /// Returns a slice of the buffered data for the given global byte range.
    #[inline(always)]
    fn slice(&self, range: Range<usize>) -> &[u8] {
        let off = self.buffer_offset();
        &self.buffer()[range.start - off..range.end - off]
    }

    /// Returns `true` if we are reaching the end of the buffer.
    ///
    /// This is only relevant for reader-based implementations.
    fn is_end_of_buffer(&self) -> bool;

    /// Shift buffered data from `keep_from` to the start of the buffer, then refill from the reader.
    ///
    /// This is only relevant for reader-based implementations.
    #[doc(hidden)]
    #[inline(always)]
    fn make_room(&mut self, _keep_from: usize) {}

    /// Seek internal buffer position to the start of the block after `block_idx`.
    #[doc(hidden)]
    #[inline(always)]
    fn seek_block(&mut self, _block_idx: usize) {}
}

pub trait FromInputData<'a, I: InputData<'a>>: Sized {
    /// Build the struct from a type implementing [`InputData`].
    fn from_input(input: I) -> io::Result<Self>;
}

/// Slice input.
/// It does not support transparent decompression.
pub struct SliceInput<'a> {
    data: &'a [u8],
    pos: usize,
    first_byte: u8,
    last_block: [u8; 64],
}

impl<'a> SliceInput<'a> {
    pub fn new(data: &'a [u8]) -> io::Result<Self> {
        if data.is_empty() {
            return Err(io::Error::other("Input data is empty"));
        }
        let mut last_block = [0; 64];
        last_block[..data.len() % 64].copy_from_slice(&data[(data.len() / 64) * 64..]);
        Ok(Self {
            data,
            pos: 0,
            first_byte: data[0],
            last_block,
        })
    }
}

impl<'a> Iterator for SliceInput<'a> {
    type Item = io::Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos;
        self.pos += 64;
        if pos + 64 <= self.data.len() {
            unsafe {
                Some(Ok(core::slice::from_raw_parts(
                    self.data.as_ptr().add(pos),
                    64,
                )))
            }
        } else if pos < self.data.len() {
            unsafe {
                Some(Ok(core::slice::from_raw_parts(
                    self.last_block.as_ptr(),
                    self.data.len() % 64,
                )))
            }
        } else {
            None
        }
    }
}

impl<'a> InputData<'a> for SliceInput<'a> {
    const RANDOM_ACCESS: bool = true;

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.first_byte
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        if 64 <= self.pos && self.pos <= self.data.len() {
            unsafe { core::slice::from_raw_parts(self.data.as_ptr().add(self.pos - 64), 64) }
        } else {
            unsafe { core::slice::from_raw_parts(self.last_block.as_ptr(), self.data.len() % 64) }
        }
    }

    #[inline(always)]
    fn data(&self) -> &'a [u8] {
        self.data
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
    fn seek_block(&mut self, block_idx: usize) {
        self.pos = (block_idx + 1) * 64;
    }
}

pub trait FromSlice<'a>: FromInputData<'a, SliceInput<'a>> {
    /// Build the struct from a slice.
    /// It does not support transparent decompression.
    #[inline(always)]
    fn from_slice(data: &'a [u8]) -> io::Result<Self> {
        Self::from_input(SliceInput::new(data)?)
    }
}

impl<'a, F: FromInputData<'a, SliceInput<'a>>> FromSlice<'a> for F {}

/// Memory-mapped file.
/// It does not support transparent decompression.
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
        let data = unsafe { core::slice::from_raw_parts(_mmap.as_ptr(), _mmap.len()) };
        Ok(Self {
            slice: SliceInput::new(data)?,
            _mmap,
        })
    }
}

impl<'a> Iterator for MmapInput<'a> {
    type Item = io::Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.slice.next()
    }
}

impl<'a> InputData<'a> for MmapInput<'a> {
    const RANDOM_ACCESS: bool = true;

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.slice.first_byte()
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        self.slice.current_block()
    }

    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        self.slice.compression_format()
    }

    #[inline(always)]
    fn data(&self) -> &'a [u8] {
        self.slice.data()
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        self.slice.buffer()
    }

    #[inline(always)]
    fn buffer_offset(&self) -> usize {
        self.slice.buffer_offset()
    }

    #[inline(always)]
    fn is_end_of_buffer(&self) -> bool {
        self.slice.is_end_of_buffer()
    }

    #[inline(always)]
    fn make_room(&mut self, keep_from: usize) {
        self.slice.make_room(keep_from);
    }

    #[inline(always)]
    fn seek_block(&mut self, block_idx: usize) {
        self.slice.seek_block(block_idx);
    }
}

pub trait FromMmap<'a>: FromInputData<'a, MmapInput<'a>> {
    /// Build the struct from a memory-mapped file.
    /// It does not support transparent decompression.
    ///
    /// Also see [`from_file`](FromFile::from_file) to use a file reader.
    #[inline(always)]
    fn from_file_mmap<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_input(MmapInput::new(path)?)
    }
}

impl<'a, F: FromInputData<'a, MmapInput<'a>>> FromMmap<'a> for F {}

#[cfg(feature = "sirocco")]
const SIROCCO_BUF_LEN: usize = 256 * 1024;

#[cfg(feature = "sirocco")]
struct SiroccoMultiGzReader<R: BufRead> {
    inner: R,
    decoder: sirocco::inflate::Decoder,
    buf: Vec<u8>,
    pending: Range<usize>,
    prev_total_out: usize,
    full: bool,
    done: bool,
    pending_err: Option<io::Error>,
}

#[cfg(feature = "sirocco")]
impl<R: BufRead> SiroccoMultiGzReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            decoder: sirocco::inflate::Decoder::new(sirocco::Format::Gzip),
            buf: Vec::new(),
            pending: 0..0,
            prev_total_out: 0,
            full: false,
            done: false,
            pending_err: None,
        }
    }

    fn fill(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            self.buf = vec![0; SIROCCO_BUF_LEN];
        }

        while self.pending.is_empty() && self.pending_err.is_none() {
            if self.done {
                let peek = self.inner.fill_buf()?;
                if peek.is_empty() {
                    break;
                }
                let is_next_gz = if peek.len() >= 2 {
                    peek[0] == 0x1F && peek[1] == 0x8B
                } else {
                    peek[0] == 0x1F
                };
                if is_next_gz {
                    self.decoder.reset();
                    self.prev_total_out = 0;
                    self.pending = 0..0;
                    self.full = false;
                    self.done = false;
                } else {
                    break;
                }
            }

            if self.full {
                self.decoder.make_room(&mut self.buf);
                self.full = false;
            }

            let input = self.inner.fill_buf()?;
            let eof = input.is_empty();
            let buf_ptr = self.buf.as_ptr() as usize;

            match self.decoder.decompress(input, &mut self.buf, eof) {
                Ok(p) => {
                    let out_start = p.output.as_ptr() as usize - buf_ptr;
                    let out_end = out_start + p.output.len();
                    (self.done, self.full) = (p.done, p.full);
                    self.inner.consume(p.consumed);
                    self.pending = out_start..out_end;
                    self.prev_total_out = self.decoder.total_out() as usize;
                }
                Err(e) => {
                    let new_total = self.decoder.total_out() as usize;
                    let new_bytes = new_total.saturating_sub(self.prev_total_out);
                    if new_bytes > 0 {
                        let out_start = self.pending.end;
                        self.pending = out_start..out_start + new_bytes;
                        self.prev_total_out = new_total;
                        self.done = true;
                        self.pending_err = Some(io::Error::new(io::ErrorKind::InvalidData, e));
                    } else {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, e));
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(feature = "sirocco")]
impl<R: BufRead> Read for SiroccoMultiGzReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        if self.pending.is_empty() {
            if let Some(err) = self.pending_err.take() {
                return Err(err);
            }
            self.fill()?;
        }

        if self.pending.is_empty() {
            if let Some(err) = self.pending_err.take() {
                return Err(err);
            }
            return Ok(0);
        }

        let n = self.pending.len().min(buf.len());
        buf[..n].copy_from_slice(&self.buf[self.pending.start..self.pending.start + n]);
        self.pending.start += n;
        Ok(n)
    }
}

enum ReaderDecoder<R: Read> {
    #[cfg(feature = "sirocco")]
    Sirocco(Box<SiroccoMultiGzReader<BufReader<R>>>),
    Deko(Box<AnyDecoder<BufReader<R>>>),
}

impl<R: Read> Read for ReaderDecoder<R> {
    #[inline(always)]
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            #[cfg(feature = "sirocco")]
            ReaderDecoder::Sirocco(r) => r.read(buf),
            ReaderDecoder::Deko(r) => r.read(buf),
        }
    }
}

/// Reader input.
/// It supports transparent decompression.
pub struct ReaderInput<'a, R: Read + Send + 'a> {
    data: Vec<u8>,
    len: usize,
    pos: usize,
    offset: usize,
    first_byte: u8,
    decoder: ReaderDecoder<R>,
    pending_err: Option<io::Error>,
    errored: bool,
    _phantom: PhantomData<&'a ()>,
}

impl<'a, R: Read + Send + 'a> ReaderInput<'a, R> {
    fn fill_buffer(&mut self) {
        while self.pending_err.is_none() && self.len < self.data.len() {
            match self.decoder.read(&mut self.data[self.len..]) {
                Ok(0) => break,
                Ok(n) => self.len += n,
                Err(e) => {
                    self.pending_err = Some(e);
                    break;
                }
            }
        }
        let padded = self.len.next_multiple_of(64);
        self.data[self.len..padded].fill(0);
    }

    pub fn new(reader: R) -> io::Result<Self> {
        let mut buf_reader = BufReader::with_capacity(64 * 1024, reader);
        let is_gz = {
            let peek = buf_reader.fill_buf()?;
            peek.len() >= 2 && peek[0] == 0x1F && peek[1] == 0x8B
        };

        #[cfg(feature = "sirocco")]
        let decoder = if is_gz {
            ReaderDecoder::Sirocco(Box::new(SiroccoMultiGzReader::new(buf_reader)))
        } else {
            ReaderDecoder::Deko(Box::new(AnyDecoder::new(buf_reader)))
        };

        #[cfg(not(feature = "sirocco"))]
        let decoder = ReaderDecoder::Deko(Box::new(AnyDecoder::new(buf_reader)));

        let mut s = Self {
            data: vec![0; BUFFER_DEFAULT_SIZE],
            len: 0,
            pos: 0,
            offset: 0,
            first_byte: 0,
            decoder,
            pending_err: None,
            errored: false,
            _phantom: PhantomData,
        };
        s.fill_buffer();
        if s.len == 0
            && let Some(e) = s.pending_err.take()
        {
            return Err(e);
        }
        if s.len > 0 {
            s.first_byte = s.data[0];
        }
        Ok(s)
    }

    #[cold]
    #[inline(never)]
    fn next_slow(&mut self) -> Option<io::Result<&'a [u8]>> {
        if self.errored {
            return None;
        }
        if self.pos >= self.len {
            if let Some(err) = self.pending_err.take() {
                self.errored = true;
                return Some(Err(err));
            }
            // Save state so we can restore it on EOF, keeping the buffer valid for
            // zero-copy callers (e.g. get_dna_string) that run after next() returns None.
            let (saved_offset, saved_pos, saved_len) = (self.offset, self.pos, self.len);
            self.offset += self.len;
            self.pos = 0;
            self.len = 0;
            self.fill_buffer();
            if self.len == 0 {
                // EOF or error: roll back so buffer() still returns the last valid window.
                (self.offset, self.pos, self.len) = (saved_offset, saved_pos, saved_len);
                if let Some(err) = self.pending_err.take() {
                    self.errored = true;
                    return Some(Err(err));
                }
                return None;
            }
        }
        let pos = self.pos;
        self.pos += 64;
        let len = (self.len - pos).min(64);
        unsafe {
            Some(Ok(core::slice::from_raw_parts(
                self.data.as_ptr().add(pos),
                len,
            )))
        }
    }
}

impl<'a, R: Read + Send + 'a> Iterator for ReaderInput<'a, R> {
    type Item = io::Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos;
        if pos + 64 <= self.len {
            self.pos = pos + 64;
            return unsafe {
                Some(Ok(core::slice::from_raw_parts(
                    self.data.as_ptr().add(pos),
                    64,
                )))
            };
        }
        self.next_slow()
    }
}

impl<'a, R: Read + Send + 'a> InputData<'a> for ReaderInput<'a, R> {
    const RANDOM_ACCESS: bool = false;

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.first_byte
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        unsafe { core::slice::from_raw_parts(self.data.as_ptr().add(self.pos - 64), 64) }
    }

    #[inline(always)]
    fn compression_format(&mut self) -> io::Result<Option<deko::Format>> {
        match &mut self.decoder {
            #[cfg(feature = "sirocco")]
            ReaderDecoder::Sirocco(_) => {
                // In deko 0.6, deko::Format::Gz is only defined under deko's flate2 feature.
                // Sirocco handles decompression directly; is_compressed() returns true.
                Ok(None)
            }
            ReaderDecoder::Deko(decoder) => {
                let format = decoder.kind()?;
                if format == deko::Format::Verbatim {
                    Ok(None)
                } else {
                    Ok(Some(format))
                }
            }
        }
    }

    #[inline(always)]
    fn is_compressed(&mut self) -> io::Result<bool> {
        match &mut self.decoder {
            #[cfg(feature = "sirocco")]
            ReaderDecoder::Sirocco(_) => Ok(true),
            ReaderDecoder::Deko(decoder) => {
                let format = decoder.kind()?;
                Ok(format != deko::Format::Verbatim)
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

    #[cold]
    #[inline(never)]
    fn make_room(&mut self, keep_from: usize) {
        // Round the shift down to a 64-byte multiple so that `pos` (which is
        // always a multiple of 64) stays aligned after the copy.
        let aligned_shift = (keep_from.saturating_sub(self.offset) / 64) * 64;
        if aligned_shift > 0 {
            self.data.copy_within(aligned_shift..self.len, 0);
            self.len -= aligned_shift;
            self.pos -= aligned_shift;
            self.offset += aligned_shift;
        } else if self.len == self.data.len() {
            let new_size = if self.data.len() < BUFFER_DOUBLE_UNTIL {
                self.data.len() * 2
            } else {
                self.data.len() + BUFFER_DOUBLE_UNTIL
            };
            self.data.resize(new_size, 0);
        }
        self.fill_buffer();
    }

    #[inline(always)]
    fn seek_block(&mut self, block_idx: usize) {
        self.pos = (block_idx + 1) * 64 - self.offset;
    }
}

pub trait FromReader<'a, R: Read + Send + 'a>: FromInputData<'a, ReaderInput<'a, R>> {
    /// Build the struct from a reader.
    /// It supports transparent decompression.
    #[inline(always)]
    fn from_reader(reader: R) -> io::Result<Self> {
        Self::from_input(ReaderInput::new(reader)?)
    }
}

impl<'a, R: Read + Send + 'a, F: FromInputData<'a, ReaderInput<'a, R>>> FromReader<'a, R> for F {}

/// File input alias.
/// It supports transparent decompression.
pub type FileInput = ReaderInput<'static, File>;

impl ReaderInput<'static, File> {
    /// Open a file for reading with transparent decompression.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::new(File::open(path)?)
    }
}

pub trait FromFile: FromInputData<'static, FileInput> {
    /// Build the struct from a file.
    /// It supports transparent decompression.
    ///
    /// Also see [`from_file_mmap`](FromMmap::from_file_mmap) to use a memory-mapped file.
    #[inline(always)]
    fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_input(FileInput::open(path)?)
    }
}

impl<F: FromInputData<'static, FileInput>> FromFile for F {}
