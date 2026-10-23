use super::{CompressionFormat, FromInputData, InputData, SliceInput};
use crate::parser::Result;
use core::marker::PhantomData;
#[cfg(feature = "sirocco")]
use core::ops::Range;
use deko::read::AnyDecoder;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::mem::{replace, take};
use std::path::Path;
use std::sync::mpsc::{Receiver, Sender, SyncSender, channel, sync_channel};

// Bytes read into a buffer at once, and initial size of the direct reader's buffer.
const BUFFER_LEN: usize = 1 << 18;
// Buffers for longer records double until this size, then grow linearly.
const BUFFER_DOUBLE_UNTIL: usize = 1 << 24;
// The read-ahead reader (`FileInput`) keeps up to `BUFFERS_AHEAD` buffers ready, each read after
// `BUFFER_ROOM` free bytes. The room is used to prepend the end of the previous buffer without
// copying the new one; longer ends of records get the new buffer appended instead.
const BUFFERS_AHEAD: usize = 2;
const BUFFER_ROOM: usize = BUFFER_LEN / 2;

impl CompressionFormat {
    /// Detects compression format from a file path by peeking magic bytes and checking file extension.
    pub fn from_path(path: impl AsRef<Path>) -> io::Result<Option<Self>> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let mut buf = [0u8; 64];
        let n = (&file).read(&mut buf)?;
        if n > 0 {
            let format = Self::from_bytes(&buf[..n]);
            if format == Some(Self::Zstd) {
                #[cfg(feature = "zstd")]
                {
                    if let Ok(mmap) = unsafe { Mmap::map(&file) } {
                        let frames = crate::parallel::index_zstd_frames(&mmap);
                        if frames.len() > 1 {
                            return Ok(Some(Self::Pzstd));
                        }
                    }
                }
                return Ok(Some(Self::Zstd));
            }
            if format.is_some() {
                return Ok(format);
            }
        }

        // Fallback: check file extension
        if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
            match ext.to_ascii_lowercase().as_str() {
                "bgz" | "bgzf" => return Ok(Some(Self::Bgzf)),
                "gz" | "gzip" => return Ok(Some(Self::Gzip)),
                "zst" | "zstd" => return Ok(Some(Self::Zstd)),
                "pzst" | "pzstd" => return Ok(Some(Self::Pzstd)),
                "bz2" => return Ok(Some(Self::Bzip2)),
                "xz" => return Ok(Some(Self::Xz)),
                _ => {}
            }
        }

        Ok(None)
    }
}

/// Memory-mapped file.
/// It does not support transparent decompression.
pub struct MmapInput<'a> {
    slice: SliceInput<'a>,
    _mmap: Mmap,
}

impl<'a> MmapInput<'a> {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        // Unsafe: mmap are intrisically unsafe.
        // Here we add on the top of that a slice of [u8] that will live as long as file is not dropped.
        // To enforce this, we keep the file in the struct.
        let file = File::open(path)?;
        let _mmap = unsafe { Mmap::map(&file)? };
        let data = unsafe { core::slice::from_raw_parts(_mmap.as_ptr(), _mmap.len()) };
        Ok(Self {
            slice: SliceInput::new(data),
            _mmap,
        })
    }
}

impl<'a> Iterator for MmapInput<'a> {
    type Item = Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.slice.next()
    }
}

impl crate::sealed::Sealed for MmapInput<'_> {}

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
    fn compression_format(&mut self) -> Result<Option<CompressionFormat>> {
        self.slice.compression_format()
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
    fn from_file_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
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
    const fn new(inner: R) -> Self {
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
            if self.pending.is_empty() {
                return self.pending_err.take().map_or(Ok(0), Err);
            }
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
    #[cfg(feature = "zstd")]
    Zstd(Box<zstd::stream::read::Decoder<'static, BufReader<R>>>),
    Deko(Box<AnyDecoder<BufReader<R>>>),
}

impl<R: Read> Read for ReaderDecoder<R> {
    #[inline(always)]
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            #[cfg(feature = "sirocco")]
            Self::Sirocco(r) => r.read(buf),
            #[cfg(feature = "zstd")]
            Self::Zstd(r) => r.read(buf),
            Self::Deko(r) => r.read(buf),
        }
    }
}

enum Source<R: Read> {
    Direct(ReaderDecoder<R>),
    ReadAhead(BufferReceiver),
}

// Buffers filled by the background thread, with their data at `buf[BUFFER_ROOM..][..len]`,
// the channel to give used buffers back, and a buffer for records longer than the room.
struct BufferReceiver {
    filled: Receiver<io::Result<(Vec<u8>, usize)>>,
    free: Sender<Vec<u8>>,
    long: Vec<u8>,
    in_long: bool,
}

impl BufferReceiver {
    fn recycle(&self, buf: Vec<u8>) {
        if buf.len() == BUFFER_ROOM + BUFFER_LEN + 64 {
            let _ = self.free.send(buf);
        }
    }
}

/// Reader input.
/// It supports transparent decompression.
pub struct ReaderInput<'a, R: Read + Send + 'a> {
    data: Vec<u8>,
    start: usize,
    len: usize,
    pos: usize,
    offset: usize,
    first_byte: u8,
    source: Source<R>,
    pending_err: Option<io::Error>,
    errored: bool,
    compression: Option<CompressionFormat>,
    _phantom: PhantomData<&'a ()>,
}

impl<'a, R: Read + Send + 'a> ReaderInput<'a, R> {
    // Keeps the buffer from `keep_start` and appends new data after it.
    // Returns whether new data was read.
    fn refill(&mut self, keep_start: usize) -> bool {
        let refilled = match &mut self.source {
            Source::Direct(_) => self.refill_direct(keep_start),
            Source::ReadAhead(_) => self.refill_read_ahead(keep_start),
        };
        let end = self.start + self.len;
        self.data[end..end.next_multiple_of(64)].fill(0);
        refilled
    }

    fn refill_direct(&mut self, keep_start: usize) -> bool {
        let Source::Direct(decoder) = &mut self.source else {
            unreachable!()
        };
        if keep_start > 0 {
            self.data.copy_within(keep_start..self.len, 0);
            self.len -= keep_start;
            self.pos -= keep_start;
            self.offset += keep_start;
        } else if self.len == self.data.len() {
            self.data.resize(grown_len(self.data.len(), 0), 0);
        }
        let prev_len = self.len;
        while self.pending_err.is_none() && self.len < self.data.len() {
            match decoder.read(&mut self.data[self.len..]) {
                Ok(0) => break,
                Ok(n) => self.len += n,
                Err(e) => {
                    self.pending_err = Some(e);
                    break;
                }
            }
        }
        self.len > prev_len
    }

    fn refill_read_ahead(&mut self, keep_start: usize) -> bool {
        let Source::ReadAhead(receiver) = &mut self.source else {
            unreachable!()
        };
        let (buffer, filled_len) = match receiver.filled.recv() {
            Ok(Ok(buffer)) => buffer,
            Ok(Err(e)) => {
                self.pending_err = Some(e);
                return false;
            }
            Err(_) => return false,
        };
        let kept = self.start + keep_start..self.start + self.len;
        let kept_len = kept.len();
        let room = BUFFER_ROOM;
        if kept_len <= room.min(filled_len) {
            let mut buffer = buffer;
            self.start = room - kept_len;
            buffer[self.start..room].copy_from_slice(&self.data[kept]);
            let used = replace(&mut self.data, buffer);
            if receiver.in_long {
                receiver.long = used;
                receiver.in_long = false;
            } else {
                receiver.recycle(used);
            }
        } else {
            let needed = kept_len + filled_len + 64;
            if receiver.in_long {
                if self.data.len() < needed {
                    self.data.resize(grown_len(self.data.len(), needed), 0);
                }
                if keep_start > 0 {
                    self.data.copy_within(kept, 0);
                }
            } else {
                let mut long = take(&mut receiver.long);
                if long.len() < needed {
                    long.resize(grown_len(long.len(), needed), 0);
                }
                long[..kept_len].copy_from_slice(&self.data[kept]);
                receiver.recycle(replace(&mut self.data, long));
                receiver.in_long = true;
            }
            self.data[kept_len..kept_len + filled_len]
                .copy_from_slice(&buffer[room..room + filled_len]);
            self.start = 0;
            receiver.recycle(buffer);
        }
        self.offset += keep_start;
        self.pos -= keep_start;
        self.len = kept_len + filled_len;
        true
    }

    pub fn new(reader: R) -> io::Result<Self> {
        let (decoder, compression) = Self::decoder(reader)?;
        Self::with_source(Source::Direct(decoder), compression)
    }

    fn decoder(reader: R) -> io::Result<(ReaderDecoder<R>, Option<CompressionFormat>)> {
        let mut buf_reader = BufReader::with_capacity(64 * 1024, reader);
        let compression = {
            let peek = buf_reader.fill_buf()?;
            CompressionFormat::from_bytes(peek)
        };

        let decoder = match compression {
            #[cfg(feature = "sirocco")]
            Some(CompressionFormat::Gzip | CompressionFormat::Bgzf) => {
                ReaderDecoder::Sirocco(Box::new(SiroccoMultiGzReader::new(buf_reader)))
            }
            #[cfg(feature = "zstd")]
            Some(CompressionFormat::Zstd | CompressionFormat::Pzstd) => ReaderDecoder::Zstd(
                Box::new(zstd::stream::read::Decoder::with_buffer(buf_reader)?),
            ),
            _ => ReaderDecoder::Deko(Box::new(AnyDecoder::new(buf_reader))),
        };

        Ok((decoder, compression))
    }

    fn with_source(source: Source<R>, compression: Option<CompressionFormat>) -> io::Result<Self> {
        let data = match source {
            Source::Direct(_) => vec![0; BUFFER_LEN],
            Source::ReadAhead(_) => vec![0; 64],
        };
        let mut s = Self {
            data,
            start: 0,
            len: 0,
            pos: 0,
            offset: 0,
            first_byte: 0,
            source,
            pending_err: None,
            errored: false,
            compression,
            _phantom: PhantomData,
        };
        s.refill(0);
        if s.len == 0
            && let Some(e) = s.pending_err.take()
        {
            return Err(e);
        }
        if s.len > 0 {
            s.first_byte = s.data[s.start];
        }
        Ok(s)
    }

    #[cold]
    #[inline(never)]
    fn next_slow(&mut self) -> Option<Result<&'a [u8]>> {
        if self.errored {
            return None;
        }
        if self.pos >= self.len {
            if let Some(err) = self.pending_err.take() {
                self.errored = true;
                return Some(Err(err.into()));
            }
            // On EOF, the buffer stays valid for zero-copy callers (e.g. get_dna_string)
            // that run after next() returns None.
            let saved = (self.offset, self.pos, self.len);
            if !self.refill(self.len) {
                (self.offset, self.pos, self.len) = saved;
                if let Some(err) = self.pending_err.take() {
                    self.errored = true;
                    return Some(Err(err.into()));
                }
                return None;
            }
        }
        let pos = self.pos;
        self.pos += 64;
        let len = (self.len - pos).min(64);
        unsafe {
            Some(Ok(core::slice::from_raw_parts(
                self.data.as_ptr().add(self.start + pos),
                len,
            )))
        }
    }
}

impl<R: Read + Send + 'static> ReaderInput<'static, R> {
    /// Build a reader input that reads and decompresses ahead in a background thread.
    pub fn with_read_ahead(reader: R) -> io::Result<Self> {
        let (decoder, compression) = Self::decoder(reader)?;
        let (filled_sender, filled) = sync_channel(BUFFERS_AHEAD);
        let (free, free_receiver) = channel();
        std::thread::spawn(move || read_buffers(decoder, filled_sender, free_receiver));
        let receiver = BufferReceiver {
            filled,
            free,
            long: Vec::new(),
            in_long: false,
        };
        Self::with_source(Source::ReadAhead(receiver), compression)
    }
}

impl<'a, R: Read + Send + 'a> Iterator for ReaderInput<'a, R> {
    type Item = Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos;
        if pos + 64 <= self.len {
            self.pos = pos + 64;
            return unsafe {
                Some(Ok(core::slice::from_raw_parts(
                    self.data.as_ptr().add(self.start + pos),
                    64,
                )))
            };
        }
        self.next_slow()
    }
}

impl<'a, R: Read + Send + 'a> crate::sealed::Sealed for ReaderInput<'a, R> {}

impl<'a, R: Read + Send + 'a> InputData<'a> for ReaderInput<'a, R> {
    const RANDOM_ACCESS: bool = false;

    #[inline(always)]
    fn first_byte(&self) -> u8 {
        self.first_byte
    }

    #[inline(always)]
    fn current_block(&self) -> &[u8] {
        unsafe {
            core::slice::from_raw_parts(self.data.as_ptr().add(self.start + self.pos - 64), 64)
        }
    }

    #[inline(always)]
    fn compression_format(&mut self) -> Result<Option<CompressionFormat>> {
        if self.compression.is_some() {
            return Ok(self.compression);
        }
        if let Source::Direct(ReaderDecoder::Deko(decoder)) = &mut self.source {
            let format = decoder.kind()?;
            let mapped = match format {
                #[cfg(feature = "bz2")]
                deko::Format::Bz => Some(CompressionFormat::Bzip2),
                #[cfg(feature = "xz")]
                deko::Format::Xz => Some(CompressionFormat::Xz),
                _ => None,
            };
            self.compression = mapped;
            return Ok(mapped);
        }
        Ok(None)
    }

    #[inline(always)]
    fn buffer(&self) -> &[u8] {
        &self.data[self.start..self.start + self.len]
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
        self.refill(aligned_shift);
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
    fn from_reader(reader: R) -> Result<Self> {
        Self::from_input(ReaderInput::new(reader)?)
    }
}

impl<'a, R: Read + Send + 'a, F: FromInputData<'a, ReaderInput<'a, R>>> FromReader<'a, R> for F {}

/// File input alias.
/// It supports transparent decompression.
pub type FileInput = ReaderInput<'static, File>;

impl ReaderInput<'static, File> {
    /// Open a file for reading with transparent decompression.
    /// It reads ahead in a background thread when more than one CPU is available.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        if std::thread::available_parallelism().map_or(1, |n| n.get()) > 1 {
            Self::with_read_ahead(file)
        } else {
            Self::new(file)
        }
    }
}

pub trait FromFile: FromInputData<'static, FileInput> {
    /// Build the struct from a file.
    /// It supports transparent decompression.
    ///
    /// Also see [`from_file_mmap`](FromMmap::from_file_mmap) to use a memory-mapped file.
    #[inline(always)]
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_input(FileInput::open(path)?)
    }
}

impl<F: FromInputData<'static, FileInput>> FromFile for F {}

fn grown_len(len: usize, needed: usize) -> usize {
    let grown = if len < BUFFER_DOUBLE_UNTIL {
        2 * len
    } else {
        len + BUFFER_DOUBLE_UNTIL
    };
    grown.max(needed)
}

// Background thread of `ReaderInput::with_read_ahead`, stops at the end of input,
// on error, or when the input is dropped.
fn read_buffers<R: Read>(
    mut decoder: ReaderDecoder<R>,
    filled: SyncSender<io::Result<(Vec<u8>, usize)>>,
    free: Receiver<Vec<u8>>,
) {
    loop {
        let mut buf = free
            .try_recv()
            .unwrap_or_else(|_| vec![0; BUFFER_ROOM + BUFFER_LEN + 64]);
        let room = BUFFER_ROOM;
        let mut len = 0;
        let mut error = None;
        while len < BUFFER_LEN {
            match decoder.read(&mut buf[room + len..room + BUFFER_LEN]) {
                Ok(0) => break,
                Ok(n) => len += n,
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }
        if len > 0 && filled.send(Ok((buf, len))).is_err() {
            return;
        }
        if let Some(e) = error {
            let _ = filled.send(Err(e));
            return;
        }
        if len < BUFFER_LEN {
            return;
        }
    }
}
