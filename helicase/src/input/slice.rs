#[allow(deprecated)]
use super::FromInputData;
use super::InputBuffer;
use crate::parser::Result;

#[cfg(feature = "std")]
use {
    super::CompressionFormat, crate::parser::Error, memmap2::Mmap, std::fs::File, std::path::Path,
};

/// Slice input.
/// It does not support transparent decompression.
pub struct SliceInput<'a> {
    data: &'a [u8],
    pos: usize,
    first_byte: u8,
    last_block: [u8; 64],
}

impl<'a> SliceInput<'a> {
    /// Creates an input reading `data`.
    pub(crate) fn new(data: &'a [u8]) -> Self {
        let mut last_block = [0; 64];
        last_block[..data.len() % 64].copy_from_slice(&data[(data.len() / 64) * 64..]);
        Self {
            data,
            pos: 0,
            first_byte: data.first().copied().unwrap_or(0),
            last_block,
        }
    }
}

impl<'a> Iterator for SliceInput<'a> {
    type Item = Result<&'a [u8]>;

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

impl<'a> InputBuffer<'a> for SliceInput<'a> {
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
            &self.last_block
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
    fn compression_format(&mut self) -> Result<Option<super::CompressionFormat>> {
        Ok(super::CompressionFormat::from_bytes(self.data))
    }

    #[inline(always)]
    fn seek_block(&mut self, block_idx: usize) {
        self.pos = (block_idx + 1) * 64;
    }
}

/// Types that can be built from a slice, such as parsers.
#[deprecated(note = "use the `from_slice` method of the parsers")]
#[allow(deprecated)]
pub trait FromSlice<'a>: FromInputData<'a, SliceInput<'a>> {
    /// Build the struct from a slice.
    /// It does not support transparent decompression.
    #[inline(always)]
    fn from_slice(data: &'a [u8]) -> Result<Self> {
        Self::from_input(SliceInput::new(data))
    }
}

#[allow(deprecated)]
impl<'a, F: FromInputData<'a, SliceInput<'a>>> FromSlice<'a> for F {}

/// Memory-mapped file.
///
/// Provides zero-copy parsing over memory-mapped files without decompression overhead.
/// For transparent decompression on compressed inputs, use [`FastxParser::from_file`](crate::parser::FastxParser::from_file)
/// or [`FastxParser::from_file_mmap`](crate::parser::FastxParser::from_file_mmap).
#[cfg(feature = "std")]
pub struct MmapInput<'a> {
    slice: SliceInput<'a>,
    _mmap: Mmap,
}

#[cfg(feature = "std")]
impl<'a> MmapInput<'a> {
    /// Creates an input from an existing memory map.
    pub(crate) fn from_mmap(mmap: Mmap) -> Self {
        let data = unsafe { core::slice::from_raw_parts(mmap.as_ptr(), mmap.len()) };
        Self {
            slice: SliceInput::new(data),
            _mmap: mmap,
        }
    }

    /// Memory-maps the file at `path`, which must not be compressed.
    pub(crate) fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if let Some(format) = CompressionFormat::from_bytes(&mmap) {
            return Err(Error::custom(format!(
                "input is compressed ({format}), use `from_file` or `FastxParser::from_file_mmap`"
            )));
        }
        Ok(Self::from_mmap(mmap))
    }
}

#[cfg(feature = "std")]
impl<'a> Iterator for MmapInput<'a> {
    type Item = Result<&'a [u8]>;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.slice.next()
    }
}

#[cfg(feature = "std")]
impl<'a> InputBuffer<'a> for MmapInput<'a> {
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

/// Types that can be built from a memory-mapped file, such as parsers.
#[cfg(feature = "std")]
#[deprecated(note = "use the `from_file_mmap` method of the parsers")]
#[allow(deprecated)]
pub trait FromMmap<'a>: FromInputData<'a, MmapInput<'a>> {
    /// Build the struct from a memory-mapped file.
    ///
    /// Memory maps the file for zero-copy parsing on uncompressed inputs.
    #[inline(always)]
    fn from_file_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_input(MmapInput::new(path)?)
    }
}

#[cfg(feature = "std")]
#[allow(deprecated)]
impl<'a, F: FromInputData<'a, MmapInput<'a>>> FromMmap<'a> for F {}
