//! Input types and helpers.

use crate::parser::Result;
use core::ops::Range;

#[cfg(feature = "std")]
mod reader;
#[cfg(feature = "std")]
pub use reader::*;

/// Supported compression formats for FASTX files and streams.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CompressionFormat {
    /// Standard single-member Gzip stream (.gz)
    Gzip,
    /// Blocked GNU Zip Format (.bgz, .bgzf)
    Bgzf,
    /// Standard single-frame Zstandard stream (.zst)
    Zstd,
    /// Parallel / multi-frame Zstandard stream (.zst, pzstd)
    Pzstd,
    /// Bzip2 stream (.bz2)
    Bzip2,
    /// XZ stream (.xz)
    Xz,
}

impl CompressionFormat {
    /// Returns the canonical name of the compression format.
    #[inline]
    pub const fn as_str(&self) -> &'static str {
        match self {
            Self::Gzip => "gzip",
            Self::Bgzf => "bgzf",
            Self::Zstd => "zstd",
            Self::Pzstd => "pzstd",
            Self::Bzip2 => "bzip2",
            Self::Xz => "xz",
        }
    }

    /// Detects compression format from the initial magic bytes of a buffer/stream.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
            if bytes.len() >= 18
                && bytes[2] == 0x08
                && (bytes[3] & 0x04) != 0
                && bytes[12] == b'B'
                && bytes[13] == b'C'
                && bytes[14] == 0x02
                && bytes[15] == 0x00
            {
                return Some(Self::Bgzf);
            }
            return Some(Self::Gzip);
        }

        // Zstandard: skippable frame used by pzstd (0x184D2A50..0x184D2A5F)
        if bytes.len() >= 4
            && bytes[3] == 0x18
            && bytes[2] == 0x4d
            && bytes[1] == 0x2a
            && (bytes[0] & 0xf0) == 0x50
        {
            return Some(Self::Pzstd);
        }

        // Standard Zstandard frame: 0xFD2FB528 (little-endian: 0x28, 0xb5, 0x2f, 0xfd)
        if bytes.len() >= 4 && bytes[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
            return Some(Self::Zstd);
        }

        // Bzip2: 'B', 'Z'
        if bytes.len() >= 2 && bytes[0] == b'B' && bytes[1] == b'Z' {
            return Some(Self::Bzip2);
        }

        // XZ: 0xFD, '7', 'z', 'X', 'Z', 0x00
        if bytes.len() >= 3 && bytes[0] == 0xfd && bytes[1] == b'7' && bytes[2] == b'z' {
            return Some(Self::Xz);
        }

        None
    }
}

impl core::fmt::Display for CompressionFormat {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.as_str())
    }
}

pub trait InputData<'a>: Iterator<Item = Result<&'a [u8]>> + crate::sealed::Sealed {
    const RANDOM_ACCESS: bool;

    /// Returns the first byte of the (uncompressed when possible) input.
    fn first_byte(&self) -> u8;

    /// Get a reference to the current block of up to 64 bytes.
    ///
    /// If the length is smaller than 64, the following bytes are guaranteed to be zeros.
    fn current_block(&self) -> &[u8];

    /// Returns the type of compression format detected, `None` if uncompressed.
    #[inline(always)]
    fn compression_format(&mut self) -> Result<Option<CompressionFormat>> {
        Ok(None)
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
        // A line feed without `\r` in a CRLF input gives an inverted range: return it empty.
        &self.buffer()[range.start - off..range.end.max(range.start) - off]
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
    fn from_input(input: I) -> Result<Self>;
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
    pub fn new(data: &'a [u8]) -> Self {
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

impl crate::sealed::Sealed for SliceInput<'_> {}

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
    fn compression_format(&mut self) -> Result<Option<CompressionFormat>> {
        Ok(CompressionFormat::from_bytes(self.data))
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
    fn from_slice(data: &'a [u8]) -> Result<Self> {
        Self::from_input(SliceInput::new(data))
    }
}

impl<'a, F: FromInputData<'a, SliceInput<'a>>> FromSlice<'a> for F {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compression_format_detection_from_bytes() {
        // Uncompressed FASTA
        assert_eq!(CompressionFormat::from_bytes(b">seq1\nACGT\n"), None);
        // Uncompressed FASTQ
        assert_eq!(
            CompressionFormat::from_bytes(b"@read1\nACGT\n+\nIIII\n"),
            None
        );

        // Standard Gzip
        let gz_bytes = [0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff];
        assert_eq!(
            CompressionFormat::from_bytes(&gz_bytes),
            Some(CompressionFormat::Gzip)
        );

        // BGZF Gzip with BC extra subfield
        let bgz_bytes = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C',
            0x02, 0x00, 0xe2, 0x3d,
        ];
        assert_eq!(
            CompressionFormat::from_bytes(&bgz_bytes),
            Some(CompressionFormat::Bgzf)
        );

        // Standard Zstd
        let zstd_bytes = [0x28, 0xb5, 0x2f, 0xfd, 0x00, 0x00];
        assert_eq!(
            CompressionFormat::from_bytes(&zstd_bytes),
            Some(CompressionFormat::Zstd)
        );

        // PZSTD skippable frame
        let pzstd_bytes = [0x50, 0x2a, 0x4d, 0x18, 0x04, 0x00, 0x00, 0x00];
        assert_eq!(
            CompressionFormat::from_bytes(&pzstd_bytes),
            Some(CompressionFormat::Pzstd)
        );

        // Bzip2
        let bz2_bytes = *b"BZh9";
        assert_eq!(
            CompressionFormat::from_bytes(&bz2_bytes),
            Some(CompressionFormat::Bzip2)
        );

        // Xz
        let xz_bytes = [0xfd, b'7', b'z', b'X', b'Z', 0x00];
        assert_eq!(
            CompressionFormat::from_bytes(&xz_bytes),
            Some(CompressionFormat::Xz)
        );
    }

    #[test]
    fn test_compression_format_properties() {
        assert_eq!(CompressionFormat::Gzip.as_str(), "gzip");
        assert_eq!(CompressionFormat::Bgzf.as_str(), "bgzf");
        assert_eq!(CompressionFormat::Zstd.as_str(), "zstd");
        assert_eq!(CompressionFormat::Pzstd.as_str(), "pzstd");
        assert_eq!(CompressionFormat::Bzip2.as_str(), "bzip2");
        assert_eq!(CompressionFormat::Xz.as_str(), "xz");

        assert_eq!(format!("{}", CompressionFormat::Bgzf), "bgzf");
        assert_eq!(format!("{}", CompressionFormat::Pzstd), "pzstd");
    }
}
