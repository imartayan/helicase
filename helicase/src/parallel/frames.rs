//! Parallel decompression of independent compressed frames (BGZF blocks, zstd frames).

use super::{FramesSource, Job, with_thread_file};
use crate::config::*;
use crate::parser::*;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

#[cfg(all(test, feature = "zstd"))]
mod tests;

#[cfg(feature = "gz")]
const MAX_BGZF_BLOCK_LEN: usize = 1 << 16;
#[cfg(feature = "gz")]
const WINDOW_LEN: usize = 1 << 18;
// A block, then the header of the next one
#[cfg(feature = "gz")]
const RESYNC_LEN: usize = 2 * MAX_BGZF_BLOCK_LEN + 18;
// The BGZF EOF marker, an empty block
#[cfg(feature = "gz")]
const EOF_BLOCK_LEN: usize = 28;

#[derive(Clone, Copy)]
pub struct Frame {
    pub(super) comp_offset: usize,
    pub(super) comp_len: usize,
    // Exact for BGZF, upper bound for zstd frames without content size.
    pub(super) uncomp_len: usize,
}

pub(super) trait FrameSource: Sized + Send {
    // Next non-empty frame starting before the limit, with its compressed bytes.
    fn next_frame(&mut self) -> io::Result<Option<(Frame, &[u8])>>;
}

pub(super) trait FrameDecoder: Sized + Send {
    fn new() -> io::Result<Self>;

    // Replaces the content of `out` with the decompressed frame.
    fn decode(&mut self, frame: &[u8], uncomp_len: usize, out: &mut Vec<u8>) -> io::Result<()>;

    // Like `decode`, but may stop after the first `len` bytes.
    fn decode_prefix(
        &mut self,
        frame: &[u8],
        uncomp_len: usize,
        _len: usize,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        self.decode(frame, uncomp_len, out)
    }
}

#[cfg(feature = "gz")]
pub(super) struct BgzfDecoder(sirocco::inflate::Decoder);

#[cfg(feature = "gz")]
impl FrameDecoder for BgzfDecoder {
    fn new() -> io::Result<Self> {
        Ok(Self(sirocco::inflate::Decoder::new(sirocco::Format::Gzip)))
    }

    fn decode(&mut self, frame: &[u8], uncomp_len: usize, out: &mut Vec<u8>) -> io::Result<()> {
        out.resize(uncomp_len, 0);
        let len = self.0.decompress(frame, out, true).map(|p| p.output.len());
        self.0.reset();
        out.truncate(
            len.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:?}")))?,
        );
        Ok(())
    }
}

// Compressed length of the BGZF block starting at `h`, from its header.
#[cfg(feature = "gz")]
fn bgzf_block_len(h: &[u8]) -> Option<usize> {
    (h.len() >= 18 && h[..4] == [0x1f, 0x8b, 8, 4] && h[10..16] == [6, 0, b'B', b'C', 2, 0])
        .then(|| u16::from_le_bytes([h[16], h[17]]) as usize + 1)
}

// Block at the start of `h`, possibly empty (the BGZF EOF marker is).
#[cfg(feature = "gz")]
#[inline(always)]
fn parse_bgzf_frame(h: &[u8], comp_offset: usize) -> Option<Frame> {
    let comp_len = bgzf_block_len(h)?;
    if comp_len > h.len() {
        return None;
    }
    let isize_bytes = &h[comp_len - 4..comp_len];
    let uncomp_len = u32::from_le_bytes(isize_bytes.try_into().unwrap()) as usize;
    Some(Frame {
        comp_offset,
        comp_len,
        uncomp_len,
    })
}

// Appends up to `len` bytes of the file from `offset` to `buf`.
fn read_at(path: &Path, buf: &mut Vec<u8>, offset: usize, len: usize) -> Result<()> {
    with_thread_file(path, |file| {
        file.seek(SeekFrom::Start(offset as u64))?;
        file.take(len as u64).read_to_end(buf)?;
        Ok(())
    })
}

// Length of the input if it ends with the BGZF EOF marker.
#[cfg(feature = "gz")]
pub(super) fn bgzf_input_len(source: FramesSource<'_>) -> Result<Option<usize>> {
    let is_eof_block = |tail: &[u8]| {
        parse_bgzf_frame(tail, 0).is_some_and(|b| b.comp_len == EOF_BLOCK_LEN && b.uncomp_len == 0)
    };
    match source {
        FramesSource::Slice(slice) => {
            let tail = &slice[slice.len().saturating_sub(EOF_BLOCK_LEN)..];
            Ok(is_eof_block(tail).then_some(slice.len()))
        }
        FramesSource::File(path) => {
            let len = with_thread_file(path, |file| Ok(file.metadata()?.len() as usize))?;
            let mut tail = Vec::new();
            read_at(
                path,
                &mut tail,
                len.saturating_sub(EOF_BLOCK_LEN),
                EOF_BLOCK_LEN,
            )?;
            Ok(is_eof_block(&tail).then_some(len))
        }
    }
}

// BGZF blocks from `pos` to `limit`, parsed on the fly. A file is read through a window.
#[cfg(feature = "gz")]
pub(super) struct BgzfBlocks<'a> {
    source: FramesSource<'a>,
    pos: usize,
    limit: usize,
    window: Vec<u8>,
    window_start: usize,
}

#[cfg(feature = "gz")]
impl<'a> BgzfBlocks<'a> {
    // Starts at the first block header at or after `pos`, accepted if the next block (or the
    // end of the input) follows it. `pos` can be anywhere.
    pub(super) fn new(source: FramesSource<'a>, pos: usize, limit: usize) -> io::Result<Self> {
        let mut blocks = Self {
            source,
            pos,
            limit,
            window: Vec::new(),
            window_start: 0,
        };
        blocks.fill(pos, RESYNC_LEN)?;
        let h = blocks.bytes_from(pos);
        // The input starts with a block: a bad one is reported when read, never skipped
        let found = (0..h.len().min(MAX_BGZF_BLOCK_LEN)).find(|&i| {
            pos + i == 0
                || bgzf_block_len(&h[i..]).is_some_and(|len| {
                    h.get(i + len..)
                        .is_some_and(|next| next.is_empty() || bgzf_block_len(next).is_some())
                })
        });
        // The previous chunk reports the invalid data it will meet
        blocks.pos = found.map_or(limit, |i| pos + i);
        Ok(blocks)
    }

    // Makes the window hold at least `len` bytes from `pos`, unless the file ends first.
    fn fill(&mut self, pos: usize, len: usize) -> Result<()> {
        let FramesSource::File(path) = self.source else {
            return Ok(());
        };
        let window_end = self.window_start + self.window.len();
        if pos >= self.window_start && pos + len <= window_end {
            return Ok(());
        }
        if pos >= self.window_start && pos <= window_end {
            self.window.drain(..pos - self.window_start);
        } else {
            self.window.clear();
        }
        self.window_start = pos;
        let kept = self.window.len();
        read_at(path, &mut self.window, pos + kept, WINDOW_LEN - kept)
    }

    fn bytes_from(&self, pos: usize) -> &[u8] {
        match self.source {
            FramesSource::Slice(slice) => &slice[pos..],
            FramesSource::File(_) => &self.window[pos - self.window_start..],
        }
    }
}

#[cfg(feature = "gz")]
impl FrameSource for BgzfBlocks<'_> {
    fn next_frame(&mut self) -> io::Result<Option<(Frame, &[u8])>> {
        loop {
            let pos = self.pos;
            if pos >= self.limit {
                return Ok(None);
            }
            self.fill(pos, MAX_BGZF_BLOCK_LEN)?;
            let Some(frame) = parse_bgzf_frame(self.bytes_from(pos), pos) else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid BGZF block",
                ));
            };
            self.pos += frame.comp_len;
            if frame.uncomp_len > 0 {
                return Ok(Some((frame, &self.bytes_from(pos)[..frame.comp_len])));
            }
        }
    }
}

#[cfg(feature = "zstd")]
pub(super) struct ZstdDecoder(zefiro::Decoder);

#[cfg(feature = "zstd")]
impl FrameDecoder for ZstdDecoder {
    fn new() -> io::Result<Self> {
        let mut decoder = zefiro::Decoder::new();
        decoder.set_verify_checksum(false);
        Ok(Self(decoder))
    }

    fn decode(&mut self, frame: &[u8], uncomp_len: usize, out: &mut Vec<u8>) -> io::Result<()> {
        if out.capacity() < uncomp_len {
            out.reserve(uncomp_len - out.capacity());
        }
        self.0.decompress_into_vec(frame, out)?;
        Ok(())
    }

    // Streams the frame block by block (at most 128 KiB each) until `len` bytes are out.
    fn decode_prefix(
        &mut self,
        frame: &[u8],
        uncomp_len: usize,
        len: usize,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        let len = len.min(uncomp_len);
        out.resize(len, 0);
        self.0.reset();
        let mut in_pos = 0;
        let mut out_pos = 0;
        while out_pos < len {
            let progress = self
                .0
                .decompress_stream(&frame[in_pos..], &mut out[out_pos..])?;
            in_pos += progress.consumed();
            out_pos += progress.written();
            if progress.frame_done() || (progress.consumed() == 0 && progress.written() == 0) {
                break;
            }
            if in_pos == frame.len() && out_pos < len {
                return Err(io::ErrorKind::UnexpectedEof.into());
            }
        }
        self.0.reset();
        out.truncate(out_pos);
        Ok(())
    }
}

// Frames of the file, without the skippable frames written by pzstd.
// `None` if the input ends with something else than zstd frames.
#[cfg(feature = "zstd")]
pub(super) fn index_zstd_frames(reader: &mut (impl Read + Seek)) -> Option<Vec<Frame>> {
    let mut frames = Vec::new();
    let mut comp_offset = 0;
    // `zstd_frame_sizes` leaves the reader at the end of the frame
    while let Some((comp_len, uncomp_len)) = zstd_frame_sizes(reader) {
        if uncomp_len > 0 {
            frames.push(Frame {
                comp_offset,
                comp_len,
                uncomp_len,
            });
        }
        comp_offset += comp_len;
    }
    (reader.seek(SeekFrom::End(0)).ok()? == comp_offset as u64).then_some(frames)
}

// The indexed frames from `limit` on, the bytes of a frame are read when it is returned.
#[cfg(feature = "zstd")]
pub(super) struct ZstdFrames<'a> {
    source: FramesSource<'a>,
    frames: &'a [Frame],
    limit: usize,
    buf: Vec<u8>,
}

#[cfg(feature = "zstd")]
impl<'a> ZstdFrames<'a> {
    pub(super) fn new(
        source: FramesSource<'a>,
        frames: &'a [Frame],
        pos: usize,
        limit: usize,
    ) -> io::Result<Self> {
        Ok(Self {
            source,
            // The frame whose middle is the first at or after `pos`, to round to the nearest
            frames: &frames[frames.partition_point(|f| f.comp_offset + f.comp_len / 2 < pos)..],
            limit,
            buf: Vec::new(),
        })
    }
}

#[cfg(feature = "zstd")]
impl FrameSource for ZstdFrames<'_> {
    fn next_frame(&mut self) -> io::Result<Option<(Frame, &[u8])>> {
        let Some((&frame, rest)) = self
            .frames
            .split_first()
            .filter(|(frame, _)| frame.comp_offset < self.limit)
        else {
            return Ok(None);
        };
        self.frames = rest;
        let range = frame.comp_offset..frame.comp_offset + frame.comp_len;
        let bytes = match self.source {
            FramesSource::Slice(slice) => &slice[range],
            FramesSource::File(path) => {
                self.buf.clear();
                read_at(path, &mut self.buf, frame.comp_offset, frame.comp_len)?;
                &self.buf
            }
        };
        Ok(Some((frame, bytes)))
    }
}

// Decompressed data from position `start` to `end`, as (frame offset, offset in the frame).
struct FrameRangeReader<S, D> {
    frames: S,
    start: (usize, usize),
    end: (usize, usize),
    decoder: D,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
}

impl<S: FrameSource, D: FrameDecoder> Read for FrameRangeReader<S, D> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.len {
            let Some((frame, bytes)) = self.frames.next_frame()? else {
                return Ok(0);
            };
            let is_last = frame.comp_offset == self.end.0;
            if is_last {
                // Only the bytes before the end of the range are needed
                let len = self.end.1;
                self.decoder
                    .decode_prefix(bytes, frame.uncomp_len, len, &mut self.buf)?;
            } else {
                self.decoder
                    .decode(bytes, frame.uncomp_len, &mut self.buf)?;
            }
            self.pos = if frame.comp_offset == self.start.0 {
                self.start.1
            } else {
                0
            };
            self.len = if is_last { self.end.1 } else { self.buf.len() };
        }
        let n = (self.len - self.pos).min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

// First record starting in the frames, as (frame offset, offset in the frame).
fn find_record_start_in_frames<D: FrameDecoder>(
    mut frames: impl FrameSource,
    format: Format,
) -> Option<(usize, usize)> {
    let mut data = Vec::new();
    let mut decoder = D::new().ok()?;
    while let Some((frame, bytes)) = frames.next_frame().ok()? {
        // The record start is usually close to the frame start, so decode a doubling prefix
        let mut len = 64 * 1024;
        while decoder
            .decode_prefix(bytes, frame.uncomp_len, len, &mut data)
            .is_ok()
        {
            if let Some(cand) = format.find_record_start(&data, 0) {
                return Some((frame.comp_offset, cand));
            }
            if data.len() < len {
                break;
            }
            len *= 2;
        }
    }
    None
}

impl<const CONFIG: Config, T, State, Init, Process, Finish> Job<CONFIG, Init, Process, Finish>
where
    State: Send,
    Init: Fn(usize) -> T + Sync,
    Process: Fn(&mut T, usize, FastxParser<'_, CONFIG>) -> Result<()> + Sync,
    Finish: Fn(T) -> Result<State> + Sync,
{
    // `open(pos, limit)` gives the frames starting in `pos..limit` of the `len` compressed bytes.
    pub(super) fn run_frames<D: FrameDecoder, S: FrameSource>(
        &self,
        len: usize,
        open: impl Fn(usize, usize) -> io::Result<S> + Sync,
        threads: usize,
        custom_chunks: Option<usize>,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        let mut frames = open(0, len)?;
        let Some((first, first_bytes)) = frames.next_frame()? else {
            return self.run_single(FastxParser::from_slice(&[]));
        };
        let mut head = Vec::new();
        D::new()?.decode_prefix(first_bytes, first.uncomp_len, 1, &mut head)?;
        let format = Format::try_from_byte(head.first().copied().unwrap_or(0))?;
        // Estimates from the first frame
        let uncomp_len = len / first.comp_len * first.uncomp_len;
        let num_chunks = custom_chunks
            .unwrap_or(threads)
            .max(uncomp_len.div_ceil(max_chunk_len))
            .min(len.div_ceil(first.comp_len));
        self.run_chunks(
            num_chunks,
            threads,
            ((0, 0), (len, 0)),
            |i| {
                open(i * len / num_chunks, len)
                    .ok()
                    .and_then(|frames| find_record_start_in_frames::<D>(frames, format))
                    .unwrap_or((len, 0))
            },
            |state, ticket, start, end| {
                let reader = FrameRangeReader {
                    frames: open(start.0, end.0 + (end.1 > 0) as usize)?,
                    start,
                    end,
                    decoder: D::new()?,
                    buf: Vec::new(),
                    pos: 0,
                    len: 0,
                };
                (self.process)(state, ticket, FastxParser::from_reader(reader)?)
            },
        )
    }
}

// Compressed size and decompressed size (or an upper bound if the header has none) of the zstd
// frame at the current position, by walking its block headers (RFC 8878). Skippable frames have
// a decompressed size of 0. On success, the reader is left at the end of the frame.
#[cfg(feature = "zstd")]
pub fn zstd_frame_sizes(reader: &mut (impl Read + Seek)) -> Option<(usize, usize)> {
    const MAX_BLOCK_LEN: usize = 128 * 1024;
    let start = reader.stream_position().ok()?;
    let mut read_at = |pos: usize, buf: &mut [u8]| {
        reader.seek(SeekFrom::Start(start + pos as u64)).ok()?;
        reader.read_exact(buf).ok()
    };
    let mut head = [0u8; 18];
    read_at(0, &mut head[..5])?;
    let (len, uncomp_len) = if head[1..4] == [0x2a, 0x4d, 0x18] && head[0] & 0xf0 == 0x50 {
        read_at(4, &mut head[..4])?;
        (
            8 + u32::from_le_bytes(head[..4].try_into().unwrap()) as usize,
            0,
        )
    } else if head[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
        let descriptor = head[4];
        let single_segment = descriptor & 0x20 != 0;
        let fcs_len = match descriptor >> 6 {
            0 => single_segment as usize,
            1 => 2,
            2 => 4,
            _ => 8,
        };
        let dict_len = [0, 1, 2, 4][(descriptor & 3) as usize];
        let fcs_start = 5 + !single_segment as usize + dict_len;
        let mut pos = fcs_start + fcs_len;
        read_at(5, &mut head[5..pos])?;
        let content_len = head[fcs_start..pos]
            .iter()
            .rev()
            .fold(0, |acc, &b| (acc << 8) | b as usize)
            + if fcs_len == 2 { 256 } else { 0 };
        let mut num_blocks = 0;
        loop {
            read_at(pos, &mut head[..3])?;
            let header = u32::from_le_bytes([head[0], head[1], head[2], 0]);
            pos += 3 + match (header >> 1) & 3 {
                1 => 1,
                3 => return None,
                _ => (header >> 3) as usize,
            };
            num_blocks += 1;
            if header & 1 != 0 {
                break;
            }
        }
        pos += if descriptor & 4 != 0 { 4 } else { 0 };
        let uncomp_len = if fcs_len > 0 {
            content_len
        } else {
            num_blocks * MAX_BLOCK_LEN
        };
        (pos, uncomp_len)
    } else {
        return None;
    };
    // Reading the last byte checks that the frame is complete
    read_at(len - 1, &mut head[..1])?;
    Some((len, uncomp_len))
}
