//! Parallel decompression of independent compressed frames (BGZF blocks, zstd frames).

use super::{FramesSource, Job, with_thread_file};
use crate::config::*;
use crate::parser::*;
use std::io::{self, Read, Seek, SeekFrom};

#[cfg(all(test, feature = "zstd"))]
mod tests;

#[cfg(feature = "gz")]
const MAX_BGZF_BLOCK_LEN: usize = 1 << 16;
#[cfg(feature = "gz")]
const WINDOW_LEN: usize = 1 << 18;

#[derive(Clone, Copy)]
pub struct Frame {
    pub(super) comp_offset: usize,
    pub(super) comp_len: usize,
    // Exact for BGZF, upper bound for zstd frames without content size.
    pub(super) uncomp_len: usize,
}

impl Frame {
    fn bytes<'a>(&self, compressed: &'a [u8]) -> &'a [u8] {
        &compressed[self.comp_offset..self.comp_offset + self.comp_len]
    }
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

// Block at the start of `h`, possibly empty (the BGZF EOF marker is).
#[cfg(feature = "gz")]
#[inline(always)]
fn parse_bgzf_frame(h: &[u8], comp_offset: usize) -> Option<Frame> {
    if h.len() >= 18
        && h[0] == 0x1f
        && h[1] == 0x8b
        && (h[3] & 0x04) != 0
        && h[12] == b'B'
        && h[13] == b'C'
    {
        let comp_len = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
        if (18..=h.len()).contains(&comp_len) {
            let isize_bytes = &h[comp_len - 4..comp_len];
            let uncomp_len = u32::from_le_bytes(isize_bytes.try_into().unwrap()) as usize;
            return Some(Frame {
                comp_offset,
                comp_len,
                uncomp_len,
            });
        }
    }
    None
}

// `None` if the input ends with something else than BGZF blocks.
#[cfg(feature = "gz")]
pub(super) fn index_bgzf_blocks(compressed: &[u8]) -> Option<Vec<Frame>> {
    let mut blocks = Vec::new();
    let mut comp_offset = 0;
    while let Some(frame) = parse_bgzf_frame(&compressed[comp_offset..], comp_offset) {
        comp_offset += frame.comp_len;
        if frame.uncomp_len > 0 {
            blocks.push(frame);
        }
    }
    compressed[comp_offset..].is_empty().then_some(blocks)
}

#[cfg(feature = "gz")]
pub(super) fn index_bgzf_file(path: &std::path::Path) -> io::Result<Option<Vec<Frame>>> {
    let mut file = std::fs::File::open(path)?;
    let mut blocks = Vec::new();
    let mut comp_offset = 0;
    // The window holds the file bytes from `window_start`
    let mut window = Vec::with_capacity(WINDOW_LEN);
    let mut window_start = 0;
    loop {
        let mut rel = comp_offset - window_start;
        if window.len() - rel < MAX_BGZF_BLOCK_LEN {
            window.drain(..rel);
            window_start = comp_offset;
            rel = 0;
            (&mut file)
                .take((WINDOW_LEN - window.len()) as u64)
                .read_to_end(&mut window)?;
        }
        let Some(frame) = parse_bgzf_frame(&window[rel..], comp_offset) else {
            return Ok(window[rel..].is_empty().then_some(blocks));
        };
        comp_offset += frame.comp_len;
        if frame.uncomp_len > 0 {
            blocks.push(frame);
        }
    }
}

#[cfg(feature = "zstd")]
pub(super) struct ZstdDecoder(zstd::zstd_safe::DCtx<'static>);

#[cfg(feature = "zstd")]
impl FrameDecoder for ZstdDecoder {
    fn new() -> io::Result<Self> {
        let context = zstd::zstd_safe::DCtx::try_create();
        context
            .map(Self)
            .ok_or_else(|| io::Error::other("cannot create zstd context"))
    }

    fn decode(&mut self, frame: &[u8], uncomp_len: usize, out: &mut Vec<u8>) -> io::Result<()> {
        out.resize(uncomp_len, 0);
        let len = self
            .0
            .decompress(out.as_mut_slice(), frame)
            .map_err(zstd_error)?;
        out.truncate(len);
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
        use zstd::zstd_safe::{InBuffer, OutBuffer, ResetDirective};
        let len = len.min(uncomp_len);
        out.resize(len, 0);
        self.0
            .reset(ResetDirective::SessionOnly)
            .map_err(zstd_error)?;
        let mut input = InBuffer::around(frame);
        let mut output = OutBuffer::around(out.as_mut_slice());
        while output.pos() < len {
            let remaining = self
                .0
                .decompress_stream(&mut output, &mut input)
                .map_err(zstd_error)?;
            if remaining == 0 {
                break;
            }
            if input.pos() == frame.len() {
                return Err(io::ErrorKind::UnexpectedEof.into());
            }
        }
        let decoded = output.pos();
        out.truncate(decoded);
        Ok(())
    }
}

#[cfg(feature = "zstd")]
fn zstd_error(code: usize) -> io::Error {
    io::Error::other(zstd::zstd_safe::get_error_name(code))
}

// Frames of the file, without the skippable frames written by pzstd.
#[cfg(feature = "zstd")]
pub(super) fn index_zstd_frames(reader: &mut (impl Read + Seek)) -> Vec<Frame> {
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
    frames
}

// Decompressed data from position `start` to `end`, as (frame index, offset in the frame).
struct FrameRangeReader<'a, D> {
    source: FramesSource<'a>,
    frames: &'a [Frame],
    start: (usize, usize),
    end: (usize, usize),
    next_frame: usize,
    decoder: D,
    comp_buf: Vec<u8>,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
}

impl<'a, D: FrameDecoder> FrameRangeReader<'a, D> {
    const fn new(
        source: FramesSource<'a>,
        frames: &'a [Frame],
        start: (usize, usize),
        end: (usize, usize),
        decoder: D,
    ) -> Self {
        Self {
            source,
            frames,
            start,
            end,
            next_frame: start.0,
            decoder,
            comp_buf: Vec::new(),
            buf: Vec::new(),
            pos: 0,
            len: 0,
        }
    }
}

impl<D: FrameDecoder> Read for FrameRangeReader<'_, D> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.len {
            let idx = self.next_frame;
            if idx > self.end.0
                || (idx == self.end.0 && self.end.1 == 0)
                || idx >= self.frames.len()
            {
                return Ok(0);
            }
            let frame = self.frames[idx];
            let bytes = frame_bytes(self.source, &frame, &mut self.comp_buf)?;
            if idx == self.end.0 {
                // Only the bytes before the end of the range are needed
                let len = self.end.1;
                self.decoder
                    .decode_prefix(bytes, frame.uncomp_len, len, &mut self.buf)?;
            } else {
                self.decoder
                    .decode(bytes, frame.uncomp_len, &mut self.buf)?;
            }
            self.pos = if idx == self.start.0 { self.start.1 } else { 0 };
            self.len = if idx == self.end.0 {
                self.end.1
            } else {
                self.buf.len()
            };
            self.next_frame += 1;
        }
        let n = (self.len - self.pos).min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

// The compressed bytes of `frame`, read into `buf` if the source is a file.
fn frame_bytes<'b>(
    source: FramesSource<'b>,
    frame: &Frame,
    buf: &'b mut Vec<u8>,
) -> io::Result<&'b [u8]> {
    match source {
        FramesSource::Slice(slice) => Ok(frame.bytes(slice)),
        FramesSource::File(path) => {
            buf.resize(frame.comp_len, 0);
            with_thread_file(path, |file| {
                file.seek(SeekFrom::Start(frame.comp_offset as u64))?;
                Ok(file.read_exact(buf)?)
            })?;
            Ok(buf)
        }
    }
}

// First record starting in frame `frame_idx` or after it.
fn find_record_start_in_frames<D: FrameDecoder>(
    source: FramesSource<'_>,
    frames: &[Frame],
    mut frame_idx: usize,
    format: Format,
) -> (usize, usize) {
    let (mut data, mut comp_buf) = (Vec::new(), Vec::new());
    let Ok(mut decoder) = D::new() else {
        return (frames.len(), 0);
    };
    while frame_idx < frames.len() {
        let frame = frames[frame_idx];
        if let Ok(bytes) = frame_bytes(source, &frame, &mut comp_buf) {
            // The record start is usually close to the frame start, so decode a doubling prefix
            let mut len = 64 * 1024;
            while decoder
                .decode_prefix(bytes, frame.uncomp_len, len, &mut data)
                .is_ok()
            {
                if let Some(cand) = format.find_record_start(&data, 0) {
                    return (frame_idx, cand);
                }
                if data.len() < len {
                    break;
                }
                len *= 2;
            }
        }
        frame_idx += 1;
    }
    (frames.len(), 0)
}

impl<const CONFIG: Config, T, State, Init, Process, Finish> Job<CONFIG, Init, Process, Finish>
where
    State: Send,
    Init: Fn(usize) -> T + Sync,
    Process: Fn(&mut T, usize, FastxParser<'_, CONFIG>) -> Result<()> + Sync,
    Finish: Fn(T) -> Result<State> + Sync,
{
    pub(super) fn run_frames<D: FrameDecoder>(
        &self,
        source: FramesSource<'_>,
        frames: &[Frame],
        threads: usize,
        custom_chunks: Option<usize>,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        let Some(first) = frames.first() else {
            return self.run_single(FastxParser::from_slice(&[]));
        };
        let (mut first_bytes, mut comp_buf) = (Vec::new(), Vec::new());
        D::new()?.decode_prefix(
            frame_bytes(source, first, &mut comp_buf)?,
            first.uncomp_len,
            1,
            &mut first_bytes,
        )?;
        let format = Format::try_from_byte(first_bytes.first().copied().unwrap_or(0))?;
        let uncomp_len: usize = frames.iter().map(|frame| frame.uncomp_len).sum();
        let num_chunks = custom_chunks
            .unwrap_or(threads)
            .max(uncomp_len.div_ceil(max_chunk_len))
            .min(frames.len());
        self.run_chunks(
            num_chunks,
            threads,
            ((0, 0), (frames.len(), 0)),
            |i| {
                find_record_start_in_frames::<D>(
                    source,
                    frames,
                    i * frames.len() / num_chunks,
                    format,
                )
            },
            |state, ticket, start, end| {
                let reader = FrameRangeReader::new(source, frames, start, end, D::new()?);
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
