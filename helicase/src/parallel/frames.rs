//! Parallel decompression of independent compressed frames (BGZF blocks, zstd frames).

use super::Job;
use crate::config::*;
use crate::input::*;
use crate::parser::*;
use std::io::{self, Read};

#[cfg(all(test, feature = "zstd"))]
mod tests;

#[derive(Clone, Copy)]
pub struct Frame {
    comp_offset: usize,
    comp_len: usize,
    // Exact for BGZF, upper bound for zstd frames without content size.
    uncomp_len: usize,
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

#[cfg(feature = "gz")]
pub(super) fn index_bgzf_blocks(compressed: &[u8]) -> Vec<Frame> {
    let mut blocks = Vec::with_capacity(compressed.len() / 2048);
    let mut comp_offset = 0;
    while comp_offset + 18 <= compressed.len() {
        let header = &compressed[comp_offset..];
        if header[0] != 0x1f || header[1] != 0x8b || (header[3] & 0x04) == 0 {
            break;
        }
        if header[12] != b'B' || header[13] != b'C' {
            break;
        }
        let comp_len = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        if comp_offset + comp_len > compressed.len() {
            break;
        }
        let isize_bytes = &compressed[comp_offset + comp_len - 4..comp_offset + comp_len];
        let uncomp_len = u32::from_le_bytes(isize_bytes.try_into().unwrap()) as usize;
        if uncomp_len > 0 {
            blocks.push(Frame {
                comp_offset,
                comp_len,
                uncomp_len,
            });
        }
        comp_offset += comp_len;
    }
    blocks
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
pub fn index_zstd_frames(compressed: &[u8]) -> Vec<Frame> {
    let mut frames = Vec::with_capacity(compressed.len() / (1024 * 1024) + 1);
    let mut comp_offset = 0;
    while comp_offset + 4 <= compressed.len() {
        let slice = &compressed[comp_offset..];
        // Skippable frame: magic 0x184D2A50..=0x184D2A5F, then a u32 payload size.
        if slice.len() >= 8
            && slice[3] == 0x18
            && slice[2] == 0x4d
            && slice[1] == 0x2a
            && (slice[0] & 0xf0) == 0x50
        {
            let payload_len = u32::from_le_bytes(slice[4..8].try_into().unwrap()) as usize;
            comp_offset += 8 + payload_len;
            continue;
        }
        match zstd_frame_sizes(slice) {
            Some((comp_len, uncomp_len)) => {
                frames.push(Frame {
                    comp_offset,
                    comp_len,
                    uncomp_len,
                });
                comp_offset += comp_len;
            }
            _ => break,
        }
    }
    frames
}

// Decompressed data from position `start` to `end`, as (frame index, offset in the frame).
struct FrameRangeReader<'a, D> {
    compressed: &'a [u8],
    frames: &'a [Frame],
    start: (usize, usize),
    end: (usize, usize),
    next_frame: usize,
    decoder: D,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
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
            let bytes = frame.bytes(self.compressed);
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

// First record starting in frame `frame_idx` or after it.
fn find_record_start_in_frames<D: FrameDecoder>(
    compressed: &[u8],
    frames: &[Frame],
    mut frame_idx: usize,
    format: Format,
) -> (usize, usize) {
    let mut data = Vec::new();
    if let Ok(mut decoder) = D::new() {
        while frame_idx < frames.len() {
            let frame = frames[frame_idx];
            // The record start is usually close to the frame start, so decode a doubling prefix
            let mut len = 64 * 1024;
            while decoder
                .decode_prefix(frame.bytes(compressed), frame.uncomp_len, len, &mut data)
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
            frame_idx += 1;
        }
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
        compressed: &[u8],
        frames: &[Frame],
        threads: usize,
        custom_chunks: Option<usize>,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        let Some(first) = frames.first() else {
            return self.run_single(FastxParser::from_slice(&[]));
        };
        let mut first_bytes = Vec::new();
        D::new()?.decode_prefix(
            first.bytes(compressed),
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
                    compressed,
                    frames,
                    i * frames.len() / num_chunks,
                    format,
                )
            },
            |state, ticket, start, end| {
                let reader = FrameRangeReader {
                    compressed,
                    frames,
                    start,
                    end,
                    next_frame: start.0,
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
// frame starting `data`, by walking its block headers (RFC 8878).
#[cfg(feature = "zstd")]
fn zstd_frame_sizes(data: &[u8]) -> Option<(usize, usize)> {
    const MAX_BLOCK_LEN: usize = 128 * 1024;
    if data.get(..4)? != [0x28, 0xb5, 0x2f, 0xfd] {
        return None;
    }
    let descriptor = *data.get(4)?;
    let single_segment = descriptor & 0x20 != 0;
    let fcs_len = match descriptor >> 6 {
        0 => single_segment as usize,
        1 => 2,
        2 => 4,
        _ => 8,
    };
    let dict_len = [0, 1, 2, 4][(descriptor & 3) as usize];
    let fcs_start = 5 + !single_segment as usize + dict_len;
    let fcs = data.get(fcs_start..fcs_start + fcs_len)?;
    let content_len = fcs.iter().rev().fold(0, |acc, &b| (acc << 8) | b as usize)
        + if fcs_len == 2 { 256 } else { 0 };
    let mut pos = fcs_start + fcs_len;
    let mut num_blocks = 0;
    loop {
        let header = data.get(pos..pos + 3)?;
        let header = u32::from_le_bytes([header[0], header[1], header[2], 0]);
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
    (pos <= data.len()).then_some((pos, uncomp_len))
}
