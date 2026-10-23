//! Parallel decompression of independent compressed frames (BGZF blocks, zstd frames).

use super::{align_to_record, run_chunks};
use crate::config::*;
use crate::input::*;
use crate::parser::*;
use std::io::{self, Read};

#[derive(Clone, Copy)]
pub struct Frame {
    comp_offset: usize,
    comp_len: usize,
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

    // First decompressed byte of the frame, without decompressing all of it.
    fn first_byte(frame: &[u8]) -> io::Result<u8>;
}

#[cfg(feature = "sirocco")]
pub(super) struct BgzfDecoder(sirocco::inflate::Decoder);

#[cfg(feature = "sirocco")]
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

    fn first_byte(frame: &[u8]) -> io::Result<u8> {
        let mut prefix = [0u8; 64];
        let mut decoder = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);
        let first = decoder
            .decompress(frame, &mut prefix, false)
            .ok()
            .and_then(|p| p.output.first().copied());
        Ok(first.unwrap_or(0))
    }
}

#[cfg(feature = "sirocco")]
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
pub(super) struct ZstdDecoder(zstd::bulk::Decompressor<'static>);

#[cfg(feature = "zstd")]
impl FrameDecoder for ZstdDecoder {
    fn new() -> io::Result<Self> {
        Ok(Self(zstd::bulk::Decompressor::new()?))
    }

    fn decode(&mut self, frame: &[u8], uncomp_len: usize, out: &mut Vec<u8>) -> io::Result<()> {
        if uncomp_len > 0 {
            out.resize(uncomp_len, 0);
            let len = self.0.decompress_to_buffer(frame, out)?;
            out.truncate(len);
        } else {
            out.clear();
            zstd::stream::read::Decoder::new(frame)?.read_to_end(out)?;
        }
        Ok(())
    }

    fn first_byte(frame: &[u8]) -> io::Result<u8> {
        let mut first = [0u8; 1];
        zstd::stream::read::Decoder::new(frame)?.read_exact(&mut first)?;
        Ok(first[0])
    }
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
        match zstd::zstd_safe::find_frame_compressed_size(slice) {
            Ok(comp_len) if comp_len > 0 && comp_offset + comp_len <= compressed.len() => {
                let uncomp_len = zstd::zstd_safe::get_frame_content_size(slice)
                    .ok()
                    .flatten()
                    .unwrap_or(0) as usize;
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
            self.decoder.decode(
                frame.bytes(self.compressed),
                frame.uncomp_len,
                &mut self.buf,
            )?;
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
            if decoder
                .decode(frame.bytes(compressed), frame.uncomp_len, &mut data)
                .is_ok()
            {
                let cand = align_to_record(&data, 1, format);
                if cand < data.len() {
                    return (frame_idx, cand);
                }
            }
            frame_idx += 1;
        }
    }
    (frames.len(), 0)
}

pub(super) fn run_frames<D, const CONFIG: Config, State, Process, Reduce>(
    compressed: &[u8],
    frames: &[Frame],
    threads: usize,
    process: Process,
    reduce: Reduce,
) -> Result<State>
where
    D: FrameDecoder,
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> Result<State> + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let Some(first) = frames.first() else {
        return Ok(State::default());
    };
    let format = Format::try_from_byte(D::first_byte(first.bytes(compressed))?)?;
    let num_chunks = threads.min(frames.len());
    run_chunks(
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
        |start, end| {
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
            process(FastxParser::from_reader(reader)?)
        },
        reduce,
    )
}
