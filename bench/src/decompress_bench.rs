//! Benchmark evaluating Sirocco (pure-Rust v0.2.0) vs flate2 (zlib-rs baseline):
//! - Evaluates Standard Single-Member Gzip (150bp Illumina & ~15kb PacBio)
//! - Evaluates Multi-Member BGZF Gzip with table-preserving Decoder reuse & parallel Rayon block chunking
//! - Evaluates End-to-End FASTX Parsing with Helicase SIMD

use std::fs::File;
use std::hint::black_box;
#[cfg(feature = "sirocco")]
use std::io::BufRead;
use std::io::{BufReader, Read};
use std::path::Path;
use std::time::Instant;

use colored::Colorize;
use helicase::HelicaseParser;
use helicase::config::advanced::DEFAULT_CONFIG;
use helicase::input::{FromReader, FromSlice};
use helicase::parser::FastxParser;

// ============================================================================
// Sirocco Streaming Multi-Member Reader & BGZF Strategy
// ============================================================================

#[cfg(feature = "sirocco")]
pub struct SiroccoGzipReader<R: BufRead> {
    current: Option<sirocco::inflate::Reader<R>>,
}

#[cfg(feature = "sirocco")]
impl<R: BufRead> SiroccoGzipReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            current: Some(sirocco::inflate::Reader::new(sirocco::Format::Gzip, reader)),
        }
    }
}

#[cfg(feature = "sirocco")]
impl<R: BufRead> Read for SiroccoGzipReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        loop {
            let Some(mut r) = self.current.take() else {
                return Ok(0);
            };
            let n = r.read(buf)?;
            if n > 0 {
                self.current = Some(r);
                return Ok(n);
            }
            // Stream ended (n == 0). Check if inner reader has another gzip member (BGZF/multi-member)
            let mut inner = r.into_inner();
            if inner.fill_buf()?.is_empty() {
                return Ok(0);
            }
            self.current = Some(sirocco::inflate::Reader::new(sirocco::Format::Gzip, inner));
        }
    }
}

#[cfg(feature = "sirocco")]
pub const fn is_bgzf(bytes: &[u8]) -> bool {
    bytes.len() >= 18
        && bytes[0] == 0x1f
        && bytes[1] == 0x8b
        && bytes[2] == 0x08
        && (bytes[3] & 0x04) != 0 // FLG.FEXTRA
        && bytes[12] == b'B'
        && bytes[13] == b'C'
        && bytes[14] == 0x02
        && bytes[15] == 0x00
}

#[cfg(feature = "sirocco")]
#[derive(Clone, Copy)]
pub struct BgzfBlockIndex {
    pub comp_offset: usize,
    pub comp_len: usize,
    pub uncomp_offset: usize,
    pub uncomp_len: usize,
}

#[cfg(feature = "sirocco")]
fn index_bgzf_blocks(compressed: &[u8]) -> (Vec<BgzfBlockIndex>, usize) {
    let mut blocks = Vec::with_capacity(compressed.len() / 2048);
    let mut comp_offset = 0;
    let mut uncomp_offset = 0;

    while comp_offset + 18 <= compressed.len() {
        let block_header = &compressed[comp_offset..];
        if block_header[0] != 0x1f || block_header[1] != 0x8b || (block_header[3] & 0x04) == 0 {
            break;
        }
        if block_header[12] != b'B' || block_header[13] != b'C' {
            break;
        }
        let bsize = u16::from_le_bytes([block_header[16], block_header[17]]) as usize;
        let comp_len = bsize + 1;
        if comp_offset + comp_len > compressed.len() {
            break;
        }
        let block = &compressed[comp_offset..comp_offset + comp_len];
        let uncomp_len = u32::from_le_bytes([
            block[comp_len - 4],
            block[comp_len - 3],
            block[comp_len - 2],
            block[comp_len - 1],
        ]) as usize;

        if uncomp_len > 0 {
            blocks.push(BgzfBlockIndex {
                comp_offset,
                comp_len,
                uncomp_offset,
                uncomp_len,
            });
            uncomp_offset += uncomp_len;
        }
        comp_offset += comp_len;
    }
    (blocks, uncomp_offset)
}

#[cfg(feature = "sirocco")]
pub fn decompress_bgzf_sirocco_1t(compressed: &[u8]) -> Vec<u8> {
    let (blocks, total_uncomp) = index_bgzf_blocks(compressed);
    let mut out = vec![0u8; total_uncomp];
    let mut decoder = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);

    for b in &blocks {
        let block = &compressed[b.comp_offset..b.comp_offset + b.comp_len];
        let out_slice = &mut out[b.uncomp_offset..b.uncomp_offset + b.uncomp_len];
        let p = decoder
            .decompress(block, out_slice, true)
            .expect("sirocco block decode failed");
        assert_eq!(p.output.len(), b.uncomp_len);
        decoder.reset();
    }
    out
}

#[cfg(feature = "sirocco")]
pub fn decompress_bgzf_sirocco_parallel(compressed: &[u8], num_threads: usize) -> Vec<u8> {
    use rayon::prelude::*;
    use std::cell::RefCell;

    let (blocks, total_uncomp) = index_bgzf_blocks(compressed);
    let mut out = vec![0u8; total_uncomp];
    let out_ptr = out.as_mut_ptr() as usize;

    thread_local! {
        static DECODER: RefCell<sirocco::inflate::Decoder> = const { RefCell::new(sirocco::inflate::Decoder::new(sirocco::Format::Gzip)) };
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    pool.install(|| {
        blocks.par_chunks(16).for_each(|chunk| {
            DECODER.with(|dec| {
                let mut decoder = dec.borrow_mut();
                for b in chunk {
                    let block = &compressed[b.comp_offset..b.comp_offset + b.comp_len];
                    let out_slice = unsafe {
                        core::slice::from_raw_parts_mut(
                            (out_ptr + b.uncomp_offset) as *mut u8,
                            b.uncomp_len,
                        )
                    };
                    let p = decoder
                        .decompress(block, out_slice, true)
                        .expect("sirocco block decode failed");
                    assert_eq!(p.output.len(), b.uncomp_len);
                    decoder.reset();
                }
            });
        });
    });

    out
}

#[cfg(feature = "sirocco")]
pub struct BgzfBlockRangeReader<'a> {
    compressed: &'a [u8],
    blocks: &'a [BgzfBlockIndex],
    current_block: usize,
    start_block: usize,
    end_block: usize,
    start_skip: usize,
    end_take: usize,

    decoder: sirocco::inflate::Decoder,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
}

#[cfg(feature = "sirocco")]
impl<'a> BgzfBlockRangeReader<'a> {
    pub fn new(
        compressed: &'a [u8],
        blocks: &'a [BgzfBlockIndex],
        start_block: usize,
        end_block: usize,
        start_skip: usize,
        end_take: usize,
    ) -> Self {
        Self {
            compressed,
            blocks,
            current_block: start_block,
            start_block,
            end_block,
            start_skip,
            end_take,
            decoder: sirocco::inflate::Decoder::new(sirocco::Format::Gzip),
            buf: Vec::with_capacity(65536),
            pos: 0,
            len: 0,
        }
    }
}

#[cfg(feature = "sirocco")]
impl<'a> Read for BgzfBlockRangeReader<'a> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        while self.pos >= self.len {
            if self.current_block > self.end_block
                || (self.current_block == self.end_block && self.end_take == 0)
                || self.current_block >= self.blocks.len()
            {
                return Ok(0);
            }
            let b = &self.blocks[self.current_block];
            let block_bytes = &self.compressed[b.comp_offset..b.comp_offset + b.comp_len];
            self.buf.resize(b.uncomp_len, 0);
            let p = self
                .decoder
                .decompress(block_bytes, &mut self.buf, true)
                .map_err(|e| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{e:?}"))
                })?;
            self.decoder.reset();

            let mut start = 0;
            let mut end = p.output.len();

            if self.current_block == self.start_block {
                start = self.start_skip;
            }
            if self.current_block == self.end_block {
                end = self.end_take;
            }

            self.pos = start;
            self.len = end;
            self.current_block += 1;
        }

        let n = (self.len - self.pos).min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

#[cfg(feature = "sirocco")]
#[inline(always)]
fn is_valid_fastq_record_start(data: &[u8], cand: usize) -> bool {
    let rem = &data[cand..];
    let Some(header_nl) = memchr::memchr(b'\n', rem) else {
        return false;
    };
    let seq_start = cand + header_nl + 1;
    if seq_start >= data.len() || data[seq_start] == b'@' || data[seq_start] == b'+' {
        return false;
    }
    let Some(seq_nl) = memchr::memchr(b'\n', &data[seq_start..]) else {
        return false;
    };
    let seq_end = seq_start + seq_nl;
    let plus_start = seq_end + 1;
    if plus_start >= data.len() || data[plus_start] != b'+' {
        return false;
    }
    let Some(plus_nl) = memchr::memchr(b'\n', &data[plus_start..]) else {
        return false;
    };
    let qual_start = plus_start + plus_nl + 1;
    if qual_start >= data.len() {
        return false;
    }
    let mut seq_len = seq_end - seq_start;
    if seq_len > 0 && data[seq_start + seq_len - 1] == b'\r' {
        seq_len -= 1;
    }
    let qual_target = qual_start + seq_len;
    if qual_target >= data.len() {
        return false;
    }
    let (is_match, next_rec) = if data[qual_target] == b'\n' {
        (true, qual_target + 1)
    } else if qual_target + 1 < data.len()
        && data[qual_target] == b'\r'
        && data[qual_target + 1] == b'\n'
    {
        (true, qual_target + 2)
    } else {
        (false, 0)
    };
    if !is_match {
        return false;
    }
    next_rec >= data.len() || data[next_rec] == b'@'
}

#[cfg(feature = "sirocco")]
fn find_fasta_record_start(data: &[u8], approx_pos: usize) -> usize {
    if approx_pos == 0 || approx_pos >= data.len() {
        return approx_pos.min(data.len());
    }
    if data[approx_pos] == b'>' && data[approx_pos - 1] == b'\n' {
        return approx_pos;
    }
    let mut pos = approx_pos - 1;
    while let Some(newline) = memchr::memchr(b'\n', &data[pos..]) {
        let cand = pos + newline + 1;
        if cand >= data.len() {
            return data.len();
        }
        if data[cand] == b'>' {
            return cand;
        }
        pos = cand;
    }
    data.len()
}

#[cfg(feature = "sirocco")]
fn find_fastq_record_start(data: &[u8], approx_pos: usize) -> usize {
    if approx_pos == 0 || approx_pos >= data.len() {
        return approx_pos.min(data.len());
    }
    let mut pos = approx_pos.saturating_sub(1);
    while let Some(newline) = memchr::memchr(b'\n', &data[pos..]) {
        let cand = pos + newline + 1;
        if cand >= data.len() {
            return data.len();
        }
        if data[cand] == b'@' && is_valid_fastq_record_start(data, cand) {
            return cand;
        }
        pos = cand;
    }
    data.len()
}

#[cfg(feature = "sirocco")]
#[inline(always)]
fn find_record_start(data: &[u8], approx_pos: usize, is_fastq: bool) -> usize {
    if is_fastq {
        find_fastq_record_start(data, approx_pos)
    } else {
        find_fasta_record_start(data, approx_pos)
    }
}

#[cfg(feature = "sirocco")]
fn find_record_start_in_bgzf(
    compressed: &[u8],
    blocks: &[BgzfBlockIndex],
    mut block_idx: usize,
    is_fastq: bool,
) -> (usize, usize) {
    let mut block_data = Vec::new();
    let mut decoder = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);
    while block_idx < blocks.len() {
        let b = &blocks[block_idx];
        block_data.resize(b.uncomp_len, 0);
        let block_bytes = &compressed[b.comp_offset..b.comp_offset + b.comp_len];
        let _ = decoder
            .decompress(block_bytes, &mut block_data, true)
            .unwrap();
        decoder.reset();

        let cand = find_record_start(&block_data, 1, is_fastq);
        if cand < block_data.len() {
            return (block_idx, cand);
        }
        block_idx += 1;
    }
    (blocks.len(), 0)
}

#[cfg(feature = "sirocco")]
pub fn run_bgzf_multi_reader_sirocco<
    const CONFIG: helicase::config::Config,
    State,
    Process,
    Reduce,
>(
    compressed: &[u8],
    num_threads: usize,
    is_fastq: bool,
    process: Process,
    reduce: Reduce,
) -> State
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    use std::sync::OnceLock;

    let (blocks, _total_uncomp) = index_bgzf_blocks(compressed);
    if blocks.is_empty() {
        return State::default();
    }
    let num_threads = num_threads.max(1).min(blocks.len());
    let process = &process;
    let reduce = &reduce;
    let blocks_ref = &blocks;

    let split_points: Vec<OnceLock<(usize, usize)>> =
        (0..num_threads + 1).map(|_| OnceLock::new()).collect();
    split_points[0].set((0, 0)).unwrap();
    split_points[num_threads]
        .set((blocks_ref.len(), 0))
        .unwrap();
    let splits_ref = &split_points;

    std::thread::scope(|s| {
        let handles: Vec<_> = (0..num_threads)
            .map(|i| {
                s.spawn(move || {
                    let (start_block, start_skip) = if i == 0 {
                        (0, 0)
                    } else {
                        let target_block = i * blocks_ref.len() / num_threads;
                        let split = find_record_start_in_bgzf(
                            compressed,
                            blocks_ref,
                            target_block,
                            is_fastq,
                        );
                        let _ = splits_ref[i].set(split);
                        split
                    };

                    let (end_block, end_take) = loop {
                        if let Some(&split) = splits_ref[i + 1].get() {
                            break split;
                        }
                        std::thread::yield_now();
                    };

                    if start_block >= blocks_ref.len()
                        || (start_block == end_block && start_skip >= end_take)
                    {
                        return State::default();
                    }

                    let reader = BgzfBlockRangeReader::new(
                        compressed,
                        blocks_ref,
                        start_block,
                        end_block,
                        start_skip,
                        end_take,
                    );
                    let parser = FastxParser::<CONFIG>::from_reader(reader).unwrap();
                    process(parser)
                })
            })
            .collect();

        let mut total = State::default();
        for h in handles {
            let thread_res = h.join().unwrap();
            total = reduce(total, thread_res);
        }
        total
    })
}

// ============================================================================
// Benchmark Harness
// ============================================================================

#[allow(dead_code)]
struct BenchResult {
    label: String,
    elapsed_ms: f64,
    uncompressed_mb_s: f64,
    compressed_mb_s: f64,
    records: usize,
    dna_len: usize,
}

fn median_time(runs: &mut [f64]) -> f64 {
    runs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    runs[runs.len() / 2]
}

fn print_results_table(
    title: &str,
    uncompressed_size: usize,
    compressed_size: usize,
    results: &[BenchResult],
) {
    println!("\n{}", "══════════════════════════════════════════════════════════════════════════════════════════════════".cyan());
    println!(
        "{} | Uncompressed: {:.2} MB | Compressed: {:.2} MB ({:.1}x ratio)",
        title.bold(),
        uncompressed_size as f64 / 1_000_000.0,
        compressed_size as f64 / 1_000_000.0,
        uncompressed_size as f64 / compressed_size as f64
    );
    println!("{}", "──────────────────────────────────────────────────────────────────────────────────────────────────".cyan());
    println!(
        "{:<38} {:>10} {:>16} {:>16} {:>14}",
        "Backend / Architecture", "Time (ms)", "Uncomp (MB/s)", "Comp (MB/s)", "Seq GB/s"
    );
    println!("{}", "──────────────────────────────────────────────────────────────────────────────────────────────────".cyan());

    for r in results {
        let seq_gb_s = (r.dna_len as f64 / (r.elapsed_ms / 1000.0)) / 1_000_000_000.0;
        let label_colored = if r.label.contains("BGZF stream 8T") {
            r.label.green().bold()
        } else if r.label.contains("BGZF stream") {
            r.label.bright_green().bold()
        } else if r.label.contains("sirocco BGZF 8T") {
            r.label.yellow().bold()
        } else if r.label.contains("sirocco") {
            r.label.yellow()
        } else if r.label.contains("SIMD") {
            r.label.bright_blue().bold()
        } else {
            r.label.white()
        };

        println!(
            "{:<38} {:>10.2} {:>16.2} {:>16.2} {:>14.2}",
            label_colored, r.elapsed_ms, r.uncompressed_mb_s, r.compressed_mb_s, seq_gb_s
        );
    }
    println!("{}", "══════════════════════════════════════════════════════════════════════════════════════════════════".cyan());
}

fn evaluate_dataset(
    name: &str,
    compressed_path: &Path,
    uncompressed_path: Option<&Path>,
    repeat: usize,
    warmup: usize,
) {
    println!(
        "\n{} {}",
        "▶ Evaluating Dataset:".bold().blue(),
        name.bold()
    );
    let compressed_bytes = std::fs::read(compressed_path).expect("Failed to read compressed file");
    let comp_size = compressed_bytes.len();

    // Baseline uncompressed validation
    let (uncomp_size, expected_records, expected_dna_len) =
        uncompressed_path.map_or((0, 0, 0), |uncomp_p| {
            let uncomp_bytes = std::fs::read(uncomp_p).expect("Failed to read uncompressed file");
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_slice(&uncomp_bytes).unwrap();
            let mut recs = 0;
            let mut len = 0;
            while parser.next().is_some() {
                recs += 1;
                len += parser.get_dna_string().len();
            }
            (uncomp_bytes.len(), recs, len)
        });
    #[cfg(not(feature = "sirocco"))]
    let _ = (expected_records, expected_dna_len);

    let mut results = Vec::new();

    // --- 0. Helicase Pure SIMD Baseline (Zero I/O, already in RAM) ---
    if let Some(uncomp_p) = uncompressed_path {
        let uncomp_bytes = std::fs::read(uncomp_p).unwrap();
        let mut runs = Vec::new();
        let mut last_recs = 0;
        let mut last_len = 0;

        for _ in 0..warmup {
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_slice(&uncomp_bytes).unwrap();
            while parser.next().is_some() {
                black_box(parser.get_dna_string());
            }
        }
        for _ in 0..repeat {
            let t0 = Instant::now();
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_slice(&uncomp_bytes).unwrap();
            let mut recs = 0;
            let mut len = 0;
            while parser.next().is_some() {
                recs += 1;
                len += parser.get_dna_string().len();
            }
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            last_recs = recs;
            last_len = len;
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "Helicase Pure SIMD (RAM slice)".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (uncomp_bytes.len() as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: last_recs,
            dna_len: last_len,
        });
    }

    // --- 1. flate2 (zlib-rs baseline) ---
    {
        // 1a. Pure decompression
        let mut runs = Vec::new();
        let mut out_len = 0;
        for _ in 0..warmup {
            let mut d = flate2::read::MultiGzDecoder::new(&compressed_bytes[..]);
            let mut out = Vec::new();
            let _ = d.read_to_end(&mut out);
        }
        for _ in 0..repeat {
            let t0 = Instant::now();
            let mut d = flate2::read::MultiGzDecoder::new(&compressed_bytes[..]);
            let mut out = Vec::new();
            d.read_to_end(&mut out).unwrap();
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            out_len = out.len();
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "flate2 (zlib-rs) [decompress only]".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (out_len as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: 0,
            dna_len: 0,
        });

        // 1b. Helicase + flate2 streaming
        runs.clear();
        let mut last_recs = 0;
        let mut last_len = 0;
        for _ in 0..repeat {
            let t0 = Instant::now();
            let f = File::open(compressed_path).unwrap();
            let d = flate2::read::MultiGzDecoder::new(BufReader::new(f));
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_reader(d).unwrap();
            let mut recs = 0;
            let mut len = 0;
            while parser.next().is_some() {
                recs += 1;
                len += parser.get_dna_string().len();
            }
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            last_recs = recs;
            last_len = len;
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "Helicase + flate2 (stream)".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (uncomp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: last_recs,
            dna_len: last_len,
        });
    }

    // --- 2. Sirocco ---
    #[cfg(feature = "sirocco")]
    {
        // 2a. Pure decompression (one-shot buffer)
        let mut runs = Vec::new();
        let mut out_len = 0;
        for _ in 0..warmup {
            let _ = sirocco::decompress(sirocco::Format::Gzip, &compressed_bytes, usize::MAX);
        }
        for _ in 0..repeat {
            let t0 = Instant::now();
            let out =
                sirocco::decompress(sirocco::Format::Gzip, &compressed_bytes, usize::MAX).unwrap();
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            out_len = out.len();
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "sirocco 1T [decompress only]".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (out_len as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: 0,
            dna_len: 0,
        });

        // 2b. Helicase + sirocco streaming reader
        runs.clear();
        let mut last_recs = 0;
        let mut last_len = 0;
        for _ in 0..repeat {
            let t0 = Instant::now();
            let f = File::open(compressed_path).unwrap();
            let reader = SiroccoGzipReader::new(BufReader::new(f));
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_reader(reader).unwrap();
            let mut recs = 0;
            let mut len = 0;
            while parser.next().is_some() {
                recs += 1;
                len += parser.get_dna_string().len();
            }
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            last_recs = recs;
            last_len = len;
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "Helicase + sirocco (stream)".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (uncomp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: last_recs,
            dna_len: last_len,
        });

        // 2c. Helicase + sirocco (decompress to RAM + SIMD slice parse)
        runs.clear();
        last_recs = 0;
        last_len = 0;
        for _ in 0..repeat {
            let t0 = Instant::now();
            let decompressed =
                sirocco::decompress(sirocco::Format::Gzip, &compressed_bytes, usize::MAX).unwrap();
            let mut parser = FastxParser::<{ DEFAULT_CONFIG }>::from_slice(&decompressed).unwrap();
            let mut recs = 0;
            let mut len = 0;
            while parser.next().is_some() {
                recs += 1;
                len += parser.get_dna_string().len();
            }
            runs.push(t0.elapsed().as_secs_f64() * 1000.0);
            last_recs = recs;
            last_len = len;
        }
        let med = median_time(&mut runs);
        results.push(BenchResult {
            label: "Helicase + sirocco (decode->RAM)".to_string(),
            elapsed_ms: med,
            uncompressed_mb_s: (uncomp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
            records: last_recs,
            dna_len: last_len,
        });

        // 2d. If BGZF, evaluate block-reused and parallel Rayon Sirocco!
        if is_bgzf(&compressed_bytes) {
            // 1T block reuse
            runs.clear();
            let mut out_len = 0;
            for _ in 0..warmup {
                let _ = decompress_bgzf_sirocco_1t(&compressed_bytes);
            }
            for _ in 0..repeat {
                let t0 = Instant::now();
                let out = decompress_bgzf_sirocco_1t(&compressed_bytes);
                runs.push(t0.elapsed().as_secs_f64() * 1000.0);
                out_len = out.len();
            }
            let med = median_time(&mut runs);
            results.push(BenchResult {
                label: "sirocco BGZF 1T (block reuse)".to_string(),
                elapsed_ms: med,
                uncompressed_mb_s: (out_len as f64 / (med / 1000.0)) / 1_000_000.0,
                compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                records: 0,
                dna_len: 0,
            });

            // Parallel threads (2T, 4T, 8T)
            for threads in [2, 4, 8] {
                runs.clear();
                for _ in 0..repeat {
                    let t0 = Instant::now();
                    let out = decompress_bgzf_sirocco_parallel(&compressed_bytes, threads);
                    runs.push(t0.elapsed().as_secs_f64() * 1000.0);
                    out_len = out.len();
                }
                let med = median_time(&mut runs);
                results.push(BenchResult {
                    label: format!("sirocco BGZF {}T (parallel)", threads),
                    elapsed_ms: med,
                    uncompressed_mb_s: (out_len as f64 / (med / 1000.0)) / 1_000_000.0,
                    compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                    records: 0,
                    dna_len: 0,
                });
            }

            // End-to-end: parallel 8T decode -> RAM + Helicase SIMD parse
            runs.clear();
            last_recs = 0;
            last_len = 0;
            for _ in 0..repeat {
                let t0 = Instant::now();
                let decompressed = decompress_bgzf_sirocco_parallel(&compressed_bytes, 8);
                let mut parser =
                    FastxParser::<{ DEFAULT_CONFIG }>::from_slice(&decompressed).unwrap();
                let mut recs = 0;
                let mut len = 0;
                while parser.next().is_some() {
                    recs += 1;
                    len += parser.get_dna_string().len();
                }
                runs.push(t0.elapsed().as_secs_f64() * 1000.0);
                last_recs = recs;
                last_len = len;
            }
            let med = median_time(&mut runs);
            results.push(BenchResult {
                label: "Helicase + sirocco BGZF 8T (decode->RAM)".to_string(),
                elapsed_ms: med,
                uncompressed_mb_s: (uncomp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                records: last_recs,
                dna_len: last_len,
            });

            // 2e. Streaming Multi-Reader (O(1) RAM: ~1.5 MB total across all threads!)
            let is_fastq_file = {
                let mut block0 = vec![0u8; 64];
                let mut dec = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);
                dec.decompress(&compressed_bytes, &mut block0, false)
                    .map(|p| p.output.first().copied() == Some(b'@'))
                    .unwrap_or(false)
            };

            for threads in [1, 2, 4, 8] {
                runs.clear();
                last_recs = 0;
                last_len = 0;
                for _ in 0..repeat {
                    let t0 = Instant::now();
                    let (recs, len) =
                        run_bgzf_multi_reader_sirocco::<{ DEFAULT_CONFIG }, (usize, usize), _, _>(
                            &compressed_bytes,
                            threads,
                            is_fastq_file,
                            |mut parser| {
                                let mut r = 0;
                                let mut l = 0;
                                while parser.next().is_some() {
                                    r += 1;
                                    l += parser.get_dna_string().len();
                                }
                                (r, l)
                            },
                            |a, b| (a.0 + b.0, a.1 + b.1),
                        );
                    runs.push(t0.elapsed().as_secs_f64() * 1000.0);
                    last_recs = recs;
                    last_len = len;
                    if expected_records > 0 {
                        assert_eq!(
                            recs, expected_records,
                            "Record count mismatch in streaming reader for {} threads!",
                            threads
                        );
                        assert_eq!(
                            len, expected_dna_len,
                            "DNA length mismatch in streaming reader for {} threads!",
                            threads
                        );
                    }
                }
                let med = median_time(&mut runs);
                results.push(BenchResult {
                    label: format!("Helicase BGZF stream {}T [O(1) RAM]", threads),
                    elapsed_ms: med,
                    uncompressed_mb_s: (uncomp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                    compressed_mb_s: (comp_size as f64 / (med / 1000.0)) / 1_000_000.0,
                    records: last_recs,
                    dna_len: last_len,
                });
            }
        }
    }

    print_results_table(name, uncomp_size, comp_size, &results);
}

fn main() {
    let repeat = 3;
    let warmup = 1;

    println!("{}", "==================================================================================================".bold());
    println!(
        "{}",
        "   Helicase Decompression Evaluation: Sirocco vs flate2"
            .bold()
            .green()
    );
    println!("{}", "==================================================================================================".bold());

    // 1. Short Reads — Standard Gzip (Illumina FASTQ)
    evaluate_dataset(
        "Short Reads (Standard Gzip, 150bp Illumina)",
        Path::new("/tmp/test_real_standard.fastq.gz"),
        Some(Path::new("/tmp/test_real.fastq")),
        repeat,
        warmup,
    );

    // 2. Short Reads — BGZF (Illumina FASTQ)
    evaluate_dataset(
        "Short Reads (BGZF Gzip, 150bp Illumina)",
        Path::new("/tmp/test_real_bgzip.fastq.gz"),
        Some(Path::new("/tmp/test_real.fastq")),
        repeat,
        warmup,
    );

    // 3. Long Reads — Standard Gzip (PacBio FASTA)
    evaluate_dataset(
        "Long Reads (Standard Gzip, ~15kb PacBio)",
        Path::new("sample_long_reads.fa.gz"),
        Some(Path::new("sample_long_reads.fa")),
        repeat,
        warmup,
    );

    // 4. Long Reads — BGZF (PacBio FASTA)
    evaluate_dataset(
        "Long Reads (BGZF Gzip, ~15kb PacBio)",
        Path::new("/tmp/sample_long_reads_133m_bgzip.fa.gz"),
        Some(Path::new("sample_long_reads.fa")),
        repeat,
        warmup,
    );
}
