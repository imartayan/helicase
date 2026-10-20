//! Parallel multi-chunk parsing over contiguous in-memory slices and memory-mapped files.
//!
//! Provides zero-dependency, lock-free work-stealing execution across multiple threads
//! using `std::thread::scope` and `AtomicUsize`.

use super::*;
use crate::config::*;
use crate::input::*;
use core::sync::atomic::{AtomicUsize, Ordering};
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Validates whether `cand` starts a genuine 4-line FASTQ record.
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

#[inline(always)]
fn is_valid_record_start(data: &[u8], cand: usize, format: Format) -> bool {
    match format {
        Format::Fasta => data[cand] == b'>',
        Format::Fastq => data[cand] == b'@' && is_valid_fastq_record_start(data, cand),
    }
}

/// Finds the exact byte index of the first record starting at or after `approx_pos`.
fn find_record_start(data: &[u8], approx_pos: usize, format: Format) -> usize {
    if approx_pos == 0 || approx_pos >= data.len() {
        return approx_pos.min(data.len());
    }
    let mut pos = approx_pos.saturating_sub(1);
    while let Some(newline) = memchr::memchr(b'\n', &data[pos..]) {
        let cand = pos + newline + 1;
        if cand >= data.len() {
            return data.len();
        }
        if is_valid_record_start(data, cand, format) {
            return cand;
        }
        pos = cand;
    }
    data.len()
}

/// Partitions a contiguous buffer into record-aligned `(start, end)` byte slices.
fn partition_chunks(data: &[u8], num_chunks: usize, format: Format) -> Vec<(usize, usize)> {
    if data.is_empty() {
        return Vec::new();
    }
    if num_chunks <= 1 || data.len() < 32 * 1024 {
        return vec![(0, data.len())];
    }
    let approx_chunk_size = data.len() / num_chunks;
    let mut chunks = Vec::with_capacity(num_chunks);
    let mut start = 0;

    for chunk_idx in 1..num_chunks {
        let approx_pos = chunk_idx * approx_chunk_size;
        let end = if approx_pos <= start {
            start
        } else {
            find_record_start(data, approx_pos, format)
        };
        if start < end {
            chunks.push((start, end));
        }
        start = end;
    }
    if start < data.len() {
        chunks.push((start, data.len()));
    }
    chunks
}

/// Executes parallel multi-chunk processing over `data` across `num_threads`.
///
/// Uses an atomic work-stealing counter over fine-grained virtual chunks ($4\times$ oversubscribed)
/// to achieve dynamic load balancing without cross-thread lock contention.
fn par_for_each_chunk<const CONFIG: Config, State, Process, Reduce>(
    data: &[u8],
    num_threads: usize,
    process: Process,
    reduce: Reduce,
) -> io::Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    if data.is_empty() {
        return Ok(State::default());
    }

    let format = Format::try_from_byte(data.first().copied().unwrap_or(0))?;

    let threads = num_threads.max(1);
    if threads == 1 || data.len() < 32 * 1024 {
        let parser = FastxParser::from_slice(data)?;
        return Ok(process(parser));
    }

    let num_chunks = (threads * 4).min((data.len() / (16 * 1024)).max(1));
    let chunks = partition_chunks(data, num_chunks, format);
    if chunks.is_empty() {
        return Ok(State::default());
    }

    if chunks.len() == 1 {
        let parser = FastxParser::from_slice(data)?;
        return Ok(process(parser));
    }

    let worker_threads = threads.min(chunks.len());
    let next_chunk = AtomicUsize::new(0);
    let process = &process;
    let reduce = &reduce;

    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..worker_threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut thread_state: Option<State> = None;
                    loop {
                        let idx = next_chunk.fetch_add(1, Ordering::Relaxed);
                        if idx >= chunks.len() {
                            break;
                        }
                        let (start, end) = chunks[idx];
                        if let Ok(parser) = FastxParser::from_slice(&data[start..end]) {
                            let chunk_res = process(parser);
                            thread_state = Some(match thread_state {
                                Some(acc) => reduce(acc, chunk_res),
                                None => chunk_res,
                            });
                        }
                    }
                    thread_state
                })
            })
            .collect();

        handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .collect::<Vec<_>>()
    });

    let mut total: Option<State> = None;
    for result in results.into_iter().flatten() {
        total = Some(match total {
            Some(acc) => reduce(acc, result),
            None => result,
        });
    }
    Ok(total.unwrap_or_default())
}

/// Probes a file to determine the record boundary starting at or after `approx_pos`.
fn probe_file_record_start(
    file: &mut File,
    approx_pos: u64,
    file_len: u64,
    format: Format,
    buf: &mut [u8],
) -> io::Result<u64> {
    if approx_pos == 0 || approx_pos >= file_len {
        return Ok(approx_pos.min(file_len));
    }

    let mut probe_pos = approx_pos.saturating_sub(1);

    loop {
        file.seek(SeekFrom::Start(probe_pos))?;
        let to_read = ((file_len - probe_pos) as usize).min(buf.len());
        if to_read == 0 {
            return Ok(file_len);
        }
        file.read_exact(&mut buf[..to_read])?;
        let slice = &buf[..to_read];

        let approx_in_slice = if approx_pos > probe_pos {
            (approx_pos - probe_pos) as usize
        } else {
            1
        };

        let cand = find_record_start(slice, approx_in_slice, format);

        if cand < to_read {
            return Ok(probe_pos + cand as u64);
        }
        if probe_pos + to_read as u64 >= file_len {
            return Ok(file_len);
        }
        probe_pos += to_read as u64 - 1;
    }
}

#[inline]
fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1)
}

enum ParallelSource<'a> {
    Slice(&'a [u8]),
    File(PathBuf),
    Mmap(PathBuf),
    Reader(Box<dyn Read + Send + 'a>),
}

/// A parallel multi-threaded parser for FASTA and FASTQ data.
/// Configured via [`ParserOptions`].
///
/// Partitions files or in-memory buffers into record-aligned chunks and distributes them
/// to worker threads. The caller drives the parser loop inside [`reduce`](ParallelFastxParser::reduce)
/// or [`for_each`](ParallelFastxParser::for_each).
///
/// # Example
///
/// ```rust,no_run
/// use helicase::config::*;
/// use helicase::*;
///
/// // set the options of the parser (at compile-time)
/// const CONFIG: Config = ParserOptions::default().config();
///
/// fn main() -> std::io::Result<()> {
///     let path = "...";
///
///     // process records in parallel with manual loop and reduce
///     let total_len = ParallelFastxParser::<CONFIG>::from_file(path)?
///         .threads(8)
///         .reduce(
///             |mut parser| {
///                 let mut len = 0;
///                 while let Some(_event) = parser.next() {
///                     len += parser.get_dna_string().len();
///                 }
///                 len
///             },
///             |acc, val| acc + val,
///         )?;
///     # Ok(())
/// }
/// ```
pub struct ParallelFastxParser<'a, const CONFIG: Config> {
    source: ParallelSource<'a>,
    threads: usize,
    chunks: Option<usize>,
}

impl<'a, const CONFIG: Config> ParallelFastxParser<'a, CONFIG> {
    /// Creates a parallel parser over an in-memory byte slice.
    #[inline]
    pub fn from_slice(slice: &'a [u8]) -> Self {
        Self {
            source: ParallelSource::Slice(slice),
            threads: default_threads(),
            chunks: None,
        }
    }

    /// Creates a parallel parser for a file on disk using the streaming multi-reader engine.
    ///
    /// Reads records in parallel using positional reads with a 128 KiB buffer per thread,
    /// bounding memory to ~1 MB total RAM across 8 threads.
    ///
    /// Also see [`from_file_mmap`](Self::from_file_mmap) for zero-copy memory mapping on local SSD/NVMe.
    pub fn from_file(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let _ = File::open(&path)?;
        Ok(Self {
            source: ParallelSource::File(path),
            threads: default_threads(),
            chunks: None,
        })
    }

    /// Creates a parallel parser for a file on disk using zero-copy memory mapping (`mmap`).
    ///
    /// Yields maximum throughput on local SSD/NVMe storage.
    ///
    /// Also see [`from_file`](Self::from_file) for streaming multi-reader mode with bounded memory.
    pub fn from_file_mmap(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let _ = File::open(&path)?;
        Ok(Self {
            source: ParallelSource::Mmap(path),
            threads: default_threads(),
            chunks: None,
        })
    }

    /// Creates a parallel parser for an arbitrary streaming reader.
    ///
    /// Non-seekable streaming readers currently execute in single-threaded streaming mode.
    #[inline]
    pub fn from_reader<Reader: Read + Send + 'a>(reader: Reader) -> Self {
        Self {
            source: ParallelSource::Reader(Box::new(reader)),
            threads: default_threads(),
            chunks: None,
        }
    }

    /// Sets the number of worker threads.
    #[inline]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads.max(1);
        self
    }

    /// Sets a custom number of chunks to partition the workload into.
    ///
    /// Defaults to `threads * 4` for dynamic work stealing.
    #[inline]
    pub fn chunks(mut self, num_chunks: usize) -> Self {
        self.chunks = Some(num_chunks.max(1));
        self
    }

    /// Executes parallel processing across worker threads, allowing the caller to drive
    /// the parser loop directly on each chunk's [`FastxParser`], and combining thread results with `reduce`.
    ///
    /// The reduction `State` type is inferred automatically from the return type of `process`,
    /// and must implement [`Send`] and [`Default`]. Because empty inputs and initial thread
    /// accumulators use [`Default::default()`], no initial state closure argument is required.
    ///
    /// # Examples
    ///
    /// **Summing sequence lengths:**
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    ///
    /// const CONFIG: Config = ParserOptions::default().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let total_len = ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .reduce(
    ///         |mut parser| {
    ///             let mut len = 0;
    ///             while let Some(_event) = parser.next() {
    ///                 len += parser.get_dna_len();
    ///             }
    ///             len
    ///         },
    ///         |acc, val| acc + val,
    ///     )?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// **Collecting per-record lengths in parallel:**
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    ///
    /// const CONFIG: Config = ParserOptions::default().ignore_headers().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let lengths: Vec<usize> = ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .reduce(
    ///         |mut parser| {
    ///             let mut local = Vec::new();
    ///             while let Some(_event) = parser.next() {
    ///                 local.push(parser.get_dna_len());
    ///             }
    ///             local
    ///         },
    ///         |mut acc, chunk| {
    ///             acc.extend(chunk);
    ///             acc
    ///         },
    ///     )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn reduce<State, Process, Reduce>(
        self,
        process: Process,
        reduce: Reduce,
    ) -> io::Result<State>
    where
        State: Send + Default,
        Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
        Reduce: Fn(State, State) -> State + Sync,
    {
        match self.source {
            ParallelSource::Slice(slice) => {
                par_for_each_chunk::<CONFIG, State, _, _>(slice, self.threads, process, reduce)
            }
            ParallelSource::Mmap(path) => {
                let file = File::open(&path)?;
                let mmap = unsafe { memmap2::Mmap::map(&file)? };
                par_for_each_chunk::<CONFIG, State, _, _>(&mmap, self.threads, process, reduce)
            }
            ParallelSource::File(path) => run_multi_reader::<CONFIG, State, _, _>(
                &path,
                self.threads,
                self.chunks,
                process,
                reduce,
            ),
            ParallelSource::Reader(reader) => {
                if self.threads > 1 {
                    eprintln!(
                        "Warning: ParallelFastxParser currently processes streaming readers in single-threaded mode."
                    );
                }
                let parser = FastxParser::<CONFIG>::from_reader(reader)?;
                Ok(process(parser))
            }
        }
    }

    /// Executes parallel processing across chunks without reduction (a reduce on unit `()`).
    ///
    /// The caller receives each chunk's [`FastxParser`] and drives the parser loop directly.
    pub fn for_each<Process>(self, process: Process) -> io::Result<()>
    where
        Process: Fn(FastxParser<'_, CONFIG>) + Sync,
    {
        self.reduce(process, |(), ()| ())
    }

    /// Executes parallel processing across worker threads using a Paraseq [`ParallelProcessor`](paraseq::parallel::ParallelProcessor).
    ///
    /// Drives worker threads via Helicase's lock-free work-stealing engine, passing each
    /// record as a zero-copy [`paraseq::Record`] reference into [`process_record`](paraseq::parallel::ParallelProcessor::process_record).
    #[cfg(feature = "paraseq")]
    pub fn process_parallel<P>(self, processor: &mut P) -> io::Result<()>
    where
        P: for<'b, 'c> paraseq::parallel::ParallelProcessor<&'b FastxParser<'c, CONFIG>> + Sync,
    {
        let proc = processor.clone();
        self.for_each(move |mut parser| {
            let mut local_proc = proc.clone();
            while let Some(_event) = parser.next() {
                if let Err(e) = local_proc.process_record(&parser) {
                    eprintln!("Paraseq processor error: {e}");
                    break;
                }
            }
            let _ = local_proc.on_thread_complete();
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FileKind {
    Uncompressed,
    Gzip,
    Zstd,
    OtherCompressed,
}

fn detect_file_kind(path: &Path) -> io::Result<FileKind> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 4];
    let n = file.read(&mut magic)?;
    if n == 0 {
        return Ok(FileKind::Uncompressed);
    }
    if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        return Ok(FileKind::Gzip);
    }
    // Zstandard: 0xFD2FB528 (little-endian: 0x28, 0xb5, 0x2f, 0xfd) or skippable frame (0x184D2A50..0x184D2A5F)
    if n >= 4
        && (magic == [0x28, 0xb5, 0x2f, 0xfd]
            || (magic[3] == 0x18
                && magic[2] == 0x4d
                && magic[1] == 0x2a
                && (magic[0] & 0xf0) == 0x50))
    {
        return Ok(FileKind::Zstd);
    }
    // Check file extension as fallback for Zstandard
    if let Some(ext) = path.extension().and_then(|e| e.to_str())
        && (ext.eq_ignore_ascii_case("zst") || ext.eq_ignore_ascii_case("zstd"))
    {
        return Ok(FileKind::Zstd);
    }
    // Bzip2: 'B', 'Z'
    if n >= 2 && magic[0] == b'B' && magic[1] == b'Z' {
        return Ok(FileKind::OtherCompressed);
    }
    // Xz: 0xFD, '7', 'z', 'X', 'Z', 0x00
    if n >= 3 && magic[0] == 0xfd && magic[1] == b'7' && magic[2] == b'z' {
        return Ok(FileKind::OtherCompressed);
    }
    Ok(FileKind::Uncompressed)
}

fn run_multi_reader<const CONFIG: Config, State, Process, Reduce>(
    path: &Path,
    num_threads: usize,
    custom_chunks: Option<usize>,
    process: Process,
    reduce: Reduce,
) -> io::Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let num_threads = num_threads.max(1);
    match detect_file_kind(path)? {
        FileKind::Uncompressed => {
            let mut file = File::open(path)?;
            let file_len = file.metadata()?.len();
            if file_len == 0 {
                return Ok(State::default());
            }

            let mut first_byte = [0u8; 1];
            file.read_exact(&mut first_byte)?;
            let format = Format::try_from_byte(first_byte[0])?;

            if num_threads == 1 || file_len < 32 * 1024 {
                file.seek(SeekFrom::Start(0))?;
                let parser = FastxParser::<CONFIG>::from_reader(file)?;
                return Ok(process(parser));
            }

            let num_chunks = custom_chunks.unwrap_or(num_threads * 4).max(1);
            let approx_chunk_size = file_len / num_chunks as u64;
            let worker_threads = num_threads.min(num_chunks);
            let next_chunk = AtomicUsize::new(0);
            let process = &process;
            let reduce = &reduce;

            let splits: Vec<std::sync::OnceLock<u64>> = (0..=num_chunks)
                .map(|_| std::sync::OnceLock::new())
                .collect();
            splits[0].set(0).unwrap();
            splits[num_chunks].set(file_len).unwrap();
            let splits_ref = &splits;

            let results = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..worker_threads)
                    .map(|_| {
                        scope.spawn(|| -> io::Result<Option<State>> {
                            let mut thread_state: Option<State> = None;
                            let mut file = File::open(path)?;
                            let mut probe_buf = vec![0u8; 64 * 1024];

                            loop {
                                let idx = next_chunk.fetch_add(1, Ordering::Relaxed);
                                if idx >= num_chunks {
                                    break;
                                }
                                let get_split =
                                    |chunk_idx: usize, file: &mut File, buf: &mut [u8]| -> u64 {
                                        *splits_ref[chunk_idx].get_or_init(|| {
                                            let approx_pos = chunk_idx as u64 * approx_chunk_size;
                                            probe_file_record_start(
                                                file, approx_pos, file_len, format, buf,
                                            )
                                            .unwrap_or(approx_pos)
                                        })
                                    };

                                let start = get_split(idx, &mut file, &mut probe_buf);
                                let end = get_split(idx + 1, &mut file, &mut probe_buf);

                                if start >= end {
                                    continue;
                                }
                                let take_len = end - start;
                                file.seek(SeekFrom::Start(start))?;
                                let reader = (&mut file).take(take_len);
                                let parser = FastxParser::<CONFIG>::from_reader(reader)?;
                                let chunk_res = process(parser);
                                thread_state = Some(match thread_state {
                                    Some(acc) => reduce(acc, chunk_res),
                                    None => chunk_res,
                                });
                            }
                            Ok(thread_state)
                        })
                    })
                    .collect();

                handles
                    .into_iter()
                    .map(|handle| handle.join().unwrap())
                    .collect::<Vec<_>>()
            });

            let mut total: Option<State> = None;
            for result in results {
                if let Some(res) = result? {
                    total = Some(match total {
                        Some(acc) => reduce(acc, res),
                        None => res,
                    });
                }
            }
            Ok(total.unwrap_or_default())
        }
        FileKind::Gzip => {
            #[cfg(feature = "sirocco")]
            {
                let file = File::open(path)?;
                let mmap = unsafe { memmap2::Mmap::map(&file)? };
                if is_bgzf_slice(&mmap) {
                    return run_bgzf_multi_reader_sirocco::<CONFIG, State, _, _>(
                        &mmap,
                        num_threads,
                        process,
                        reduce,
                    );
                }
                let reader = sirocco::inflate::Reader::new(
                    sirocco::Format::Gzip,
                    std::io::BufReader::new(file),
                );
                let parser = FastxParser::<CONFIG>::from_reader(reader)?;
                Ok(process(parser))
            }
            #[cfg(not(feature = "sirocco"))]
            {
                let file = File::open(path)?;
                let parser = FastxParser::<CONFIG>::from_reader(file)?;
                Ok(process(parser))
            }
        }
        FileKind::Zstd => {
            #[cfg(feature = "zstd")]
            {
                let file = File::open(path)?;
                let mmap = unsafe { memmap2::Mmap::map(&file)? };
                let frames = index_zstd_frames(&mmap);
                if frames.len() > 1 && num_threads > 1 {
                    return run_zstd_multi_reader::<CONFIG, State, _, _>(
                        &mmap,
                        &frames,
                        num_threads,
                        process,
                        reduce,
                    );
                }
                let parser = FastxParser::<CONFIG>::from_file(path)?;
                Ok(process(parser))
            }
            #[cfg(not(feature = "zstd"))]
            {
                let file = File::open(path)?;
                let parser = FastxParser::<CONFIG>::from_reader(file)?;
                Ok(process(parser))
            }
        }
        FileKind::OtherCompressed => {
            let file = File::open(path)?;
            let parser = FastxParser::<CONFIG>::from_reader(file)?;
            Ok(process(parser))
        }
    }
}

#[cfg(feature = "sirocco")]
#[inline]
pub(crate) fn is_bgzf_slice(bytes: &[u8]) -> bool {
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
struct BgzfBlock {
    comp_offset: usize,
    comp_len: usize,
    uncomp_len: usize,
}

#[cfg(feature = "sirocco")]
fn index_bgzf_blocks(compressed: &[u8]) -> Vec<BgzfBlock> {
    let mut blocks = Vec::with_capacity(compressed.len() / 2048);
    let mut comp_offset = 0;

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
            blocks.push(BgzfBlock {
                comp_offset,
                comp_len,
                uncomp_len,
            });
        }
        comp_offset += comp_len;
    }
    blocks
}

#[cfg(feature = "sirocco")]
struct BgzfBlockRangeReader<'a> {
    compressed: &'a [u8],
    blocks: &'a [BgzfBlock],
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
    fn new(
        compressed: &'a [u8],
        blocks: &'a [BgzfBlock],
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
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
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
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:?}")))?;
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
fn find_record_start_in_bgzf(
    compressed: &[u8],
    blocks: &[BgzfBlock],
    mut block_idx: usize,
    format: Format,
) -> (usize, usize) {
    let mut block_data = Vec::new();
    let mut decoder = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);
    while block_idx < blocks.len() {
        let b = &blocks[block_idx];
        block_data.resize(b.uncomp_len, 0);
        let block_bytes = &compressed[b.comp_offset..b.comp_offset + b.comp_len];
        let _ = decoder.decompress(block_bytes, &mut block_data, true);
        decoder.reset();

        let cand = find_record_start(&block_data, 1, format);
        if cand < block_data.len() {
            return (block_idx, cand);
        }
        block_idx += 1;
    }
    (blocks.len(), 0)
}

#[cfg(feature = "sirocco")]
fn run_bgzf_multi_reader_sirocco<const CONFIG: Config, State, Process, Reduce>(
    compressed: &[u8],
    num_threads: usize,
    process: Process,
    reduce: Reduce,
) -> io::Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let blocks = index_bgzf_blocks(compressed);
    if blocks.is_empty() {
        return Ok(State::default());
    }

    let format = {
        let mut block0 = vec![0u8; 64];
        let mut dec = sirocco::inflate::Decoder::new(sirocco::Format::Gzip);
        let first_byte = dec
            .decompress(&compressed[..blocks[0].comp_len], &mut block0, false)
            .ok()
            .and_then(|p| p.output.first().copied())
            .unwrap_or(0);
        Format::try_from_byte(first_byte)?
    };

    let num_threads = num_threads.max(1).min(blocks.len());
    let process = &process;
    let reduce = &reduce;
    let blocks_ref = &blocks;

    let split_points: Vec<std::sync::OnceLock<(usize, usize)>> = (0..=num_threads)
        .map(|_| std::sync::OnceLock::new())
        .collect();
    split_points[0].set((0, 0)).unwrap();
    split_points[num_threads]
        .set((blocks_ref.len(), 0))
        .unwrap();
    let splits_ref = &split_points;

    let get_split = |idx: usize| -> (usize, usize) {
        *splits_ref[idx].get_or_init(|| {
            let target_block = idx * blocks_ref.len() / num_threads;
            find_record_start_in_bgzf(compressed, blocks_ref, target_block, format)
        })
    };

    let results = std::thread::scope(|s| {
        let handles: Vec<_> = (0..num_threads)
            .map(|i| {
                s.spawn(move || -> io::Result<State> {
                    let (start_block, start_skip) = get_split(i);
                    let (end_block, end_take) = get_split(i + 1);

                    if start_block >= blocks_ref.len()
                        || (start_block == end_block && start_skip >= end_take)
                    {
                        return Ok(State::default());
                    }
                    let reader = BgzfBlockRangeReader::new(
                        compressed,
                        blocks_ref,
                        start_block,
                        end_block,
                        start_skip,
                        end_take,
                    );
                    let parser = FastxParser::<CONFIG>::from_reader(reader)?;
                    Ok(process(parser))
                })
            })
            .collect();

        handles
            .into_iter()
            .map(|h| h.join().unwrap())
            .collect::<Vec<_>>()
    });

    let mut total = State::default();
    for res in results {
        total = reduce(total, res?);
    }
    Ok(total)
}

#[cfg(feature = "zstd")]
#[derive(Clone, Copy)]
struct ZstdFrame {
    comp_offset: usize,
    comp_len: usize,
    uncomp_len: usize,
}

#[cfg(feature = "zstd")]
fn index_zstd_frames(compressed: &[u8]) -> Vec<ZstdFrame> {
    let mut frames = Vec::with_capacity(compressed.len() / (1024 * 1024) + 1);
    let mut comp_offset = 0;

    while comp_offset + 4 <= compressed.len() {
        let slice = &compressed[comp_offset..];
        // Check for skippable frame: 0x184D2A50 to 0x184D2A5F
        // Format: 4-byte magic + 4-byte frame size (u32 little-endian) + payload
        if slice.len() >= 8
            && slice[3] == 0x18
            && slice[2] == 0x4d
            && slice[1] == 0x2a
            && (slice[0] & 0xf0) == 0x50
        {
            let user_size = u32::from_le_bytes([slice[4], slice[5], slice[6], slice[7]]) as usize;
            comp_offset += 8 + user_size;
            continue;
        }

        match zstd::zstd_safe::find_frame_compressed_size(slice) {
            Ok(comp_len) if comp_len > 0 && comp_offset + comp_len <= compressed.len() => {
                let uncomp_len = zstd::zstd_safe::get_frame_content_size(slice)
                    .ok()
                    .flatten()
                    .unwrap_or(0) as usize;
                frames.push(ZstdFrame {
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

#[cfg(feature = "zstd")]
struct ZstdFrameRangeReader<'a> {
    compressed: &'a [u8],
    frames: &'a [ZstdFrame],
    current_frame: usize,
    start_frame: usize,
    end_frame: usize,
    start_skip: usize,
    end_take: usize,

    decompressor: zstd::bulk::Decompressor<'static>,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
}

#[cfg(feature = "zstd")]
impl<'a> ZstdFrameRangeReader<'a> {
    fn new(
        compressed: &'a [u8],
        frames: &'a [ZstdFrame],
        start_frame: usize,
        end_frame: usize,
        start_skip: usize,
        end_take: usize,
    ) -> io::Result<Self> {
        Ok(Self {
            compressed,
            frames,
            current_frame: start_frame,
            start_frame,
            end_frame,
            start_skip,
            end_take,
            decompressor: zstd::bulk::Decompressor::new()?,
            buf: Vec::new(),
            pos: 0,
            len: 0,
        })
    }
}

#[cfg(feature = "zstd")]
impl<'a> Read for ZstdFrameRangeReader<'a> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.len {
            if self.current_frame > self.end_frame
                || (self.current_frame == self.end_frame && self.end_take == 0)
                || self.current_frame >= self.frames.len()
            {
                return Ok(0);
            }
            let f = &self.frames[self.current_frame];
            let frame_bytes = &self.compressed[f.comp_offset..f.comp_offset + f.comp_len];

            // Ignore skippable frames (magic 0x184D2A50..=0x184D2A5F)
            if frame_bytes.len() >= 4
                && (frame_bytes[0] & 0xF0) == 0x50
                && frame_bytes[1] == 0x2A
                && frame_bytes[2] == 0x4D
                && frame_bytes[3] == 0x18
            {
                self.current_frame += 1;
                continue;
            }

            if f.uncomp_len > 0 {
                self.buf.resize(f.uncomp_len, 0);
                let n = self
                    .decompressor
                    .decompress_to_buffer(frame_bytes, &mut self.buf)?;
                self.buf.truncate(n);
            } else {
                self.buf.clear();
                let mut decoder = zstd::stream::read::Decoder::new(frame_bytes)?;
                decoder.read_to_end(&mut self.buf)?;
            }

            let mut start = 0;
            let mut end = self.buf.len();

            if self.current_frame == self.start_frame {
                start = self.start_skip;
            }
            if self.current_frame == self.end_frame {
                end = self.end_take;
            }

            self.pos = start;
            self.len = end;
            self.current_frame += 1;
        }

        let n = (self.len - self.pos).min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

#[cfg(feature = "zstd")]
fn find_record_start_in_zstd(
    compressed: &[u8],
    frames: &[ZstdFrame],
    mut frame_idx: usize,
    format: Format,
) -> (usize, usize) {
    let mut frame_data = Vec::new();
    let mut decompressor = match zstd::bulk::Decompressor::new() {
        Ok(d) => d,
        Err(_) => return (frames.len(), 0),
    };

    while frame_idx < frames.len() {
        let f = &frames[frame_idx];
        let frame_bytes = &compressed[f.comp_offset..f.comp_offset + f.comp_len];

        let decomp_ok = if f.uncomp_len > 0 {
            frame_data.resize(f.uncomp_len, 0);
            if let Ok(n) = decompressor.decompress_to_buffer(frame_bytes, &mut frame_data) {
                frame_data.truncate(n);
                true
            } else {
                false
            }
        } else {
            frame_data.clear();
            zstd::stream::read::Decoder::new(frame_bytes)
                .and_then(|mut r| r.read_to_end(&mut frame_data))
                .is_ok()
        };

        if decomp_ok {
            let cand = find_record_start(&frame_data, 1, format);
            if cand < frame_data.len() {
                return (frame_idx, cand);
            }
        }
        frame_idx += 1;
    }
    (frames.len(), 0)
}

#[cfg(feature = "zstd")]
fn run_zstd_multi_reader<const CONFIG: Config, State, Process, Reduce>(
    compressed: &[u8],
    frames: &[ZstdFrame],
    num_threads: usize,
    process: Process,
    reduce: Reduce,
) -> io::Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> State + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let format = {
        let f0 = &frames[0];
        let frame_bytes = &compressed[f0.comp_offset..f0.comp_offset + f0.comp_len];
        let mut header = [0u8; 64];
        let n = zstd::stream::read::Decoder::new(frame_bytes)
            .and_then(|mut r| r.read(&mut header))
            .unwrap_or(0);
        let first_byte = if n > 0 { header[0] } else { 0 };
        Format::try_from_byte(first_byte)?
    };

    let num_threads = num_threads.max(1).min(frames.len());
    let process = &process;
    let reduce = &reduce;
    let frames_ref = frames;

    let split_points: Vec<std::sync::OnceLock<(usize, usize)>> = (0..=num_threads)
        .map(|_| std::sync::OnceLock::new())
        .collect();
    split_points[0].set((0, 0)).unwrap();
    split_points[num_threads]
        .set((frames_ref.len(), 0))
        .unwrap();
    let splits_ref = &split_points;

    let get_split = |idx: usize| -> (usize, usize) {
        *splits_ref[idx].get_or_init(|| {
            let target_frame = idx * frames_ref.len() / num_threads;
            find_record_start_in_zstd(compressed, frames_ref, target_frame, format)
        })
    };

    let results = std::thread::scope(|s| {
        let handles: Vec<_> = (0..num_threads)
            .map(|i| {
                s.spawn(move || -> io::Result<State> {
                    let (start_frame, start_skip) = get_split(i);
                    let (end_frame, end_take) = get_split(i + 1);

                    if start_frame >= frames_ref.len()
                        || (start_frame == end_frame && start_skip >= end_take)
                    {
                        return Ok(State::default());
                    }
                    let reader = ZstdFrameRangeReader::new(
                        compressed,
                        frames_ref,
                        start_frame,
                        end_frame,
                        start_skip,
                        end_take,
                    )?;
                    let parser = FastxParser::<CONFIG>::from_reader(reader)?;
                    Ok(process(parser))
                })
            })
            .collect();

        handles
            .into_iter()
            .map(|h| h.join().unwrap())
            .collect::<Vec<_>>()
    });

    let mut total = State::default();
    for res in results {
        total = reduce(total, res?);
    }
    Ok(total)
}
