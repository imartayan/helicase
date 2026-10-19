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

/// Finds the exact byte index of the first FASTA record (`>`) starting at or after `approx_pos`.
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

/// Finds the exact byte index of the first FASTQ record (`@`) starting at or after `approx_pos`.
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

/// Partitions a contiguous buffer into record-aligned `(start, end)` byte slices.
fn partition_chunks(data: &[u8], num_chunks: usize, format: Format) -> Vec<(usize, usize)> {
    if data.is_empty() {
        return Vec::new();
    }
    if num_chunks <= 1 || data.len() < 32 * 1024 {
        return vec![(0, data.len())];
    }
    let mut split_points = Vec::with_capacity(num_chunks + 1);
    split_points.push(0);

    let approx_chunk_size = data.len() / num_chunks;
    let mut last_sync = 0;

    for chunk_idx in 1..num_chunks {
        let approx_pos = chunk_idx * approx_chunk_size;
        if approx_pos <= last_sync {
            split_points.push(last_sync);
            continue;
        }
        let sync_pos = match format {
            Format::Fasta => find_fasta_record_start(data, approx_pos),
            Format::Fastq => find_fastq_record_start(data, approx_pos),
        };
        last_sync = sync_pos;
        split_points.push(sync_pos);
    }
    split_points.push(data.len());

    let mut chunks = Vec::with_capacity(num_chunks);
    for window in split_points.windows(2) {
        let (start, end) = (window[0], window[1]);
        if start < end {
            chunks.push((start, end));
        }
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

    let format = match data.first() {
        Some(b'>') => Format::Fasta,
        Some(b'@') => Format::Fastq,
        _ => {
            return Err(io::Error::other(
                "Invalid record start, expected '>' or '@'",
            ));
        }
    };

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
) -> io::Result<u64> {
    if approx_pos == 0 || approx_pos >= file_len {
        return Ok(approx_pos.min(file_len));
    }

    let mut probe_pos = approx_pos.saturating_sub(1);
    let mut buf = vec![0u8; 256 * 1024];

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

        let cand = match format {
            Format::Fasta => find_fasta_record_start(slice, approx_in_slice),
            Format::Fastq => find_fastq_record_start(slice, approx_in_slice),
        };

        if cand < to_read {
            return Ok(probe_pos + cand as u64);
        }
        if probe_pos + to_read as u64 >= file_len {
            return Ok(file_len);
        }
        probe_pos += to_read as u64 - 1;
    }
}

/// Partitions a file on disk into record-aligned `(start, end)` byte ranges via microsecond probe reads.
fn partition_file_chunks(
    path: impl AsRef<Path>,
    num_chunks: usize,
) -> io::Result<(Format, Vec<(u64, u64)>)> {
    let mut file = File::open(&path)?;
    let file_len = file.metadata()?.len();
    if file_len == 0 {
        return Ok((Format::Fasta, Vec::new()));
    }

    let mut first_byte = [0u8; 1];
    file.read_exact(&mut first_byte)?;
    let format = match first_byte[0] {
        b'>' => Format::Fasta,
        b'@' => Format::Fastq,
        _ => {
            return Err(io::Error::other(
                "Invalid record start, expected '>' or '@'",
            ));
        }
    };

    if num_chunks <= 1 || file_len < 32 * 1024 {
        return Ok((format, vec![(0, file_len)]));
    }

    let mut split_points = Vec::with_capacity(num_chunks + 1);
    split_points.push(0);

    let approx_chunk_size = file_len / num_chunks as u64;
    let mut last_sync = 0;

    for chunk_idx in 1..num_chunks {
        let approx_pos = chunk_idx as u64 * approx_chunk_size;
        if approx_pos <= last_sync {
            split_points.push(last_sync);
            continue;
        }
        let sync_pos = probe_file_record_start(&mut file, approx_pos, file_len, format)?;
        last_sync = sync_pos;
        split_points.push(sync_pos);
    }
    split_points.push(file_len);

    let mut chunks = Vec::with_capacity(num_chunks);
    for window in split_points.windows(2) {
        let (start, end) = (window[0], window[1]);
        if start < end {
            chunks.push((start, end));
        }
    }
    Ok((format, chunks))
}

#[inline]
fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1)
}

enum ParallelSource<'a> {
    Slice(&'a [u8]),
    File { path: PathBuf, use_mmap: bool },
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

    /// Creates a parallel parser for a file on disk.
    ///
    /// By default, uses the multi-reader engine (positional reads with 128 KiB buffer per thread,
    /// bounded to ~1 MB total RAM across 8 threads). To enable zero-copy memory mapping instead,
    /// call [`.mmap()`](Self::mmap).
    pub fn from_file(path: impl AsRef<Path>) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let _ = File::open(&path)?;
        Ok(Self {
            source: ParallelSource::File {
                path,
                use_mmap: false,
            },
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

    /// Enables zero-copy memory mapping (`mmap`) for maximum throughput on local SSD/NVMe.
    #[inline]
    pub const fn mmap(mut self) -> Self {
        if let ParallelSource::File {
            ref mut use_mmap, ..
        } = self.source
        {
            *use_mmap = true;
        }
        self
    }

    /// Enables multi-reader mode (default for files: no mmap, 128 KiB buffer per thread).
    #[inline]
    pub const fn reader(mut self) -> Self {
        if let ParallelSource::File {
            ref mut use_mmap, ..
        } = self.source
        {
            *use_mmap = false;
        }
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
    ///             while parser.next().is_some() {
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
    ///             while parser.next().is_some() {
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
            ParallelSource::File { path, use_mmap } => {
                if use_mmap {
                    let file = File::open(&path)?;
                    let mmap = unsafe { memmap2::Mmap::map(&file)? };
                    par_for_each_chunk::<CONFIG, State, _, _>(&mmap, self.threads, process, reduce)
                } else {
                    run_multi_reader::<CONFIG, State, _, _>(
                        &path,
                        self.threads,
                        self.chunks,
                        process,
                        reduce,
                    )
                }
            }
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
    let num_chunks = custom_chunks.unwrap_or(num_threads * 4).max(1);
    let (_format, chunks) = partition_file_chunks(path, num_chunks)?;

    if chunks.is_empty() {
        return Ok(State::default());
    }

    if num_threads == 1 || chunks.len() == 1 {
        let file = File::open(path)?;
        let parser = FastxParser::<CONFIG>::from_reader(file)?;
        return Ok(process(parser));
    }

    let worker_threads = num_threads.min(chunks.len());
    let next_chunk = AtomicUsize::new(0);
    let process = &process;
    let reduce = &reduce;

    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..worker_threads)
            .map(|_| {
                scope.spawn(|| -> io::Result<Option<State>> {
                    let mut thread_state: Option<State> = None;
                    let mut file = File::open(path)?;

                    loop {
                        let idx = next_chunk.fetch_add(1, Ordering::Relaxed);
                        if idx >= chunks.len() {
                            break;
                        }
                        let (start, end) = chunks[idx];
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
