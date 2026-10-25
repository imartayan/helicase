//! Parallel multi-chunk parsing.
//!
//! Provides lock-free work-stealing execution across multiple threads,
//! for both reader and memory-mapped files.

use crate::config::*;
use crate::input::*;
use crate::parser::*;
use core::cell::RefCell;
use core::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

#[cfg(any(feature = "sirocco", feature = "zstd"))]
mod frames;
#[cfg(feature = "zstd")]
pub(crate) use frames::index_zstd_frames;

/// Processes `num_chunks` chunks across `threads` with work stealing, and reduces the results.
///
/// Chunk `i` goes from split `i` to split `i + 1`. The first and last splits are given, the others
/// are computed once by `split`, when a thread first needs them. Empty chunks are skipped.
/// Thread `t` starts with chunk `t`, so results are in input order with one chunk per thread.
/// If a chunk fails, the other threads stop and the first error is returned.
fn run_chunks<P, State, Split, Parse, Reduce>(
    num_chunks: usize,
    threads: usize,
    (first, last): (P, P),
    split: Split,
    parse_chunk: Parse,
    reduce: Reduce,
) -> Result<State>
where
    P: Copy + PartialOrd + Send + Sync,
    State: Send + Default,
    Split: Fn(usize) -> P + Sync,
    Parse: Fn(P, P) -> Result<State> + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let splits: Vec<OnceLock<P>> = (0..=num_chunks).map(|_| OnceLock::new()).collect();
    let _ = splits[0].set(first);
    let _ = splits[num_chunks].set(last);
    let split_at = &|i: usize| *splits[i].get_or_init(|| split(i));
    let workers = threads.min(num_chunks);
    let next_chunk = &AtomicUsize::new(workers);
    let stop = &AtomicBool::new(false);
    let (parse_chunk, reduce) = (&parse_chunk, &reduce);
    let results = std::thread::scope(|scope| {
        // Collected to spawn all the threads before joining them.
        #[allow(clippy::needless_collect)]
        let handles: Vec<_> = (0..workers)
            .map(|first_chunk| {
                scope.spawn(move || -> Result<Option<State>> {
                    let mut thread_state = None;
                    let mut idx = first_chunk;
                    while idx < num_chunks && !stop.load(Ordering::Relaxed) {
                        let (start, end) = (split_at(idx), split_at(idx + 1));
                        if start < end {
                            let res = parse_chunk(start, end)
                                .inspect_err(|_| stop.store(true, Ordering::Relaxed))?;
                            thread_state = Some(match thread_state {
                                Some(acc) => reduce(acc, res),
                                None => res,
                            });
                        }
                        idx = next_chunk.fetch_add(1, Ordering::Relaxed);
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
    let mut total = None;
    for thread_state in results {
        if let Some(res) = thread_state? {
            total = Some(match total {
                Some(acc) => reduce(acc, res),
                None => res,
            });
        }
    }
    Ok(total.unwrap_or_default())
}

fn run_slice<const CONFIG: Config, State, Process, Reduce>(
    data: &[u8],
    threads: usize,
    process: Process,
    reduce: Reduce,
) -> Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> Result<State> + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    let Some(&first_byte) = data.first() else {
        return Ok(State::default());
    };
    let format = Format::try_from_byte(first_byte)?;
    if threads == 1 || data.len() < 32 * 1024 {
        return process(FastxParser::from_slice(data)?);
    }
    let num_chunks = (threads * 4).min(data.len() / (16 * 1024));
    let chunk_len = data.len() / num_chunks;
    run_chunks(
        num_chunks,
        threads,
        (0, data.len()),
        |i| {
            format
                .find_record_start(data, i * chunk_len)
                .unwrap_or(data.len())
        },
        |start, end| process(FastxParser::from_slice(&data[start..end])?),
        reduce,
    )
}

/// Probes a file to determine the record boundary starting at or after `approx_pos`.
fn probe_file_record_start(
    file: &mut File,
    approx_pos: u64,
    file_len: u64,
    format: Format,
) -> Result<u64> {
    if approx_pos == 0 || approx_pos >= file_len {
        return Ok(approx_pos.min(file_len));
    }
    let mut buf = vec![0u8; 64 * 1024];
    let mut probe_pos = approx_pos.saturating_sub(1);
    loop {
        file.seek(SeekFrom::Start(probe_pos))?;
        let to_read = ((file_len - probe_pos) as usize).min(buf.len());
        if to_read == 0 {
            return Ok(file_len);
        }
        file.read_exact(&mut buf[..to_read])?;
        let approx_in_slice = approx_pos.saturating_sub(probe_pos) as usize;
        if let Some(cand) = format.find_record_start(&buf[..to_read], approx_in_slice) {
            return Ok(probe_pos + cand as u64);
        }
        if probe_pos + to_read as u64 >= file_len {
            return Ok(file_len);
        }
        probe_pos += to_read as u64 - 1;
    }
}

fn run_file<const CONFIG: Config, State, Process, Reduce>(
    path: &Path,
    threads: usize,
    custom_chunks: Option<usize>,
    process: Process,
    reduce: Reduce,
) -> Result<State>
where
    State: Send + Default,
    Process: Fn(FastxParser<'_, CONFIG>) -> Result<State> + Sync,
    Reduce: Fn(State, State) -> State + Sync,
{
    match CompressionFormat::from_path(path)? {
        None => {}
        #[cfg(feature = "sirocco")]
        Some(CompressionFormat::Bgzf) if threads > 1 => {
            let mmap = unsafe { memmap2::Mmap::map(&File::open(path)?)? };
            let blocks = frames::index_bgzf_blocks(&mmap);
            return frames::run_frames::<frames::BgzfDecoder, CONFIG, _, _, _>(
                &mmap, &blocks, threads, process, reduce,
            );
        }
        #[cfg(feature = "zstd")]
        Some(CompressionFormat::Pzstd) if threads > 1 => {
            let mmap = unsafe { memmap2::Mmap::map(&File::open(path)?)? };
            let zstd_frames = frames::index_zstd_frames(&mmap);
            if zstd_frames.len() > 1 {
                return frames::run_frames::<frames::ZstdDecoder, CONFIG, _, _, _>(
                    &mmap,
                    &zstd_frames,
                    threads,
                    process,
                    reduce,
                );
            }
            return process(FastxParser::from_file(path)?);
        }
        Some(_) => return process(FastxParser::from_file(path)?),
    }

    let mut file = File::open(path)?;
    let file_len = file.metadata()?.len();
    if file_len == 0 {
        return Ok(State::default());
    }
    let mut first_byte = [0u8; 1];
    file.read_exact(&mut first_byte)?;
    let format = Format::try_from_byte(first_byte[0])?;
    if threads == 1 || file_len < 32 * 1024 {
        return process(FastxParser::from_file(path)?);
    }
    let num_chunks = custom_chunks.unwrap_or(threads * 4);
    let chunk_len = file_len / num_chunks as u64;
    run_chunks(
        num_chunks,
        threads,
        (0, file_len),
        |i| {
            let pos = i as u64 * chunk_len;
            with_thread_file(path, |f| probe_file_record_start(f, pos, file_len, format))
                .unwrap_or(pos)
        },
        |start, end| {
            with_thread_file(path, |file| {
                file.seek(SeekFrom::Start(start))?;
                process(FastxParser::from_reader(file.take(end - start))?)
            })
        },
        reduce,
    )
}

/// Runs `f` on `path` opened once per worker (reopening it per chunk costs 5-10% at 8 threads).
fn with_thread_file<T>(path: &Path, f: impl FnOnce(&mut File) -> Result<T>) -> Result<T> {
    thread_local!(static FILE: RefCell<Option<File>> = const { RefCell::new(None) });
    FILE.with_borrow_mut(|file| match file {
        Some(file) => f(file),
        None => f(file.insert(File::open(path)?)),
    })
}

/// Trait for return types of [`ParallelFastxParser::for_each`] closures.
///
/// Implemented for `()` (infallible) and [`Result<()>`].
pub trait ParallelResult: crate::sealed::Sealed {
    /// Converts the value returned by the closure into a [`Result`].
    fn into_result(self) -> Result<()>;
}

impl crate::sealed::Sealed for () {}
impl crate::sealed::Sealed for Result<()> {}

impl ParallelResult for () {
    #[inline(always)]
    fn into_result(self) -> Result<()> {
        Ok(())
    }
}

impl ParallelResult for Result<()> {
    #[inline(always)]
    fn into_result(self) -> Result<()> {
        self
    }
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
/// to worker threads. The caller drives the parser loop inside [`for_each`](ParallelFastxParser::for_each)
/// or [`reduce`](ParallelFastxParser::reduce).
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
///     // process records in parallel with for_each
///     ParallelFastxParser::<CONFIG>::from_file(path)?
///         .threads(8)
///         .for_each(|mut parser| {
///             while let Some(event) = parser.next() {
///                 event?;
///                 let header = parser.get_header();
///                 let seq = parser.get_dna_string();
///                 // ...
///             }
///             Ok(())
///         })?;
///     # Ok(())
/// }
/// ```
pub struct ParallelFastxParser<'a, const CONFIG: Config> {
    source: ParallelSource<'a>,
    threads: usize,
    chunks: Option<usize>,
}

impl<'a, const CONFIG: Config> ParallelFastxParser<'a, CONFIG> {
    fn with_source(source: ParallelSource<'a>) -> Self {
        Self {
            source,
            threads: std::thread::available_parallelism().map_or(1, |n| n.get()),
            chunks: None,
        }
    }

    /// Creates a parallel parser over an in-memory byte slice.
    #[inline]
    pub fn from_slice(slice: &'a [u8]) -> Self {
        Self::with_source(ParallelSource::Slice(slice))
    }

    /// Creates a parallel parser for a file on disk using the streaming multi-reader engine.
    ///
    /// Reads records in parallel using positional reads with a 256 KiB buffer per thread,
    /// bounding memory to ~2 MB total RAM across 8 threads.
    ///
    /// Also see [`from_file_mmap`](Self::from_file_mmap) for zero-copy memory mapping on local SSD/NVMe.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        File::open(&path)?;
        Ok(Self::with_source(ParallelSource::File(
            path.as_ref().to_path_buf(),
        )))
    }

    /// Creates a parallel parser for a file on disk using zero-copy memory mapping (`mmap`).
    ///
    /// Yields maximum throughput on local SSD/NVMe storage.
    ///
    /// Also see [`from_file`](Self::from_file) for streaming multi-reader mode with bounded memory.
    pub fn from_file_mmap(path: impl AsRef<Path>) -> Result<Self> {
        File::open(&path)?;
        Ok(Self::with_source(ParallelSource::Mmap(
            path.as_ref().to_path_buf(),
        )))
    }

    /// Creates a parallel parser for an arbitrary streaming reader.
    ///
    /// Non-seekable streaming readers currently execute in single-threaded streaming mode.
    #[inline]
    pub fn from_reader<Reader: Read + Send + 'a>(reader: Reader) -> Self {
        Self::with_source(ParallelSource::Reader(Box::new(reader)))
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

    /// Executes parallel processing across worker threads without reduction.
    ///
    /// Distributes chunks of the input to worker threads and allows the caller
    /// to drive the record loop directly on each worker chunk's [`FastxParser`].
    /// The closure can return `()` (infallible) or `Result<()>` (enabling `event?;` error propagation).
    ///
    /// If any worker thread encounters an error, sibling worker threads abort early
    /// and the error is returned to the caller.
    ///
    /// # Examples
    ///
    /// Simple record iteration with error propagation:
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    ///
    /// const CONFIG: Config = ParserOptions::default().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .for_each(|mut parser| {
    ///         while let Some(event) = parser.next() {
    ///             event?;
    ///             let header = parser.get_header();
    ///             let seq = parser.get_dna_string();
    ///             // process record...
    ///         }
    ///         Ok(())
    ///     })?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// Collecting per-record lengths into a shared thread-safe buffer:
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    /// use std::sync::Mutex;
    ///
    /// const CONFIG: Config = ParserOptions::default().ignore_headers().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let dest = Mutex::new(Vec::new());
    ///
    /// ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .for_each(|mut parser| {
    ///         let mut local = Vec::new();
    ///         while let Some(event) = parser.next() {
    ///             event?;
    ///             local.push(parser.get_dna_len());
    ///         }
    ///         dest.lock().unwrap().extend(local);
    ///         Ok(())
    ///     })?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn for_each<Process, Res>(self, process: Process) -> Result<()>
    where
        Res: ParallelResult,
        Process: Fn(FastxParser<'_, CONFIG>) -> Res + Sync,
    {
        self.reduce(move |parser| process(parser).into_result(), |(), ()| ())
    }

    /// Executes parallel processing across worker threads, combining thread results with a reduction function.
    ///
    /// The caller drives the parser loop on each worker chunk's [`FastxParser`], returning
    /// `Result<State>` (enabling `event?;` error propagation). The reduction `State` type is
    /// inferred automatically from the return type of `process`, and must implement [`Send`] and [`Default`].
    ///
    /// If any worker thread encounters an error, sibling worker threads abort early and the error
    /// is returned to the caller.
    ///
    /// # Examples
    ///
    /// Summing sequence lengths across worker threads:
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
    ///             while let Some(event) = parser.next() {
    ///                 event?;
    ///                 len += parser.get_dna_len();
    ///             }
    ///             Ok(len)
    ///         },
    ///         |acc, val| acc + val,
    ///     )?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// Collecting per-record lengths in parallel:
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
    ///             while let Some(event) = parser.next() {
    ///                 event?;
    ///                 local.push(parser.get_dna_len());
    ///             }
    ///             Ok(local)
    ///         },
    ///         |mut acc, chunk| {
    ///             acc.extend(chunk);
    ///             acc
    ///         },
    ///     )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn reduce<State, Process, Reduce>(self, process: Process, reduce: Reduce) -> Result<State>
    where
        State: Send + Default,
        Process: Fn(FastxParser<'_, CONFIG>) -> Result<State> + Sync,
        Reduce: Fn(State, State) -> State + Sync,
    {
        match self.source {
            ParallelSource::Slice(slice) => {
                run_slice::<CONFIG, State, _, _>(slice, self.threads, process, reduce)
            }
            ParallelSource::Mmap(path) => {
                let mmap = unsafe { memmap2::Mmap::map(&File::open(&path)?)? };
                run_slice::<CONFIG, State, _, _>(&mmap, self.threads, process, reduce)
            }
            ParallelSource::File(path) => {
                run_file::<CONFIG, State, _, _>(&path, self.threads, self.chunks, process, reduce)
            }
            ParallelSource::Reader(reader) => {
                if self.threads > 1 {
                    eprintln!(
                        "Warning: ParallelFastxParser currently processes streaming readers in single-threaded mode."
                    );
                }
                process(FastxParser::<CONFIG>::from_reader(reader)?)
            }
        }
    }

    /// Executes parallel processing across worker threads using a Paraseq [`ParallelProcessor`](paraseq::parallel::ParallelProcessor).
    ///
    /// Drives worker threads via Helicase's lock-free work-stealing engine, passing each
    /// record as a zero-copy [`paraseq::Record`] reference into [`process_record`](paraseq::parallel::ParallelProcessor::process_record).
    #[cfg(feature = "paraseq")]
    pub fn process_parallel<P>(self, processor: &mut P) -> Result<()>
    where
        P: for<'b, 'c> paraseq::parallel::ParallelProcessor<&'b FastxParser<'c, CONFIG>> + Sync,
    {
        let proc = processor.clone();
        self.for_each(move |mut parser| -> Result<()> {
            let mut local_proc = proc.clone();
            while let Some(event) = parser.next() {
                event?;
                local_proc
                    .process_record(&parser)
                    .map_err(|e| std::io::Error::other(e.to_string()))?;
            }
            local_proc
                .on_thread_complete()
                .map_err(|e| std::io::Error::other(e.to_string()))?;
            Ok(())
        })
    }
}
