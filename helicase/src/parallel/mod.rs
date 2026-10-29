//! Parallel multi-chunk parsing.
//!
//! Provides work-stealing execution across multiple threads,
//! for slices, files, and streams decompressed by a single thread.

use crate::config::*;
use crate::input::*;
use crate::parser::*;
use core::cell::RefCell;
use core::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
#[cfg(feature = "paraseq")]
use std::sync::Condvar;
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::sync::{Mutex, OnceLock};

/// Target bytes per chunk when a compressed stream is decompressed by a single thread
const STREAM_CHUNK_LEN: usize = 4 << 20;

/// Maximum bytes per chunk in [`ParallelFastxParser::process_parallel`], where a chunk is a batch
#[cfg(feature = "paraseq")]
const PARASEQ_MAX_CHUNK_LEN: usize = 16 << 20;

/// Callbacks of the worker threads: `init` creates the state of a worker from its index,
/// `process` parses a chunk given its index in stream order, `finish` turns the state into the
/// result of the worker. Running a job returns the results of the workers, in worker order.
struct Job<const CONFIG: Config, Init, Process, Finish> {
    init: Init,
    process: Process,
    finish: Finish,
}

/// Coordinates worker threads so `on_batch_complete` runs in stream order.
#[cfg(feature = "paraseq")]
struct OrderGate {
    state: Mutex<usize>,
    cond: Condvar,
}

/// Worker state of [`ParallelFastxParser::process_parallel`].
/// Poisons the gate when dropped before `finish`, so an error or a panic cannot block other workers.
#[cfg(feature = "paraseq")]
struct ParaseqWorker<'a, P> {
    processor: P,
    gate: Option<&'a OrderGate>,
}

#[cfg(feature = "paraseq")]
impl OrderGate {
    const POISONED: usize = usize::MAX;

    fn new() -> Self {
        Self {
            state: Mutex::new(0),
            cond: Condvar::new(),
        }
    }

    fn wait_turn(&self, ticket: usize) {
        let mut cur = self.state.lock().unwrap();
        while *cur != ticket && *cur != Self::POISONED {
            cur = self.cond.wait(cur).unwrap();
        }
    }

    fn advance(&self, ticket: usize) {
        let mut cur = self.state.lock().unwrap();
        if *cur != Self::POISONED && ticket + 1 > *cur {
            *cur = ticket + 1;
            self.cond.notify_all();
        }
    }

    fn poison(&self) {
        let mut cur = self.state.lock().unwrap();
        *cur = Self::POISONED;
        self.cond.notify_all();
    }
}

#[cfg(feature = "paraseq")]
impl<P> Drop for ParaseqWorker<'_, P> {
    fn drop(&mut self) {
        if let Some(gate) = self.gate {
            gate.poison();
        }
    }
}

#[cfg(any(feature = "gz", feature = "zstd"))]
mod frames;
#[cfg(feature = "zstd")]
pub(crate) use frames::index_zstd_frames;
#[cfg(test)]
mod tests;

impl<const CONFIG: Config, T, State, Init, Process, Finish> Job<CONFIG, Init, Process, Finish>
where
    State: Send,
    Init: Fn(usize) -> T + Sync,
    Process: Fn(&mut T, usize, FastxParser<'_, CONFIG>) -> Result<()> + Sync,
    Finish: Fn(T) -> Result<State> + Sync,
{
    fn run(&self, parser: ParallelFastxParser<'_, CONFIG>) -> Result<Vec<State>> {
        let (threads, chunks, max_chunk_len) =
            (parser.threads, parser.chunks, parser.max_chunk_len);
        match parser.source {
            ParallelSource::Slice(slice) => self.run_slice(slice, threads, chunks, max_chunk_len),
            ParallelSource::Mmap(path) => {
                let mmap = unsafe { memmap2::Mmap::map(&File::open(&path)?)? };
                self.run_slice(&mmap, threads, chunks, max_chunk_len)
            }
            ParallelSource::File(path) => self.run_file(&path, threads, chunks, max_chunk_len),
            ParallelSource::Reader(reader) => {
                self.run_stream(decompress(reader)?, threads, max_chunk_len)
            }
        }
    }

    fn run_single(&self, parser: Result<FastxParser<'_, CONFIG>>) -> Result<Vec<State>> {
        let mut state = (self.init)(0);
        (self.process)(&mut state, 0, parser?)?;
        Ok(vec![(self.finish)(state)?])
    }

    fn run_workers(
        &self,
        threads: usize,
        work: impl Fn(&mut T, usize) -> Option<Result<()>> + Sync,
    ) -> Result<Vec<State>> {
        let next_ticket = &AtomicUsize::new(threads);
        let stop = &AtomicBool::new(false);
        let work = &work;
        std::thread::scope(|scope| {
            // Collected to spawn all the threads before joining them.
            #[allow(clippy::needless_collect)]
            let handles: Vec<_> = (0..threads)
                .map(|worker| {
                    scope.spawn(move || -> Result<State> {
                        let mut state = (self.init)(worker);
                        let mut ticket = worker;
                        while !stop.load(Ordering::Relaxed)
                            && let Some(res) = work(&mut state, ticket)
                        {
                            res.inspect_err(|_| stop.store(true, Ordering::Relaxed))?;
                            ticket = next_ticket.fetch_add(1, Ordering::Relaxed);
                        }
                        (self.finish)(state).inspect_err(|_| stop.store(true, Ordering::Relaxed))
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        })
    }

    /// Parses `num_chunks` chunks delimited by `split`, an empty chunk gives an empty parser.
    fn run_chunks<P>(
        &self,
        num_chunks: usize,
        threads: usize,
        (first, last): (P, P),
        split: impl Fn(usize) -> P + Sync,
        parse_chunk: impl Fn(&mut T, usize, P, P) -> Result<()> + Sync,
    ) -> Result<Vec<State>>
    where
        P: Copy + Send + Sync,
    {
        let splits: Vec<OnceLock<P>> = (0..=num_chunks).map(|_| OnceLock::new()).collect();
        let _ = splits[0].set(first);
        let _ = splits[num_chunks].set(last);
        let split_at = |i: usize| *splits[i].get_or_init(|| split(i));
        self.run_workers(threads.min(num_chunks), |state, i| {
            (i < num_chunks).then(|| parse_chunk(state, i, split_at(i), split_at(i + 1)))
        })
    }

    fn run_stream(
        &self,
        reader: impl Read + Send,
        threads: usize,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        let chunk_len = STREAM_CHUNK_LEN.min(max_chunk_len);
        let (chunk_sender, chunks) = sync_channel(threads);
        let (free_sender, free) = channel();
        let chunks = Mutex::new(chunks);
        // Dropped when the workers are done, which stops the reader thread
        let work = move |state: &mut T, _| {
            let (ticket, mut chunk): (usize, Vec<u8>) = chunks.lock().unwrap().recv().ok()?;
            let res =
                FastxParser::from_slice(&chunk).and_then(|p| (self.process)(state, ticket, p));
            chunk.clear();
            chunk.shrink_to(2 * chunk_len);
            let _ = free_sender.send(chunk);
            Some(res)
        };
        std::thread::scope(|scope| {
            let sent = scope.spawn(move || send_chunks(reader, chunk_len, free, chunk_sender));
            let res = self.run_workers(threads, work);
            sent.join().unwrap()?;
            res
        })
    }

    fn run_slice(
        &self,
        data: &[u8],
        threads: usize,
        custom_chunks: Option<usize>,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        if (threads == 1 && data.len() <= max_chunk_len) || data.len() < 32 * 1024 {
            return self.run_single(FastxParser::from_slice(data));
        }
        let format = Format::try_from_byte(data[0])?;
        let num_chunks = custom_chunks
            .unwrap_or(threads * 4)
            .max(data.len().div_ceil(max_chunk_len))
            .min(data.len() / (16 * 1024));
        let chunk_len = data.len() / num_chunks;
        self.run_chunks(
            num_chunks,
            threads,
            (0, data.len()),
            |i| {
                format
                    .find_record_start(data, i * chunk_len)
                    .unwrap_or(data.len())
            },
            |state, ticket, start, end| {
                (self.process)(state, ticket, FastxParser::from_slice(&data[start..end])?)
            },
        )
    }

    fn run_file(
        &self,
        path: &Path,
        threads: usize,
        custom_chunks: Option<usize>,
        max_chunk_len: usize,
    ) -> Result<Vec<State>> {
        let compression = CompressionFormat::from_path(path)?;
        match compression {
            None => {}
            #[cfg(feature = "gz")]
            Some(CompressionFormat::Bgzf) if threads > 1 => {
                let mmap = unsafe { memmap2::Mmap::map(&File::open(path)?)? };
                let blocks = frames::index_bgzf_blocks(&mmap);
                return self.run_frames::<frames::BgzfDecoder>(
                    &mmap,
                    &blocks,
                    threads,
                    custom_chunks,
                    max_chunk_len,
                );
            }
            #[cfg(feature = "zstd")]
            Some(CompressionFormat::Pzstd) if threads > 1 => {
                let mmap = unsafe { memmap2::Mmap::map(&File::open(path)?)? };
                let zstd_frames = frames::index_zstd_frames(&mmap);
                if zstd_frames.len() > 1 {
                    return self.run_frames::<frames::ZstdDecoder>(
                        &mmap,
                        &zstd_frames,
                        threads,
                        custom_chunks,
                        max_chunk_len,
                    );
                }
            }
            Some(_) => {}
        }
        if compression.is_some() {
            // The read-ahead of `from_file` is about 3% faster than `run_stream` on one thread
            if threads == 1 {
                return self.run_single(FastxParser::from_file(path));
            }
            return self.run_stream(decompress(File::open(path)?)?, threads, max_chunk_len);
        }

        let mut file = File::open(path)?;
        let file_len = file.metadata()?.len();
        if (threads == 1 && file_len <= max_chunk_len as u64) || file_len < 32 * 1024 {
            return self.run_single(FastxParser::from_file(path));
        }
        let mut first_byte = [0u8; 1];
        file.read_exact(&mut first_byte)?;
        let format = Format::try_from_byte(first_byte[0])?;
        let num_chunks = custom_chunks
            .unwrap_or(threads * 4)
            .max(file_len.div_ceil(max_chunk_len as u64) as usize)
            .min((file_len / (16 * 1024)) as usize);
        let chunk_len = file_len / num_chunks as u64;
        self.run_chunks(
            num_chunks,
            threads,
            (0, file_len),
            |i| {
                let pos = i as u64 * chunk_len;
                with_thread_file(path, |f| probe_file_record_start(f, pos, file_len, format))
                    .unwrap_or(pos)
            },
            |state, ticket, start, end| {
                with_thread_file(path, |file| {
                    file.seek(SeekFrom::Start(start))?;
                    let parser = FastxParser::from_reader(file.take(end - start))?;
                    (self.process)(state, ticket, parser)
                })
            },
        )
    }
}

/// Sends the content of `reader` in chunks of `chunk_len` bytes or more, cut before a record start.
/// Buffers are reused from `free`, and stops when `chunks` is closed.
fn send_chunks(
    mut reader: impl Read,
    chunk_len: usize,
    free: Receiver<Vec<u8>>,
    chunks: SyncSender<(usize, Vec<u8>)>,
) -> Result<()> {
    let mut carry = Vec::new();
    let mut ticket = 0;
    loop {
        let mut buf = free.try_recv().unwrap_or_default();
        buf.append(&mut carry);
        let from = buf.len() + chunk_len;
        // Read ahead to find the next record start, doubling for long records
        let mut slack = chunk_len.div_ceil(16);
        let (cut, eof) = loop {
            let missing = from + slack - buf.len();
            if (&mut reader).take(missing as u64).read_to_end(&mut buf)? < missing {
                break (buf.len(), true);
            }
            if let Some(cut) = Format::try_from_byte(buf[0])?.find_record_start(&buf, from) {
                break (cut, false);
            }
            slack *= 2;
        };
        carry.extend_from_slice(&buf[cut..]);
        buf.truncate(cut);
        if buf.is_empty() || chunks.send((ticket, buf)).is_err() || eof {
            return Ok(());
        }
        ticket += 1;
    }
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

/// Runs `f` on `path` opened once per worker (reopening it per chunk costs 5-10% at 8 threads).
fn with_thread_file<T>(path: &Path, f: impl FnOnce(&mut File) -> Result<T>) -> Result<T> {
    thread_local!(static FILE: RefCell<Option<File>> = const { RefCell::new(None) });
    FILE.with_borrow_mut(|file| match file {
        Some(file) => f(file),
        None => f(file.insert(File::open(path)?)),
    })
}

/// Trait for return types of [`ParallelFastxParser::for_each`] and [`ParallelFastxParser::for_each_with`] closures.
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
/// or [`for_each_with`](ParallelFastxParser::for_each_with).
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
    max_chunk_len: usize,
}

impl<'a, const CONFIG: Config> ParallelFastxParser<'a, CONFIG> {
    fn with_source(source: ParallelSource<'a>) -> Self {
        Self {
            source,
            threads: std::thread::available_parallelism().map_or(1, |n| n.get()),
            chunks: None,
            max_chunk_len: usize::MAX,
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
    /// The reader is decompressed by a single thread, which distributes chunks to the workers.
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

    /// Sets the number of chunks to split the input into, for work stealing.
    ///
    /// Defaults to `threads * 4` (`threads` for BGZF and multi-frame zstd).
    /// More chunks balance uneven work better, fewer chunks cost less overhead.
    ///
    /// Ignored with one thread, and on streams (`from_reader`, plain gzip, single-frame zstd)
    /// which are cut in chunks of 4 MiB.
    #[inline]
    pub fn chunks(mut self, num_chunks: usize) -> Self {
        self.chunks = Some(num_chunks.max(1));
        self
    }

    /// Executes parallel processing across worker threads.
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
        self.for_each_with(|| (), |(), parser| process(parser))
            .map(drop)
    }

    /// Executes parallel processing across worker threads with a state per thread, and returns the states.
    ///
    /// `init` creates the state of each worker thread, and `process` parses each chunk of the input
    /// into the state of its thread, as in [`for_each`](Self::for_each).
    /// Returns one state per worker thread: between 1 and `threads` states, fewer on small inputs,
    /// and a state is left as created if its thread got no records.
    /// Which records end up in which state changes between runs, so only a merge that does not
    /// depend on the order (sum, union, counts) gives a stable result.
    ///
    /// If any worker thread encounters an error, sibling worker threads abort early and the error
    /// is returned to the caller.
    ///
    /// # Examples
    ///
    /// Counting records and bases:
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    ///
    /// const CONFIG: Config = ParserOptions::default().ignore_headers().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let counts = ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .for_each_with(
    ///         || (0, 0),
    ///         |(records, len), mut parser| {
    ///             while let Some(event) = parser.next() {
    ///                 event?;
    ///                 *records += 1;
    ///                 *len += parser.get_dna_len();
    ///             }
    ///             Ok(())
    ///         },
    ///     )?;
    /// let records: usize = counts.iter().map(|c| c.0).sum();
    /// let len: usize = counts.iter().map(|c| c.1).sum();
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// Mapping headers to sequence lengths, with one table per worker thread:
    /// ```rust,no_run
    /// use helicase::config::*;
    /// use helicase::*;
    /// use std::collections::HashMap;
    ///
    /// const CONFIG: Config = ParserOptions::default().config();
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let tables = ParallelFastxParser::<CONFIG>::from_file("reads.fastq")?
    ///     .threads(8)
    ///     .for_each_with(HashMap::new, |table, mut parser| {
    ///         while let Some(event) = parser.next() {
    ///             event?;
    ///             let header = String::from_utf8_lossy(parser.get_header()).into_owned();
    ///             table.insert(header, parser.get_dna_len());
    ///         }
    ///         Ok(())
    ///     })?;
    /// let lengths: HashMap<String, usize> = tables.into_iter().flatten().collect();
    /// # Ok(())
    /// # }
    /// ```
    pub fn for_each_with<State, Init, Process, Res>(
        self,
        init: Init,
        process: Process,
    ) -> Result<Vec<State>>
    where
        State: Send,
        Init: Fn() -> State + Sync,
        Res: ParallelResult,
        Process: Fn(&mut State, FastxParser<'_, CONFIG>) -> Res + Sync,
    {
        Job {
            init: |_| init(),
            process: |state: &mut State, _, parser: FastxParser<'_, CONFIG>| {
                process(state, parser).into_result()
            },
            finish: |state| Ok(state),
        }
        .run(self)
    }

    /// Executes parallel processing across worker threads, merging the results of the chunks with `reduce`.
    ///
    /// Creates a state per chunk, [`for_each_with`](Self::for_each_with) keeps one per thread instead.
    #[deprecated(
        since = "1.2.0",
        note = "use `for_each_with`, which keeps one state per thread instead of one per chunk"
    )]
    pub fn reduce<State, Process, Reduce>(self, process: Process, reduce: Reduce) -> Result<State>
    where
        State: Send + Default,
        Process: Fn(FastxParser<'_, CONFIG>) -> Result<State> + Sync,
        Reduce: Fn(State, State) -> State + Sync,
    {
        let states = self.for_each_with(
            || None,
            |acc: &mut Option<State>, parser| {
                let res = process(parser)?;
                *acc = Some(match acc.take() {
                    Some(acc) => reduce(acc, res),
                    None => res,
                });
                Ok(())
            },
        )?;
        Ok(states
            .into_iter()
            .flatten()
            .reduce(reduce)
            .unwrap_or_default())
    }

    /// Executes parallel processing across worker threads using a Paraseq [`ParallelProcessor`](paraseq::parallel::ParallelProcessor).
    ///
    /// Drives worker threads via Helicase's lock-free work-stealing engine, passing each
    /// record as a zero-copy [`paraseq::Record`] reference into [`process_record`](paraseq::parallel::ParallelProcessor::process_record).
    /// A chunk (at most 16 MiB of input) is a batch: after its last record, a worker calls
    /// [`on_batch_complete`](paraseq::parallel::ParallelProcessor::on_batch_complete).
    /// [`set_thread_id`](paraseq::parallel::ParallelProcessor::set_thread_id) receives the index of the worker thread (below `threads`),
    /// and [`on_thread_complete`](paraseq::parallel::ParallelProcessor::on_thread_complete) is called once per thread upon worker exit.
    ///
    /// If the processor requests ordering via [`requires_ordering`](paraseq::parallel::ParallelProcessor::requires_ordering)
    /// (or is wrapped in [`Ordered`](paraseq::parallel::Ordered)), batches are committed in stream order.
    #[cfg(feature = "paraseq")]
    pub fn process_parallel<P>(mut self, processor: &mut P) -> Result<()>
    where
        P: for<'b, 'c> paraseq::parallel::ParallelProcessor<&'b FastxParser<'c, CONFIG>> + Sync,
    {
        self.max_chunk_len = PARASEQ_MAX_CHUNK_LEN;
        let order_gate = processor.requires_ordering().then(OrderGate::new);
        let gate = order_gate.as_ref();
        let processor = &*processor;
        let to_io_error = |e: paraseq::Error| std::io::Error::other(e.to_string());
        Job {
            init: |worker_id| {
                let mut processor = processor.clone();
                processor.set_thread_id(worker_id);
                ParaseqWorker { processor, gate }
            },
            process: |worker: &mut ParaseqWorker<P>,
                      ticket,
                      mut parser: FastxParser<'_, CONFIG>| {
                let mut empty = true;
                while let Some(event) = parser.next() {
                    event?;
                    worker
                        .processor
                        .process_record(&parser)
                        .map_err(to_io_error)?;
                    empty = false;
                }
                if let Some(gate) = gate {
                    gate.wait_turn(ticket);
                }
                if !empty {
                    worker.processor.on_batch_complete().map_err(to_io_error)?;
                }
                if let Some(gate) = gate {
                    gate.advance(ticket);
                }
                Ok(())
            },
            finish: |mut worker: ParaseqWorker<P>| {
                worker.processor.on_thread_complete().map_err(to_io_error)?;
                worker.gate = None;
                Ok(())
            },
        }
        .run(self)
        .map(drop)
    }
}
