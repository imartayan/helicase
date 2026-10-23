use helicase::config::*;
use helicase::*;
use paraseq::Record;
use paraseq::parallel::ParallelProcessor;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

// Ignore headers to maximize sequence extraction throughput
const CONFIG: Config = ParserOptions::default().ignore_headers().config();

/// A Paraseq processor that collects sequence lengths across threads.
#[derive(Clone)]
struct ParaseqLenCollector {
    dest: Arc<Mutex<Vec<usize>>>,
    local_buf: Vec<usize>,
}

impl ParaseqLenCollector {
    fn new(dest: Arc<Mutex<Vec<usize>>>) -> Self {
        Self {
            dest,
            local_buf: Vec::with_capacity(1024),
        }
    }
}

impl<R: Record> ParallelProcessor<R> for ParaseqLenCollector {
    #[inline(always)]
    fn process_record(&mut self, record: R) -> Result<(), paraseq::Error> {
        // Uses Helicase's implementation of paraseq::Record!
        // `record.seq_raw()` returns the borrowed &[u8] sequence directly.
        self.local_buf.push(record.seq_raw().len());
        Ok(())
    }

    fn on_thread_complete(&mut self) -> Result<(), paraseq::Error> {
        if !self.local_buf.is_empty() {
            self.dest.lock().unwrap().extend_from_slice(&self.local_buf);
            self.local_buf.clear();
        }
        Ok(())
    }
}

fn collect_record_lengths(
    path: &Path,
    threads: usize,
    capacity_hint: usize,
) -> std::io::Result<Vec<usize>> {
    let dest = Arc::new(Mutex::new(Vec::with_capacity(capacity_hint)));
    let mut processor = ParaseqLenCollector::new(Arc::clone(&dest));

    // Uses Helicase's `process_parallel` helper to run Paraseq's ParallelProcessor directly:
    ParallelFastxParser::<CONFIG>::from_file(path)?
        .threads(threads)
        .process_parallel(&mut processor)?;

    let result = std::mem::take(&mut *dest.lock().unwrap());
    Ok(result)
}

fn main() -> std::io::Result<()> {
    let path_str = std::env::args().nth(1).expect("No input file given");
    let path = Path::new(&path_str);

    let file_size = std::fs::metadata(path)?.len();
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);

    let capacity_hint = (file_size / 300).max(1024) as usize;

    let t0 = Instant::now();
    let lengths = collect_record_lengths(path, threads, capacity_hint)?;
    let elapsed = t0.elapsed().as_secs_f64();

    let total_len: usize = lengths.iter().sum();
    let throughput = (file_size as f64 / 1e9) / elapsed;

    println!("Records:    {}", lengths.len());
    println!("Total len:  {total_len} bp");
    println!(
        "Throughput: {throughput:.2} GB/s ({:.2} ms)",
        elapsed * 1000.0
    );

    Ok(())
}
