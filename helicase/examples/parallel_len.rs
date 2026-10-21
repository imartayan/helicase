use helicase::config::*;
use helicase::*;
use std::cell::RefCell;
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

// Ignore headers to maximize sequence extraction throughput
const CONFIG: Config = ParserOptions::default().ignore_headers().config();

thread_local! {
    static LOCAL_BUF: RefCell<Vec<usize>> = const { RefCell::new(Vec::new()) };
}

fn collect_record_lengths(
    path: &Path,
    threads: usize,
    capacity_hint: usize,
) -> std::io::Result<Vec<usize>> {
    let dest = Mutex::new(Vec::with_capacity(capacity_hint));

    ParallelFastxParser::<CONFIG>::from_file(path)?
        .threads(threads)
        .for_each(|mut parser| {
            LOCAL_BUF.with_borrow_mut(|buf| -> parser::Result<()> {
                buf.clear(); // Reuses allocated capacity across chunks on this thread
                while let Some(event) = parser.next() {
                    event?;
                    buf.push(parser.get_dna_len());
                }
                dest.lock().unwrap().extend_from_slice(buf);
                Ok(())
            })
        })?;

    Ok(dest.into_inner().unwrap())
}

fn main() -> std::io::Result<()> {
    let path_str = std::env::args().nth(1).expect("No input file given");
    let path = Path::new(&path_str);

    let file_size = std::fs::metadata(path)?.len();
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);

    // Estimate record count from file size to preallocate destination vector
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
