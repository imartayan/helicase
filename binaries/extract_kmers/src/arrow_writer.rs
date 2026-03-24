use crate::TraceTimer;
use crate::arrow_types::*;
use crate::kernels::dedup_sorted_columns;
use arrow::compute::{SortColumn, concat_batches, lexsort_to_indices, take};
use arrow::datatypes::Schema;
use arrow::ipc::reader::FileReader;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;
use helicase::kmer::*;
use helicase::kmer_collection::*;
use helicase::*;
use rayon::prelude::*;
use std::fmt::Display;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread::JoinHandle;
use tracing::{debug, info, instrument, trace};

use eyre::{Result, eyre};

pub struct BucketBuilder<B: ArrowDispatch> {
    writer: FileWriter<File>,
    index: usize,
    filename: String,
    builder: B::Builder,
    capacity: usize,
    handle: Option<JoinHandle<FileWriter<File>>>,
    schema: Arc<Schema>,
    is_closed: bool,
}

impl<B: ArrowDispatch> BucketBuilder<B> {
    pub fn new(
        filename: &Path,
        index: usize,
        capacity: usize,
        schema: Arc<Schema>,
    ) -> Result<Self> {
        let file = File::create(filename)?;
        let writer = FileWriter::try_new(file, &schema)?;
        let builder = B::Builder::with_capacity(capacity);
        trace!(file = %filename.display(), "creating bucket");
        Ok(Self {
            writer,
            index,
            filename: filename.display().to_string(),
            builder,
            capacity,
            handle: None,
            schema,
            is_closed: false,
        })
    }

    pub fn append_merchunk<const K: usize>(&mut self, m: MerChunk<K, B>) -> Result<()> {
        self.builder.append(m);
        self.commit()?;
        Ok(())
    }

    fn commit(&mut self) -> Result<()> {
        if self.builder.len() >= self.capacity {
            {
                let _t = TraceTimer::skippable(format!("waiting writing {}", self.index), 2);
                self.join()?;
            }
            let batch = self.builder.finish(self.schema.clone())?;
            trace!(
                bucket = self.filename,
                columns = batch.num_columns(),
                "Writing batch"
            );
            // move writer into flush thread
            let mut writer = std::mem::replace(&mut self.writer, dummy_writer());
            self.handle = Some(std::thread::spawn(move || {
                writer.write(&batch).unwrap();
                writer
            }));
        }
        Ok(())
    }

    pub fn join(&mut self) -> Result<()> {
        if let Some(handle) = self.handle.take() {
            self.writer = handle
                .join()
                .map_err(|e| eyre::eyre!("writer thread panicked: {:?}", e))?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.join()?;
        if let Some(handle) = self.handle.take() {
            self.writer = handle
                .join()
                .map_err(|e| eyre::eyre!("writer thread panic at finishing: {:?}", e))?;
        }

        // flush remaining rows
        if self.builder.len() > 0 {
            let batch = self.builder.finish(self.schema.clone())?;
            self.writer.write(&batch)?;
        }
        self.writer.finish()?;
        self.is_closed = true;
        Ok(())
    }

    pub fn clean(&mut self, sort: bool, dedup: bool, count: bool) -> Result<()> {
        if dedup & !sort {
            return Err(eyre!("Impossible to deduplicate kmers without sorting"));
        }

        if count & !dedup {
            return Err(eyre!("Impossible to count kmers without deduplication"));
        }
        let _time = TraceTimer::new(format!("Cleaning arrow file {}", self.index));

        // Ensure everything is flushed

        if !self.is_closed {
            self.finish()?;
        }

        // Reopen file for reading
        let file = File::open(&self.filename)?;
        let reader = FileReader::try_new(file, None)?;

        let schema = reader.schema();

        // Collect all batches
        let mut batches = Vec::new();
        for batch in reader {
            batches.push(batch.unwrap());
        }

        if batches.is_empty() {
            return Ok(());
        }

        // Concatenate into a single batch
        let combined = concat_batches(&schema, &batches)?;
        let mut hash_array = combined.column(0).clone();
        let mut high_array = combined.column(1).clone();
        let mut low_array = combined.column(2).clone();
        info!("Concat done for bucket {}", self.index);

        if sort {
            let _time = TraceTimer::new(format!("Sorting ...{}", self.index));
            let sort_columns = vec![
                SortColumn {
                    values: hash_array.clone(),
                    options: None,
                },
                SortColumn {
                    values: high_array.clone(),
                    options: None,
                },
                SortColumn {
                    values: low_array.clone(),
                    options: None,
                },
            ];
            let indices = lexsort_to_indices(&sort_columns, None)?;
            hash_array = take(&hash_array, &indices, None)?;
            high_array = take(&high_array, &indices, None)?;
            low_array = take(&low_array, &indices, None)?;
            if dedup {
                let deduped = dedup_sorted_columns(&[
                    hash_array.clone(),
                    high_array.clone(),
                    low_array.clone(),
                ])?;
                hash_array = deduped[0].clone();
                high_array = deduped[1].clone();
                low_array = deduped[2].clone();
            }
            if count {
                todo!();
            }
        }

        // Rewrite file with a single batch
        let batch = RecordBatch::try_new(schema.clone(), vec![hash_array, high_array, low_array])?;
        let file = File::create(&self.filename)?;
        let mut writer = FileWriter::try_new(file, &schema)?;
        writer.write(&batch)?;
        writer.finish()?;
        Ok(())
    }
}

// dummy_writer is just a placeholder for mem::replace
fn dummy_writer() -> FileWriter<File> {
    // you can open /dev/null or create a temp file
    let f = File::create("/dev/null").unwrap();
    FileWriter::try_new(f, &Arc::new(Schema::empty())).unwrap()
}

#[instrument(
    level = "info",
    skip(data),
    fields(
        k = K,
        bucket_log_nb,
        bucket_nb = 1usize << bucket_log_nb,
        temp_dir,
        output_path,
    )
)]
/// Parse a slice and produce a kmer collection of the shape: (hash, high, low) where:
///     - hash is the u64 hash value
///     - high is the high bits of each kmer
///     - low it the low bits of each kmer
///
/// The output format is a collection of Arrow IPC file
///
/// The algorithm will produce bucket 2^bucket_log_nb buckets and dispatch kmers in each
/// buckets according to the bucket_log_nb most significant bits. By doing so, we
/// have the guarantee that each kmer in bucket i is smaller in hash order than
/// kmers in bucket j with i < j.
///
/// Optionally, each bucket can be :
///     - sort
///     - deduplicated with only one occurence of each kmer
///     - add a count value on the output with the count of each kmer
///   
/// Remark that count requires dedup which requires sort. Otherwise the function will return
/// an Error.
///
/// Finally, it is possible to concatenate each bucket into one very big arrow IPC file.
pub fn fastx_slice_to_arrow<const K: usize, B: BitStorage + ArrowDispatch + Display + Sync>(
    data: &[u8],
    output_path: &str,
    bucket_log_nb: usize,
    bucket_capacity: usize,
    sort: bool,
    dedup: bool,
    concatenate_bucket: bool,
) -> Result<()> {
    const COUNT: bool = false;
    let _t = TraceTimer::new("fastx_slice_to_arrow");
    let dir_path = if concatenate_bucket {
        PathBuf::from(&format!("{}_dir", output_path))
    } else {
        PathBuf::from(output_path)
    };
    if dir_path.exists() {
        return Err(eyre!("Path {output_path} exists"));
    } else {
        std::fs::create_dir_all(&dir_path)?;
    }
    assert!(bucket_log_nb <= 64);
    let bucket_nb: usize = 1 << bucket_log_nb;
    assert!(bucket_nb <= 4096);
    info!(
        bucket_nb,
        bucket_capacity, "initializing kmer extraction pipeline"
    );

    let schema = B::build_schema(false);
    debug!(
        datatype = ?B::ARROW_DATA_TYPE,
        "arrow schema initialized"
    );
    let _init_span = tracing::info_span!("initialize_buckets");
    let mut builders: Vec<BucketBuilder<B>> = Vec::with_capacity(bucket_nb);

    for i in 0..bucket_nb {
        let filename = Path::new(&dir_path).join(format!("{i}.arrow"));
        builders.push(BucketBuilder::new(
            &filename,
            i,
            bucket_capacity,
            schema.clone(),
        )?);
    }

    const DNA_STRING: Config = ParserOptions::default()
        .ignore_headers()
        .dna_columnar()
        .split_non_actg()
        .return_dna_chunk(true)
        .return_record(false)
        .config();
    let mut parser = FastxParser::<DNA_STRING>::from_slice(data)?;
    let mut seq_count: u64 = 0;
    let mut kmer_count: usize = 0;
    let flush_count: u64 = 0;
    let mut mer_chunk = MerChunk::<K, B>::new();
    while parser.next().is_some() {
        seq_count += 1;
        let seq = parser.get_dna_columnar();
        {
            let _time = TraceTimer::new("Computing merchunk");
            mer_chunk.append_from_columnar(seq);
        }
        kmer_count += mer_chunk.len();
        if mer_chunk.len() > bucket_capacity {
            let split_chunks;
            {
                let _time = TraceTimer::new("Splitting merchunk");
                split_chunks = mer_chunk.par_split_by_keys(bucket_nb, |mer: Mer<K, B>| {
                    (mer.hash() >> (64 - bucket_log_nb)) as usize
                });
            }
            {
                let _time = TraceTimer::new("writting merchunk");
                for (i, chunk) in split_chunks.into_iter().enumerate() {
                    builders[i].append_merchunk(chunk).unwrap();
                }
            }
            info!("Seq: {} and kmers:{}", seq_count, kmer_count);
            mer_chunk.truncate(0);
        }
    }
    for builder in builders.iter_mut() {
        builder.finish()?;
    }
    builders
        .par_iter_mut()
        .try_for_each(|builder| builder.clean(sort, dedup, COUNT))?;

    info!(
        sequences = seq_count,
        kmers = kmer_count,
        flushes = flush_count,
        buckets = bucket_nb,
        "kmer extraction completed"
    );

    if concatenate_bucket {
        let _t = TraceTimer::new("concatenating bucket");
        let output = Path::new(output_path);
        let mut writer = FileWriter::try_new(File::create(output)?, &schema)?;

        for i in 0..bucket_nb {
            info!("Concatenating {i}");
            let bucket_file = &dir_path.join(format!("{i}.arrow"));
            let file = File::open(bucket_file)?;
            let reader = FileReader::try_new(file, None)?;
            for batch in reader {
                writer.write(&batch?)?;
            }
        }

        writer.finish()?;
        std::fs::remove_dir_all(dir_path)?;
    }
    Ok(())
}
