use crate::TraceTimer;
use crate::arrow_types::*;
use arrow::array::*;
use arrow::compute::concat_batches;
use arrow::datatypes::Schema;
use arrow::ipc::reader::FileReader;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;
use helicase::kmer::*;
use helicase::kmer_collection::*;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::thread::JoinHandle;
use tracing::{info, trace};

use eyre::Result;
type ArrayType<B> = <<B as ArrowDispatch>::Builder as Builder>::ArrayType;

pub struct BucketBuilder<B: ArrowDispatch> {
    writer: FileWriter<File>,
    index: usize,
    // number of high bits already consumed by upstream clustering
    hash_bit_offset: usize,
    agg_mer_log_size_factor: u32,
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
        hash_bit_offset: usize,
        agg_mer_log_size_factor: u32,
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
            hash_bit_offset,
            agg_mer_log_size_factor,
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

    pub fn clean<const K: usize>(&mut self, sort_dedup: bool) -> Result<()> {
        let _time = TraceTimer::new(format!("Cleaning arrow file {}", self.index));

        // Ensure everything is flushed

        if !self.is_closed {
            self.finish()?;
        }

        // Reopen file for reading
        let file = File::open(&self.filename)?;
        let reader = FileReader::try_new(&file, None)?;
        let mut total_rows: usize = 0;
        for batch in reader {
            total_rows += batch?.num_rows();
        }
        let reader = FileReader::try_new(file, None)?;
        let mut schema = reader.schema();
        let batch;
        if sort_dedup {
            let mut agg_map;
            {
                let log2_capacity: usize =
                    (63 - total_rows.leading_zeros() + self.agg_mer_log_size_factor) as usize;
                let _time = TraceTimer::new(format!(
                    "init of the MerGroupedMap {} with log_capacity {}",
                    self.index, log2_capacity
                ));
                agg_map = MerGroupedMap::<K, B, u16>::new(log2_capacity, self.hash_bit_offset);
            }
            for batch in reader {
                let batch = batch?;
                let hash_array = batch
                    .column(0)
                    .as_any()
                    .downcast_ref::<UInt64Array>()
                    .unwrap();

                let high_array = batch
                    .column(1)
                    .as_any()
                    .downcast_ref::<ArrayType<B>>()
                    .unwrap();

                let low_array = batch
                    .column(2)
                    .as_any()
                    .downcast_ref::<ArrayType<B>>()
                    .unwrap();
                let hashes = hash_array.values();
                let highs = B::Builder::array_to_values(high_array);
                let lows = B::Builder::array_to_values(low_array);
                {
                    let _time = TraceTimer::new(format!("hash-sorting of bucket {}", self.index));
                    for i in 0..hashes.len() {
                        agg_map.add_hash(
                            hashes[i],
                            Mer::<K, B>(highs[i].into(), lows[i].into()),
                            1,
                        );
                    }
                }
            }
            {
                let ratio = (agg_map.len() as f64) / ((1 << agg_map.log2_capacity) as f64);

                let _time = TraceTimer::new(format!(
                    "Compactify {} of size {} (fill-ratio: {})",
                    self.index,
                    agg_map.len(),
                    ratio
                ));
                let coll_ratio = (agg_map.collision_count() as f64) / (agg_map.len() as f64);
                info!(
                    "Collision {} on {} collision-ratio:{} max-run: {}, bucket {}",
                    agg_map.collision_count(),
                    agg_map.len(),
                    coll_ratio,
                    agg_map.max_run(),
                    self.index
                );
                let (hashes, merchunk, values) = agg_map.compactify();
                schema = B::build_schema(true);
                let (high, low) = merchunk.as_slice().as_slices();
                let hash_array = UInt64Array::from(hashes);
                let high_array = B::Builder::array_from_values(high);
                let low_array = B::Builder::array_from_values(low);
                let count_array = UInt16Array::from(values);
                batch = RecordBatch::try_new(
                    schema.clone(),
                    vec![
                        Arc::new(hash_array),
                        Arc::new(high_array),
                        Arc::new(low_array),
                        Arc::new(count_array),
                    ],
                )?;
            }
        } else {
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
            let hash_array = combined.column(0).clone();
            let high_array = combined.column(1).clone();
            let low_array = combined.column(2).clone();
            info!("Concat done for bucket {}", self.index);
            batch = RecordBatch::try_new(schema.clone(), vec![hash_array, high_array, low_array])?;
        }

        // Rewrite file with a single batch
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
