use crate::TraceTimer;
use crate::arrow_builders::*;
use crate::arrow_types::*;
use arrow::ipc::reader::FileReader;
use arrow::ipc::writer::FileWriter;
use helicase::kmer::*;
use helicase::kmer_collection::*;
use helicase::*;
use rayon::prelude::*;
use std::fmt::Display;
use std::fs::File;
use std::path::{Path, PathBuf};
use tracing::{debug, info, instrument};

use eyre::{Result, eyre};

fn update_builders_from_chunk<const K: usize, B: BitStorage + ArrowDispatch>(
    builders: &mut [BucketBuilder<B>],
    mer_chunk: &mut MerChunk<K, B>,
    bucket_log_nb: usize,
    offset: usize,
) {
    let bucket_nb: usize = 1 << bucket_log_nb;
    let split_chunks = mer_chunk.par_split_by_keys(bucket_nb, |mer: Mer<K, B>| {
        let shift = offset - bucket_log_nb;
        (mer.hash() >> shift) as usize
    });
    for (i, chunk) in split_chunks.into_iter().enumerate() {
        builders[i].append_merchunk(chunk).unwrap();
    }
    mer_chunk.truncate(0);
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
/// The output format is a collection of Arrow IPC file or a concatenated version of that
///
/// The algorithm will produce bucket 2^bucket_log_nb buckets and dispatch kmers in each
/// buckets according to the bucket_log_nb most significant bits. By doing so, we
/// have the guarantee that each kmer in bucket i is smaller in hash order than
/// kmers in bucket j with i < j.
///
/// Optionally, each bucket can be :
///     - sort
///     - deduplicated with only one occurence of each kmer
///     - add a count value on the output with the count of each kmer /not implemented\
///   
///
/// Finally, it is possible to concatenate each bucket into one very big arrow IPC file.
///
/// The bucketing use the hash of the kmer. For a kmer of size K only the first min(64, 2K) bits
/// of this hash are used. We bucket by MSB first so we have a first offset bits are not used.
/// We split among the bucket using `bucket_log_nb` bits after this offset.
///
///
///  Recall that hash is 64-bit word:
///
///  ┌──────────────────────────────────────────────────────────────┐
///  │        offset bits  │ bucket_log_nb  bits │   remaining      │
///  │     (unused MSB)    │   (bucket index)    │      bits        │
///  ├─────────────────────┼─────────────────────┼──────────────────┤
///  │ 0 0 0 ... 0         │ b b b ... b         │ r r r ... r      │
///  └──────────────────────────────────────────────────────────────┘
///   ↑                   ↑
///   │                   └── used to select bucket
///   └── ignored / not used
///
pub fn fastx_slice_to_arrow<const K: usize, B: BitStorage + ArrowDispatch + Display + Sync>(
    data: &[u8],
    output_path: &str,
    bucket_log_nb: usize,
    bucket_capacity: usize,
    agg_mer_log_size_factor: u32,
    dedup: bool,
    concatenate_bucket: bool,
) -> Result<()> {
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

    let mut schema = B::build_schema(false);
    debug!(
        datatype = ?B::ARROW_DATA_TYPE,
        "arrow schema initialized"
    );
    let _init_span = tracing::info_span!("initialize_buckets");
    let mut builders: Vec<BucketBuilder<B>> = Vec::with_capacity(bucket_nb);

    let offset = if K < 64 { 2 * K } else { 64 };
    for i in 0..bucket_nb {
        let filename = Path::new(&dir_path).join(format!("{i}.arrow"));
        builders.push(BucketBuilder::new(
            &filename,
            i,
            (64 - offset) + bucket_log_nb,
            agg_mer_log_size_factor,
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
            let _time = TraceTimer::skippable("Computing merchunk", 1000);
            mer_chunk.append_from_columnar(seq);
        }
        kmer_count += mer_chunk.len();
        if mer_chunk.len() > bucket_capacity {
            update_builders_from_chunk(&mut builders, &mut mer_chunk, bucket_log_nb, offset);
        }
    }
    // do not forgot the final buffer !
    update_builders_from_chunk(&mut builders, &mut mer_chunk, bucket_log_nb, offset);

    for builder in builders.iter_mut() {
        builder.finish()?;
    }
    builders
        .par_iter_mut()
        .try_for_each(|builder| builder.clean::<K>(dedup))?;

    info!(
        sequences = seq_count,
        kmers = kmer_count,
        flushes = flush_count,
        buckets = bucket_nb,
        "kmer extraction completed"
    );

    if concatenate_bucket {
        if dedup {
            schema = B::build_schema(true);
        }
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
