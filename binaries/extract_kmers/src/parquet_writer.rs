use crate::ArrowDispatch;
use arrow::array::UInt64Array;
use arrow::record_batch::RecordBatch;
use helicase::*;
use parquet::arrow::ArrowWriter;
use std::fs::File;
use std::io;
use std::sync::Arc;

pub fn write_mer_chunk_to_parquet<const K: usize, B: ArrowDispatch, W: Send + io::Write>(
    mer_chunk: &mut MerChunk<K, B>,
    writer: &mut ArrowWriter<W>,
) -> io::Result<()> {
    let hash = Arc::new(UInt64Array::from_iter_values(
        mer_chunk.iter().map(|k| k.hash()),
    ));
    // Take ownership of the Vecs for zero-copy
    let high_vec = std::mem::take(&mut mer_chunk.0);
    let low_vec = std::mem::take(&mut mer_chunk.1);
    let high_array = B::to_arrow_array(high_vec);
    let low_array = B::to_arrow_array(low_vec);
    let batch = RecordBatch::try_new(B::build_schema(), vec![hash, high_array, low_array])
        .map_err(io::Error::other)?;
    writer.write(&batch).map_err(io::Error::other)?;
    Ok(())
}

pub fn mer_chunk_to_parquet<const K: usize, B: ArrowDispatch>(
    path: &str,
    mer_chunk: &mut MerChunk<K, B>,
) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer =
        ArrowWriter::try_new(file, B::build_schema(), None).map_err(io::Error::other)?;
    write_mer_chunk_to_parquet::<K, B, _>(mer_chunk, &mut writer)?;
    writer.close().map_err(io::Error::other)?;
    Ok(())
}

pub fn fastx_slice_to_parquet<const K: usize, B: ArrowDispatch>(
    path: &str,
    data: &[u8],
) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, B::build_schema(), None)?;
    let closure = |m: &mut MerChunk<K, B>| -> io::Result<()> {
        write_mer_chunk_to_parquet::<K, B, _>(m, &mut writer)?;
        Ok(())
    };
    chunk_process_from_fastx_slice::<K, B>(data, closure)?;
    writer.close()?;
    Ok(())
}
