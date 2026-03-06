use arrow::array::{ArrayRef, BinaryArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array};
use arrow::datatypes::*;
use arrow::record_batch::RecordBatch;
use helicase::*;
use parquet::arrow::ArrowWriter;
use std::fs::File;
use std::io;
use std::sync::Arc;

pub trait ArrowDispatch: BitStorage {
    const ARROW_DATA_TYPE: DataType;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef;
    fn build_schema() -> Arc<Schema> {
        Arc::new(Schema::new(vec![
            Field::new("hash", DataType::UInt64, false),
            Field::new("high", Self::ARROW_DATA_TYPE, false),
            Field::new("low", Self::ARROW_DATA_TYPE, false),
        ]))
    }
}

/// Implementations for primitive types
impl ArrowDispatch for u8 {
    const ARROW_DATA_TYPE: DataType = DataType::UInt8;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef {
        Arc::new(UInt8Array::from_iter_values(slice.iter().map(|b| b.value))) as ArrayRef
    }
}

impl ArrowDispatch for u16 {
    const ARROW_DATA_TYPE: DataType = DataType::UInt16;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef {
        Arc::new(UInt16Array::from_iter_values(slice.iter().map(|b| b.value))) as ArrayRef
    }
}

impl ArrowDispatch for u32 {
    const ARROW_DATA_TYPE: DataType = DataType::UInt32;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef {
        Arc::new(UInt32Array::from_iter_values(slice.iter().map(|b| b.value))) as ArrayRef
    }
}

impl ArrowDispatch for u64 {
    const ARROW_DATA_TYPE: DataType = DataType::UInt64;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef {
        Arc::new(UInt64Array::from_iter_values(slice.iter().map(|b| b.value))) as ArrayRef
    }
}

/// Special implementation for u128 → store as binary
impl ArrowDispatch for u128 {
    const ARROW_DATA_TYPE: DataType = DataType::Binary;
    fn to_arrow_array<const K: usize>(slice: &[BitString<Self, K>]) -> ArrayRef {
        Arc::new(BinaryArray::from_iter_values(
            slice.iter().map(|b| b.value.to_le_bytes().to_vec()),
        )) as ArrayRef
    }
}

pub fn write_mer_slice_to_parquet<const K: usize, B: ArrowDispatch, W: Send + io::Write>(
    mer_slice: &MerSlice<K, B>,
    writer: &mut ArrowWriter<W>,
) -> io::Result<()> {
    let hash = Arc::new(UInt64Array::from_iter_values(
        mer_slice.iter().map(|k| k.hash()),
    ));
    let high_array = B::to_arrow_array(mer_slice.0);
    let low_array = B::to_arrow_array(mer_slice.1);
    let batch = RecordBatch::try_new(B::build_schema(), vec![hash, high_array, low_array])
        .map_err(io::Error::other)?;
    writer.write(&batch).map_err(io::Error::other)?;
    Ok(())
}

pub fn mer_slice_to_parquet<const K: usize, B: ArrowDispatch>(
    path: &str,
    mer_slice: &MerSlice<K, B>,
) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer =
        ArrowWriter::try_new(file, B::build_schema(), None).map_err(io::Error::other)?;
    write_mer_slice_to_parquet::<K, B, _>(mer_slice, &mut writer)?;
    writer.close().map_err(io::Error::other)?;
    Ok(())
}

pub fn fastx_slice_to_parquet<const K: usize, B: ArrowDispatch>(
    path: &str,
    data: &[u8],
) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, B::build_schema(), None)?;
    let closure = |m: &MerSlice<K, B>| -> io::Result<()> {
        write_mer_slice_to_parquet::<K, B, _>(m, &mut writer)?;
        Ok(())
    };
    chunk_process_from_fastx_slice::<K, B>(data, closure)?;
    writer.close()?;
    Ok(())
}
