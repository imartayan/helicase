use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use helicase::kmer::*;
use helicase::kmer_collection::MerChunk;

use paste::paste;
use std::sync::{Arc, OnceLock};

use eyre::Result;

pub trait Builder: Send + Sync {
    type Value: BitStorage;
    type ArrayType;
    const ARROW_DATA_TYPE: DataType;
    fn with_capacity(capacity: usize) -> Self;
    fn append<const K: usize>(&mut self, m: MerChunk<K, Self::Value>);
    fn finish(&mut self, schema: Arc<Schema>) -> Result<RecordBatch>;
    fn len(&self) -> usize;
}

/// ArrowDispatch trait
pub trait ArrowDispatch: BitStorage + Send + Sync {
    type Builder: Builder<Value = Self>;
    const ARROW_DATA_TYPE: DataType = Self::Builder::ARROW_DATA_TYPE;

    /// Schema helper
    fn build_schema(count: bool) -> Arc<Schema> {
        static SCHEMA: OnceLock<Arc<Schema>> = OnceLock::new();
        SCHEMA
            .get_or_init(|| {
                let mut columns = vec![
                    Field::new("hash", DataType::UInt64, false),
                    Field::new("high", Self::ARROW_DATA_TYPE, false),
                    Field::new("low", Self::ARROW_DATA_TYPE, false),
                ];
                if count {
                    columns.push(Field::new("count", DataType::UInt16, false));
                }
                Arc::new(Schema::new(columns))
            })
            .clone()
    }
}

macro_rules! impl_arrow_fixed {
    ($t:ty, $arrow:ty, $dt:expr, $atype: ty) => {
        paste! {
            pub struct [<Builder $t:upper>] {
                hash: Vec<u64>,
                high: Vec<$t>,
                low: Vec<$t>,
                capacity: usize,
            }
            impl Builder for [<Builder $t:upper>] {
                type Value = $t;
                type ArrayType = $atype;
                const ARROW_DATA_TYPE : DataType = $dt;
                fn with_capacity(capacity: usize) -> Self {
                    Self {
                        hash: Vec::<u64>::with_capacity(capacity),
                        high: Vec::<$t>::with_capacity(capacity),
                        low: Vec::<$t>::with_capacity(capacity),
                        capacity,
                    }
                }

                fn append<const K: usize>(&mut self, m: MerChunk<K, $t>) {
                    let mut hashes : Vec<u64> = m.compute_hashes().collect();
                    // unsafe:: bitstring is a transparent wrapper around $t so it it fine.
                    let high: &[$t] = unsafe { std::mem::transmute(m.0.as_slice()) };
                    let low: &[$t] = unsafe { std::mem::transmute(m.1.as_slice()) };
                    self.high.extend_from_slice(high);
                    self.low.extend_from_slice(low);
                    self.hash.append(&mut hashes);
                }


                fn finish(&mut self, schema: Arc<Schema>) -> Result<RecordBatch>{
                    let v0 = std::mem::replace(&mut self.hash, Vec::with_capacity(self.capacity));
                    let v1 = std::mem::replace(&mut self.high, Vec::with_capacity(self.capacity));
                    let v2 = std::mem::replace(&mut self.low, Vec::with_capacity(self.capacity));

                    let arr0 : UInt64Array = v0.into();
                    let arr1 : $atype = v1.into();
                    let arr2 : $atype = v2.into();
                    let arrays : Vec<Arc<dyn Array>> = vec![
                        Arc::new(arr0),
                        Arc::new(arr1),
                        Arc::new(arr2)
                    ];
                    let batch = RecordBatch::try_new(schema, arrays)?;
                    Ok(batch)
                }

                fn len(&self) -> usize{
                    self.hash.len()
                }

            }
            impl ArrowDispatch for $t {
                type Builder =[<Builder $t:upper>];
            }
        }
    };
}

impl_arrow_fixed!(u8, UInt8Type, DataType::UInt8, UInt8Array);
impl_arrow_fixed!(u16, UInt16Type, DataType::UInt16, UInt16Array);
impl_arrow_fixed!(u32, UInt32Type, DataType::UInt32, UInt32Array);
impl_arrow_fixed!(u64, UInt64Type, DataType::UInt64, UInt64Array);

pub struct BuilderU128 {
    hash: Vec<u64>,
    high: Vec<u128>,
    low: Vec<u128>,
    capacity: usize,
}
impl Builder for BuilderU128 {
    type Value = u128;
    type ArrayType = FixedSizeBinaryArray;
    const ARROW_DATA_TYPE: DataType = DataType::FixedSizeBinary(16);
    fn with_capacity(capacity: usize) -> Self {
        Self {
            hash: Vec::<u64>::with_capacity(capacity),
            high: Vec::<u128>::with_capacity(capacity),
            low: Vec::<u128>::with_capacity(capacity),
            capacity,
        }
    }

    fn append<const K: usize>(&mut self, m: MerChunk<K, u128>) {
        let mut hashes: Vec<u64> = m.compute_hashes().collect();
        self.hash.append(&mut hashes);
        // unsafe:: bitstring is a transparent wrapper around $t so it it fine.
        let high: &[u128] = unsafe { std::mem::transmute(m.0.as_slice()) };
        let low: &[u128] = unsafe { std::mem::transmute(m.1.as_slice()) };
        self.high.extend_from_slice(high);
        self.low.extend_from_slice(low);
    }

    fn finish(&mut self, schema: Arc<Schema>) -> Result<RecordBatch> {
        let v0 = std::mem::replace(&mut self.hash, Vec::with_capacity(self.capacity));
        let v1 = std::mem::replace(&mut self.high, Vec::with_capacity(self.capacity));
        let v2 = std::mem::replace(&mut self.low, Vec::with_capacity(self.capacity));

        let arr0: UInt64Array = v0.into();
        let vv1: Vec<_> = v1.into_iter().map(|e| e.to_le_bytes()).collect();
        let vvv1: Vec<_> = vv1.iter().collect();
        let vv2: Vec<_> = v2.into_iter().map(|e| e.to_le_bytes()).collect();
        let vvv2: Vec<_> = vv2.iter().collect();
        let arr1 = FixedSizeBinaryArray::from(vvv1);
        let arr2 = FixedSizeBinaryArray::from(vvv2);
        let arrays: Vec<Arc<dyn Array>> = vec![Arc::new(arr0), Arc::new(arr1), Arc::new(arr2)];
        let batch = RecordBatch::try_new(schema, arrays)?;
        Ok(batch)
    }

    fn len(&self) -> usize {
        self.hash.len()
    }
}
impl ArrowDispatch for u128 {
    type Builder = BuilderU128;
}
