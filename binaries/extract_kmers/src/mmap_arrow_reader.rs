use arrow::array::RecordBatch;
use arrow::buffer::Buffer;
use arrow::error::Result as ArrowResult;
use arrow::ipc::convert::fb_to_schema;
use arrow::ipc::reader::{FileDecoder, read_footer_length};
use arrow::ipc::{Block, root_as_footer};
use bytes::Bytes;
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

/// A zero-copy IPC file reader using memory mapping.
///
/// Each `RecordBatch` references the underlying mmap buffer directly.
pub struct MmapIpcFile {
    buffer: Buffer,
    decoder: FileDecoder,
    batches: Vec<Block>,
    pos: usize,
}

impl MmapIpcFile {
    /// Open an Arrow IPC file and prepare zero-copy batch access.
    pub fn open<P: AsRef<Path>>(path: P) -> ArrowResult<Self> {
        // memory-map the file
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let bytes = Bytes::from_owner(mmap);
        let buffer = Buffer::from(bytes);

        // parse footer
        let trailer_start = buffer.len() - 10;
        let footer_len = read_footer_length(buffer[trailer_start..].try_into().unwrap())?;
        let footer = root_as_footer(&buffer[trailer_start - footer_len..trailer_start]).unwrap();

        // construct schema + decoder
        let schema = Arc::new(fb_to_schema(footer.schema().unwrap()));
        let mut decoder = FileDecoder::new(schema.clone(), footer.version());

        // read dictionaries
        for block in footer.dictionaries().iter().flatten() {
            let block_len = block.bodyLength() as usize + block.metaDataLength() as usize;
            let data = buffer.slice_with_length(block.offset() as _, block_len);
            decoder.read_dictionary(block, &data)?;
        }

        // collect record batch blocks
        let batches = footer
            .recordBatches()
            .map(|v| v.iter().copied().collect())
            .unwrap_or_default();

        Ok(Self {
            buffer,
            decoder,
            batches,
            pos: 0,
        })
    }

    /// Number of batches in the file
    pub fn num_batches(&self) -> usize {
        self.batches.len()
    }

    /// Returns the next batch or `None` if done
    pub fn next_batch(&mut self) -> ArrowResult<Option<RecordBatch>> {
        if self.pos >= self.batches.len() {
            return Ok(None);
        }
        let block = &self.batches[self.pos];
        self.pos += 1;

        let block_len = block.bodyLength() as usize + block.metaDataLength() as usize;
        let data = self
            .buffer
            .slice_with_length(block.offset() as _, block_len);
        self.decoder.read_record_batch(block, &data)
    }
}

impl Iterator for MmapIpcFile {
    type Item = ArrowResult<RecordBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_batch() {
            Ok(Some(batch)) => Some(Ok(batch)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
