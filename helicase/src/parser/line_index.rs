use crate::config::Config;
use crate::lexer::SCAN_MODE_HEADER;
use crate::simd::extract_bitmask;

// Short records are parsed from an index of line feeds instead of the block state machine,
// once enough consecutive records are shorter than `MAX_RECORD_LEN`.
const MAX_RECORD_LEN: usize = 1024;
const MIN_SHORT_RECORDS: usize = 16;
const WINDOW_BLOCKS: usize = 64;
const MAX_SCAN: usize = 4096;
const LINE_FEEDS_ONLY: Config = 0;

#[derive(Default)]
pub(crate) struct LineIndex {
    pub active: bool,
    pub next_start: usize,
    short_records: usize,
    // Entries that gave no record double the streak needed for the next one.
    backoff: u32,
    consumed: bool,
    scan_pos: usize,
    line_feeds: Vec<usize>,
    head: usize,
    len: usize,
}

impl LineIndex {
    // The next `N` line feeds from `next_start`, if they are within `MAX_SCAN` bytes
    // and before the last block.
    #[inline(always)]
    pub fn peek<const N: usize>(&mut self, data: &[u8]) -> Option<[usize; N]> {
        if self.len - self.head < N && !self.refill::<N>(data) {
            return None;
        }
        self.line_feeds[self.head..self.head + N].try_into().ok()
    }

    #[inline(always)]
    pub fn consume<const N: usize>(&mut self, next_start: usize) {
        self.head += N;
        self.next_start = next_start;
        self.consumed = true;
    }

    fn refill<const N: usize>(&mut self, data: &[u8]) -> bool {
        self.line_feeds.copy_within(self.head..self.len, 0);
        let mut len = self.len - self.head;
        self.head = 0;
        let end = scan_end(data).min(self.next_start + MAX_SCAN);
        let mut pos = self.scan_pos;
        while len < N && pos < end {
            let stop = (pos + 64 * WINDOW_BLOCKS).min(end);
            while pos < stop {
                len = push_line_feeds(&mut self.line_feeds, len, pos, scan_block(data, pos));
                pos += 64;
            }
        }
        self.scan_pos = pos;
        self.len = len;
        len >= N
    }

    // Called after the block parser emitted a record from `record_start` to `next_start`.
    #[inline(always)]
    pub fn record_parsed(&mut self, data: &[u8], record_start: usize, next_start: usize) {
        if next_start - record_start > MAX_RECORD_LEN {
            self.short_records = 0;
        } else {
            self.short_records += 1;
            if self.short_records >= MIN_SHORT_RECORDS << self.backoff {
                self.activate(data, next_start);
            }
        }
    }

    // Out of line to keep the block parser small.
    #[inline(never)]
    fn activate(&mut self, data: &[u8], next_start: usize) {
        let block_start = next_start & !63;
        if block_start >= scan_end(data) {
            return;
        }
        self.active = true;
        self.consumed = false;
        self.next_start = next_start;
        self.line_feeds.resize(64 * WINDOW_BLOCKS + 64 + 3, 0);
        let line_feeds = scan_block(data, block_start) & (!0 << (next_start % 64));
        self.head = 0;
        self.len = push_line_feeds(&mut self.line_feeds, 0, block_start, line_feeds);
        self.scan_pos = block_start + 64;
    }

    // Returns where the block parser resumes.
    #[inline(always)]
    pub fn deactivate(&mut self) -> usize {
        self.active = false;
        self.short_records = 0;
        self.backoff = if self.consumed {
            0
        } else {
            (self.backoff + 1).min(16)
        };
        self.next_start
    }
}

// The last block is never indexed: the block parser handles the end of input.
#[inline(always)]
fn scan_end(data: &[u8]) -> usize {
    data.len().saturating_sub(1) & !63
}

#[inline(always)]
pub(super) fn scan_block(data: &[u8], pos: usize) -> u64 {
    extract_bitmask::<LINE_FEEDS_ONLY, SCAN_MODE_HEADER>(&data[pos..pos + 64]).line_feeds
}

#[inline(always)]
fn push_line_feeds(out: &mut [usize], mut len: usize, pos: usize, mut line_feeds: u64) -> usize {
    while line_feeds != 0 {
        out[len] = pos + line_feeds.trailing_zeros() as usize;
        line_feeds &= line_feeds - 1;
        len += 1;
    }
    len
}
