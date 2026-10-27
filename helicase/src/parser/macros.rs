// Items with the same implementation in `FastaParser` and `FastqParser`.

// Panics if `$flag` is not set: getters of `FastxParser` check it at compile time instead.
macro_rules! check {
    ($flag:expr, $msg:literal) => {
        if flag_is_not_set(CONFIG, $flag) {
            panic!(concat!("Parser config error: ", $msg))
        }
    };
}

// Getters of `HelicaseParser` that only depend on the shared fields.
macro_rules! shared_getters {
    () => {
        #[inline(always)]
        fn get_header(&self) -> &[u8] {
            check!(COMPUTE_HEADER, "headers are ignored");
            self.lexer.input.slice(self.header_range.clone())
        }

        #[inline(always)]
        fn get_header_owned(&mut self) -> Vec<u8> {
            self.get_header().to_vec()
        }

        #[inline(always)]
        fn get_dna_columnar(&self) -> &ColumnarDNA {
            check!(COMPUTE_DNA_COLUMNAR, "dna_columnar is not enabled");
            &self.bitpackers.columnar
        }

        #[inline(always)]
        fn swap_dna_columnar(&mut self, other: &mut ColumnarDNA) {
            check!(COMPUTE_DNA_COLUMNAR, "dna_columnar is not enabled");
            core::mem::swap(&mut self.bitpackers.columnar, other);
        }

        #[inline(always)]
        fn get_dna_packed(&self) -> &PackedDNA {
            check!(COMPUTE_DNA_PACKED, "dna_packed is not enabled");
            &self.bitpackers.packed
        }

        #[inline(always)]
        fn swap_dna_packed(&mut self, other: &mut PackedDNA) {
            check!(COMPUTE_DNA_PACKED, "dna_packed is not enabled");
            core::mem::swap(&mut self.bitpackers.packed, other);
            self.bitpackers.pad::<CONFIG>();
        }

        #[inline(always)]
        fn get_mask_non_actg(&self) -> &BitMask {
            check!(COMPUTE_MASK_NON_ACTG, "mask_non_actg is not enabled");
            &self.bitpackers.mask_non_actg
        }

        #[inline(always)]
        fn swap_mask_non_actg(&mut self, other: &mut BitMask) {
            check!(COMPUTE_MASK_NON_ACTG, "mask_non_actg is not enabled");
            core::mem::swap(&mut self.bitpackers.mask_non_actg, other);
        }

        #[inline(always)]
        fn get_mask_n(&self) -> &BitMask {
            check!(COMPUTE_MASK_N, "mask_n is not enabled");
            &self.bitpackers.mask_n
        }

        #[inline(always)]
        fn swap_mask_n(&mut self, other: &mut BitMask) {
            check!(COMPUTE_MASK_N, "mask_n is not enabled");
            core::mem::swap(&mut self.bitpackers.mask_n, other);
        }

        #[inline(always)]
        fn get_dna_len(&self) -> usize {
            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                self.dna_len
            } else if flag_is_set(CONFIG, COMPUTE_DNA_STRING) {
                self.get_dna_string().len()
            } else if flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                self.bitpackers.columnar.len()
            } else if flag_is_set(CONFIG, COMPUTE_DNA_PACKED) {
                self.bitpackers.packed.len()
            } else {
                panic!("Parser config error: dna is ignored")
            }
        }

        #[cfg(feature = "paraseq")]
        #[inline(always)]
        fn record_index(&self) -> u64 {
            self.record_index
        }
    };
}

// Block-level helpers of the parsers.
macro_rules! shared_block_helpers {
    () => {
        #[inline(always)]
        fn return_error_or_none(&mut self) -> Option<Result<Event>> {
            self.pending_error.take().map(Err)
        }

        #[inline(always)]
        const fn global_pos(&self) -> usize {
            64 * self.block_counter + self.pos_in_block
        }

        #[inline(always)]
        fn next_block<const MODE: u8>(&mut self) -> bool {
            match self.lexer.next_block::<MODE>() {
                Some(Ok(b)) => {
                    self.block.set::<CONFIG>(b);
                    self.block_counter += 1;
                    true
                }
                Some(Err(e)) => {
                    self.pending_error = Some(e);
                    self.finished = true;
                    false
                }
                None => {
                    self.finished = true;
                    false
                }
            }
        }

        #[inline(always)]
        fn current_byte(&self) -> u8 {
            if self.pos_in_block < self.block.len {
                self.lexer.input.current_block()[self.pos_in_block & 63]
            } else {
                0
            }
        }

        #[inline(always)]
        fn increment_pos(&mut self) {
            if self.pos_in_block + 1 < self.block.len {
                self.pos_in_block += 1;
            } else if self.next_block::<SCAN_MODE_DNA>() {
                self.pos_in_block = 0;
            }
        }

        #[inline(always)]
        fn consume_newline(&mut self) {
            self.block.newline &= self.block.newline.wrapping_sub(1);
            self.increment_pos();
        }

        // Keeps the record from `anchor(line_start)` when refilling the buffer of a reader.
        #[inline(always)]
        fn make_room_record(&mut self, line_start: usize) {
            if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
                && !I::RANDOM_ACCESS
                && self.lexer.input.is_end_of_buffer()
            {
                self.lexer.input.make_room(self.anchor(line_start));
            }
        }

        // Same as `make_room_record`, before moving past the last byte of the block.
        #[inline(always)]
        fn make_room_record_end(&mut self, line_start: usize) {
            if flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING)
                && !I::RANDOM_ACCESS
                && self.pos_in_block + 1 >= self.block.len
                && self.lexer.input.is_end_of_buffer()
            {
                self.lexer.input.make_room(self.anchor(line_start));
            }
        }

        // Length of the `\r` before the line feed at the current position, none at the end of input.
        #[inline(always)]
        const fn crlf_before(&self) -> usize {
            if self.pos_in_block < self.block.len {
                self.crlf_len
            } else {
                0
            }
        }

        #[inline(always)]
        fn detect_crlf(&mut self, nl_pos: usize) {
            if !self.crlf_detected {
                if nl_pos > 0 && self.lexer.input.slice(nl_pos - 1..nl_pos) == b"\r" {
                    self.crlf_len = 1;
                }
                self.crlf_detected = true;
            }
        }

        // Appends the DNA up to the next line feed (or non-ACTG base when splitting) without the
        // `\r` of CRLF line endings. Returns `false` at the end of input.
        #[inline(always)]
        fn append_dna_run(&mut self, line_start: usize) -> bool {
            let mut position = !self.block.is_dna & (!0 << self.pos_in_block);
            if position != 0 {
                let first_pos = self.pos_in_block;
                self.pos_in_block = position.trailing_zeros() as usize;
                self.append_dna_range(first_pos, self.without_cr(self.pos_in_block - first_pos));
                return true;
            }
            if self.pos_in_block > 0 {
                self.make_room_record(line_start);
                self.append_dna_range(self.pos_in_block, 64 - self.pos_in_block);
                if !self.next_block::<SCAN_MODE_DNA>() {
                    return false;
                }
                self.pos_in_block = 0;
                position = !self.block.is_dna;
            }
            let mut dna_blocks: usize = 0;
            while position == 0 {
                self.make_room_record(line_start);
                self.append_dna_block();
                dna_blocks += 1;
                if dna_blocks == 16 && !Self::MERGE_DNA_STRING {
                    let anchor = self.anchor(line_start);
                    self.lexer.skip_long_dna(
                        &mut self.bitpackers,
                        &mut self.dna_len,
                        &mut self.block_counter,
                        anchor,
                        flag_is_set(CONFIG, COMPUTE_HEADER | COMPUTE_DNA_STRING),
                    );
                    dna_blocks = 15;
                }
                if !self.next_block::<SCAN_MODE_DNA>() {
                    return false;
                }
                position = !self.block.is_dna;
            }
            self.pos_in_block = position.trailing_zeros() as usize;
            if self.pos_in_block > 0 {
                self.append_dna_range(0, self.without_cr(self.pos_in_block));
            } else if self.crlf_len == 1 && flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                // The `\r` was the last byte of the previous block.
                self.bitpackers.pop_base::<CONFIG>();
                if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                    self.dna_len -= 1;
                }
            }
            true
        }

        #[inline(always)]
        const fn without_cr(&self, len: usize) -> usize {
            if flag_is_not_set(CONFIG, SPLIT_NON_ACTG) {
                len.saturating_sub(self.crlf_before())
            } else {
                len
            }
        }

        #[inline(always)]
        fn append_dna_range(&mut self, start: usize, len: usize) {
            self.bitpackers
                .append_range_from::<CONFIG>(&self.block, start, len);
            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                self.dna_len += len;
            }
        }

        #[inline(always)]
        fn append_dna_block(&mut self) {
            self.bitpackers.append_block_from::<CONFIG>(&self.block);
            if flag_is_set(CONFIG, COMPUTE_DNA_LEN) {
                self.dna_len += 64;
            }
        }

        // Ends the iteration: `next_event` returns `None` once finished with `REJECT_NON_ACTG`.
        #[cold]
        #[inline(never)]
        fn reject_base(&mut self, byte: u8, pos: usize) -> Option<Result<Event>> {
            self.finished = true;
            Some(Err(Error::invalid_base(byte, pos)))
        }

        #[cold]
        #[inline(never)]
        fn reject_first_non_actg(&mut self, range: Range<usize>) -> Option<Result<Event>> {
            let seq = self.lexer.input.slice(range.clone());
            let offset = seq
                .iter()
                .position(|c| !matches!(c & 0xDF, b'A' | b'C' | b'G' | b'T'))
                .unwrap_or(0);
            self.reject_base(seq[offset], range.start + offset)
        }

        #[cold]
        fn leave_indexed(&mut self) {
            let start = self.line_index.deactivate();
            self.block_counter = start / 64;
            self.pos_in_block = start % 64;
            self.block = self.lexer.seek_chunk(start);
        }
    };
}

// Iterator over the events of a parser, padding the packed DNA after each event.
macro_rules! shared_iterator {
    ($parser:ident) => {
        impl<'a, const CONFIG: Config, I: InputData<'a>> Iterator for $parser<'a, CONFIG, I> {
            type Item = Result<Event>;

            #[inline(always)]
            fn next(&mut self) -> Option<Self::Item> {
                let event = self.next_event();
                self.bitpackers.pad::<CONFIG>();
                event
            }
        }
    };
}
