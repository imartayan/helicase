// Items with the same implementation in `FastaParser` and `FastqParser`.

// Getters of `HelicaseParser` that only depend on the shared fields.
macro_rules! shared_getters {
    () => {
        #[inline(always)]
        fn get_header(&self) -> &[u8] {
            if flag_is_not_set(CONFIG, COMPUTE_HEADER) {
                panic!("Parser config error: headers are ignored")
            }
            self.lexer.input.slice(self.header_range.clone())
        }

        #[inline(always)]
        fn get_header_owned(&mut self) -> Vec<u8> {
            self.get_header().to_vec()
        }

        #[inline(always)]
        fn get_dna_columnar(&self) -> &ColumnarDNA {
            if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                panic!("Parser config error: dna_columnar is not enabled")
            }
            &self.bitpackers.columnar
        }

        #[inline(always)]
        fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
            if flag_is_not_set(CONFIG, COMPUTE_DNA_COLUMNAR) {
                panic!("Parser config error: dna_columnar is not enabled")
            }
            self.bitpackers.take_columnar()
        }

        #[inline(always)]
        fn get_dna_packed(&self) -> &PackedDNA {
            if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
                panic!("Parser config error: dna_packed is not enabled")
            }
            &self.bitpackers.packed
        }

        #[inline(always)]
        fn get_dna_packed_owned(&mut self) -> PackedDNA {
            if flag_is_not_set(CONFIG, COMPUTE_DNA_PACKED) {
                panic!("Parser config error: dna_packed is not enabled")
            }
            self.bitpackers.take_packed()
        }

        #[inline(always)]
        fn get_mask_non_actg(&self) -> &BitMask {
            if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                panic!("Parser config error: mask_non_actg is not enabled")
            }
            &self.bitpackers.mask_non_actg
        }

        #[inline(always)]
        fn get_mask_non_actg_owned(&mut self) -> BitMask {
            if flag_is_not_set(CONFIG, COMPUTE_MASK_NON_ACTG) {
                panic!("Parser config error: mask_non_actg is not enabled")
            }
            self.bitpackers.take_mask_non_actg()
        }

        #[inline(always)]
        fn get_mask_n(&self) -> &BitMask {
            if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
                panic!("Parser config error: mask_n is not enabled")
            }
            &self.bitpackers.mask_n
        }

        #[inline(always)]
        fn get_mask_n_owned(&mut self) -> BitMask {
            if flag_is_not_set(CONFIG, COMPUTE_MASK_N) {
                panic!("Parser config error: mask_n is not enabled")
            }
            self.bitpackers.take_mask_n()
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
            self.pending_error.take().map(Error::from).map(Err)
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
