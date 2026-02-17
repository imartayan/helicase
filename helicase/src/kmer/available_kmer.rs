//use crate::annotation::ConditionD;
use crate::kmer::*;
use crate::kmer::{gen_kmer};
use eyre::{eyre, Report, Result};
use paste::paste;
use std::any::Any;
use std::fmt;

// Record available macro to dispatch and panic if not available.
// Generate a dispatch_k macro that select the appropriate data store.
macro_rules! valid_k_t_pairs {
    ($([($($k:expr),*) $t:ty])*) => {
        paste!{
            gen_kmer!{$([($($k),*) $t])* }
            gen_sealed_shard!{ $([($($k),*) $t])* }
            gen_loose_build_shard!{ $([($($k),*) $t])* }
            gen_loose_build_view!{ $([($($k),*) $t])* }

        }
    };
}


// Define the valid mapping from k-mer lengths to types
valid_k_t_pairs! {
    [(4, 5, 6, 7, 8) u8]
    [(12, 13, 14, 15, 16) u16]
    [(21, 23, 25, 27, 29, 30, 31, 32) u32]
    [(53, 55, 57, 59, 61, 62, 63, 64) u64]
    [(101, 107, 113, 117, 123, 124, 125, 126, 127, 128) u128]
}
