use paste::paste;

#[macro_export]
macro_rules! gen_dispatch {
    ($([($($k:expr),*) $t:ty])*) => {
        #[macro_export]
        macro_rules! dispatch_k {
            ($kval:expr, |$K:ident, $T:ident| $body:expr) => {
                match $kval {
                    $(
                        $(
                            $k => {
                                const $K: usize = $k;
                                type $T = $t;
                                $body
                            }
                        ),*
                    ),*,
                    _ => panic!("unsupported k value: {}", $kval),
                }
            };
        }
    };
}

// Record available macro to dispatch and panic if not available.
// Generate a dispatch_k macro that select the appropriate data store.
macro_rules! valid_k_t_pairs {
    ($([($($k:expr),*) $t:ty])*) => {
        paste!{
            gen_dispatch! { $([($($k),*) $t])* }
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
