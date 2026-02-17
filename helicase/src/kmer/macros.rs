#[macro_export]
macro_rules! gen_kmer{
    ($([($($k:expr),*) $t:ty])*) => {
        paste!{
            #[derive(Clone, PartialEq, Debug)]
            pub enum Kmer {
                $(
                    $(
                        [<K $k>](Mer::<$k,$t>),
                    )*
                )*
            }
            $(
                $(
                    impl From<Mer<$k, $t>> for Kmer {
                        fn from(kmer: Mer<$k, $t>) -> Self {
                            Kmer::[<K $k>](kmer)
                        }
                    }
                )*
            )*

            impl<'a> TryFrom<&'a str> for Kmer {
                type Error = Report;
                fn try_from(data: &'a str) -> Result<Self> {
                    match data.len(){
                        $(
                            $(
                                $k => {
                                    let t: Mer<$k,$t> = data.as_bytes().try_into()?;
                                    Ok(Self::[<K $k>](t))
                                },
                            )*
                        )*
                        x => Err(eyre!("Unsupported k-mer length, {}", x)),
                    }
                }
            }


            impl Kmer {
                pub fn k(&self) -> usize  {
                    match &self {
                        $(
                            $(
                                Self::[<K $k>](_) => $k,
                            )*
                        )*
                    }
                }

                pub fn to_string(&self) -> String {
                    match &self {
                        $(
                            $(
                                Self::[<K $k>](mer) => mer.to_string(),
                            )*
                        )*
                    }
                }
            }

            impl fmt::Display for Kmer {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, "{}", self.to_string())
                }
            }


            pub enum KmerChunk {
                $(
                    $(
                        [<K $k>](MerChunk<$k, $t>),
                    )*
                )*
            }
            impl KmerChunk {
                pub fn new(k:usize) -> Self {
                    match k {
                        $(
                            $(
                                $k => Self::[<K $k>](MerChunk::<$k,$t>::new()),
                            )*
                        )*

                        _ => panic!("Unsupported k-mer length: {}", k),

                    }
                }
                pub fn append(&mut self, other: &mut Self) {
                    match (self, other) {
                        $(
                            $(
                                (Self::[<K $k>](ref mut me), Self::[<K $k>](ref mut them)) => { me.append(them); },
                            )*
                        )*
                            _=> panic!("Append KmerCunk of different size"),
                    }
                }

                pub fn k(&self) -> usize  {
                    match &self {
                        $(
                            $(
                                Self::[<K $k>](_) => $k,
                            )*
                        )*
                    }
                }
                pub fn len(&self) -> usize  {
                    match &self {
                        $(
                            $(
                                Self::[<K $k>](chunk) => chunk.len(),
                            )*
                        )*
                    }
                }
                pub fn to_vec(&self) -> Vec<Kmer> {
                    match &self {
                        $(
                            $(
                                Self::[<K $k>](chunk) => {
                                    let it = chunk.iter();
                                    it.map(|mer: Mer<$k,$t> | Kmer::[<K $k>](mer)).collect()
                                }
                            )*
                        )*
                    }
                }
                pub fn from_dna(dna: DNA, k: usize) -> Self {
                    match k {
                        $(
                            $(
                                $k => {
                                    let chunk : MerChunk<$k, $t> = dna.into();
                                    KmerChunk::[<K $k>](chunk)
                                },
                            )*
                        )*

                        _ => panic!("Unsupported k-mer length: {}", k),
                    }
                }
            }
            $(
                $(
                    impl From<MerChunk<$k, $t>> for KmerChunk {
                        fn from(chunk: MerChunk<$k, $t>) -> Self {
                            KmerChunk::[<K $k>](chunk)
                        }
                    }
                )*
            )*

            impl FromIterator<KmerChunk> for KmerChunk{
                fn from_iter<I>(pre_iter: I) -> Self
                    where I: IntoIterator<Item = KmerChunk>{
                        let mut iter = pre_iter.into_iter();
                        let mut first = iter.next().unwrap(); // error if empty iterator :(
                        for mut chunk in iter {
                           first.append(&mut chunk);
                        }
                        first
                }
            }
        }
    }
}
