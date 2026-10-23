// Example from needletail's README, using helicase as a drop-in replacement.
// Only the import of `parse_fastx_file` and the `normalize` call differ.
use helicase::parse_fastx_file;
use needletail::Sequence;

fn main() {
    let filename = std::env::args().nth(1).expect("No input file given");

    let mut n_bases = 0;
    let mut n_valid_kmers = 0;
    let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        // keep track of the total number of bases
        n_bases += seqrec.num_bases();
        // normalize to make sure all the bases are consistently capitalized,
        // newlines are already removed by helicase
        let norm_seq = seqrec.raw_seq().normalize(false);
        // we make a reverse complemented copy of the sequence first for
        // `canonical_kmers` to draw the complemented sequences from.
        let rc = norm_seq.reverse_complement();
        // now we keep track of the number of AAAAs (or TTTTs via
        // canonicalization) in the file
        for (_, kmer, _) in norm_seq.canonical_kmers(4, &rc) {
            if kmer == b"AAAA" {
                n_valid_kmers += 1;
            }
        }
    }
    println!("There are {n_bases} bases in your file.");
    println!("There are {n_valid_kmers} AAAAs in your file.");
}
