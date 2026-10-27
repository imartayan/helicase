use super::*;
use fastx_gen::FastxGen;
use std::io::Cursor;

const CONFIG: Config = crate::ParserOptions::default().config();

fn count(mut parser: FastxParser<'_, CONFIG>) -> Result<(usize, usize)> {
    let (mut records, mut len) = (0, 0);
    while let Some(event) = parser.next() {
        event?;
        records += 1;
        len += parser.get_dna_string().len();
    }
    Ok((records, len))
}

fn count_stream(input: &[u8], threads: usize, chunk_len: usize) -> Result<(usize, usize)> {
    run_stream(Cursor::new(input), threads, chunk_len, count, |a, b| {
        (a.0 + b.0, a.1 + b.1)
    })
}

#[test]
fn stream_chunks_match_sequential() {
    for gen_records in [FastxGen::fasta, FastxGen::fastq] {
        let records = |n, len| gen_records(n, len).with_seed(len as u64).generate_bytes();
        // Records much longer than the chunks, at the start, in the middle and at the end
        let lens = [
            (1, 20_000),
            (300, 50),
            (1, 50_000),
            (2, 30_000),
            (100, 10),
            (1, 40_000),
        ];
        let input: Vec<u8> = lens.iter().flat_map(|&(n, len)| records(n, len)).collect();
        let expected = count(FastxParser::from_slice(&input).unwrap()).unwrap();
        for threads in [1, 3] {
            for chunk_len in [100, 4096] {
                let res = count_stream(&input, threads, chunk_len).unwrap();
                assert_eq!(res, expected, "{threads} threads, chunks of {chunk_len}");
            }
        }
    }
    assert_eq!(count_stream(b"", 2, 100).unwrap(), (0, 0));
}

#[test]
fn stream_stops_on_error() {
    let input = FastxGen::fasta(1000, 100).generate_bytes();
    let res =
        run_stream::<CONFIG, (), _, _>(Cursor::new(&input), 2, 100, |_| Err("".into()), |_, _| ());
    assert!(res.is_err());
    assert!(count_stream(b"xyz\n".repeat(100).as_slice(), 2, 100).is_err());
}
