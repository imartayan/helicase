use helicase::input::*;
use helicase::*;
use std::cmp::{max, min};
use std::path::PathBuf;

fn human_readable(n: f64) -> String {
    if n < 1_000.0 {
        format!("{:.2}", n)
    } else if n < 1_000_000.0 {
        format!("{:.2}K", n / 1_000.0)
    } else if n < 1_000_000_000.0 {
        format!("{:.2}M", n / 1_000_000.0)
    } else {
        format!("{:.2}G", n / 1_000_000_000.0)
    }
}

fn column_widths(table: &[Vec<String>]) -> Vec<usize> {
    let cols = table.first().map(|r| r.len()).unwrap_or(0);
    let mut widths = vec![0; cols];

    for row in table {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.len());
        }
    }

    widths
}
fn print_table(table: &[Vec<String>]) {
    let widths = column_widths(table);

    for row in table {
        for (i, cell) in row.iter().enumerate() {
            print!("{:<width$} ", cell, width = widths[i]);
        }
        println!();
    }
}

const MINIMAL: Config = ParserOptions::default().ignore_headers().ignore_dna().compute_dna_len().config();


fn main() {
    let mut table: Vec<Vec<String>> = vec![];
    let row = vec![
        "file".to_string(),
        "format".to_string(),
        "type".to_string(),
        "num_seqs".to_string(),
        "sum_len".to_string(),
        "min_len".to_string(),
        "avg_len".to_string(),
        "max_len".to_string(),
    ];
    table.push(row);
    for arg in std::env::args().skip(1) {
        let path = PathBuf::from(&arg);
        if path.exists() {
            let mut parser =
                FastxParser::<MINIMAL>::from_file_mmap(&path).expect("Cannot open file");
            let mut min_size = usize::MAX;
            let mut max_size = 0;
            let mut total_size = 0;
            let mut record_nb = 0;
            let format = parser.format();
            while parser.next().is_some() {
                record_nb += 1;
                let dna_len = parser.get_dna_len();
                total_size += dna_len;
                min_size = min(min_size, dna_len);
                max_size = max(max_size, dna_len);
            }
            let avg = (total_size as f64) / (record_nb as f64);
            let row = vec![
                format!("{}", path.display()),
                format!("{:?}", format),
                "?".to_string(),
                human_readable(record_nb as f64),
                human_readable(total_size as f64),
                human_readable(min_size as f64),
                human_readable(avg),
                human_readable(max_size as f64),
            ];
            table.push(row);
        } else {
            let row = vec![
                format!("{}", path.display()),
                format!("{}", "ERR" ),
                "".to_string(),
                "".to_string(),
                "".to_string(),
                "".to_string(),
                "".to_string(),
                "".to_string(),
            ];
            table.push(row);
        }
    }
    print_table(&table);
}
