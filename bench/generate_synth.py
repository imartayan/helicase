
#!/usr/bin/env python3

import argparse
import random
import sys

NUCLEOTIDES = ("A", "C", "G", "T")
NON_NUCLEOTIDE = "N"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic FASTA files with controlled structure"
    )

    parser.add_argument(
        "seed",
        type=int,
        help="Random seed (deterministic output)",
    )

    parser.add_argument(
        "--records",
        type=int,
        required=True,
        help="Number of FASTA records",
    )

    parser.add_argument(
        "--lines-per-record",
        type=int,
        required=True,
        help="Number of sequence lines per record",
    )

    parser.add_argument(
        "--line-length",
        type=int,
        required=True,
        help="Length of each sequence line",
    )

    parser.add_argument(
        "--non-nucleotid",
        type=float,
        default=0.0,
        help="Fraction of non-nucleotid characters (N), between 0 and 1",
    )

    parser.add_argument(
        "--no-multiline",
        action="store_true",
        help="Emit one single sequence line per record",
    )

    return parser.parse_args()


def random_base(rng: random.Random, non_nucleotid_freq: float) -> str:
    if rng.random() < non_nucleotid_freq:
        return NON_NUCLEOTIDE
    return rng.choice(NUCLEOTIDES)


def generate_sequence(
    rng: random.Random,
    total_length: int,
    non_nucleotid_freq: float,
) -> str:
    return "".join(
        random_base(rng, non_nucleotid_freq)
        for _ in range(total_length)
    )


def main() -> None:
    args = parse_args()

    if not (0.0 <= args.non_nucleotid <= 1.0):
        sys.exit("--non-nucleotid must be between 0 and 1")

    rng = random.Random(args.seed)

    for i in range(1, args.records + 1):
        print(f">record_{i}")

        if args.no_multiline:
            total_len = args.lines_per_record * args.line_length
            seq = generate_sequence(
                rng,
                total_len,
                args.non_nucleotid,
            )
            print(seq)
        else:
            for _ in range(args.lines_per_record):
                seq = generate_sequence(
                    rng,
                    args.line_length,
                    args.non_nucleotid,
                )
                print(seq)


if __name__ == "__main__":
    main()
