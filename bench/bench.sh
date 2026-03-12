
#!/bin/sh
set -eu

./datasets.sh
DATASETS_DIR="datasets"
RESULT="benchs.csv"

flags=""

# Iterate over all files in datasets, skip checksum files
for dataset in "$DATASETS_DIR"/*; do
    case "$dataset" in
        *.sha256|*.sha1|*.gz) continue ;;  # skip checksum files
    esac

    # Extract filename
    filename=$(basename "$dataset")

    echo "[BENCH] $filename -> $RESULT"

    # Run the benchmark and save output
    RUSTFLAGS="-C target-cpu=native" \
    cargo r -r --quiet --bin bench -- "$dataset" -c $flags >> "$RESULT"
    flags="-H"
done

for dataset in "$DATASETS_DIR"/*; do
    case "$dataset" in
        *.sha256|*.sha1|*.gz) continue ;;  # skip checksum files
    esac

    # Extract filename
    filename=$(basename "$dataset")

    echo "[BENCH] $filename -> $RESULT"

    # Run the benchmark and save output
    RUSTFLAGS="-C target-cpu=native" \
    cargo r -r --quiet --features helicase/no-pdep --bin bench -- "$dataset" -c $flags >> "$RESULT"
done
