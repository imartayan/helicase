
#!/bin/sh
set -eu

./datasets.sh
DATASETS_DIR="datasets"
RESULTS_DIR="bench_results"

rm -rf $RESULTS_DIR
mkdir -p "$RESULTS_DIR"

# Iterate over all files in datasets, skip checksum files
for dataset in "$DATASETS_DIR"/*; do
    case "$dataset" in
        *.sha256|*.sha1|*.gz) continue ;;  # skip checksum files
    esac

    # Extract filename
    filename=$(basename "$dataset")
    out="$RESULTS_DIR/${filename}.bench"

    echo "[BENCH] $filename -> $out"

    # Run the benchmark and save output
    RUSTFLAGS="-C target-cpu=native -Awarning" \
    cargo r -r --quiet --bin bench -- "$dataset" > "$out" 2> "$out.err"
done
