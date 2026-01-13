
#!/bin/sh
set -eu

DATASETS_DIR="datasets"
JOBS=3
CHECKSUM=1   # 1=enable checksum, 0=disable

# ---------------------------
# Options
# ---------------------------
for arg in "$@"; do
    case "$arg" in
        --no-checksum)
            CHECKSUM=0
            ;;
        --jobs=*)
            JOBS="${arg#*=}"
            ;;
        *)
            echo "Unknown option: $arg" >&2
            exit 1
            ;;
    esac
done

mkdir -p "$DATASETS_DIR"

# ---------------------------
# Real dataset downloads
# ---------------------------
URLS="
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/NIST_Illumina_2x250bps/D1_S1_L001_R2_007.fastq.gz
https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/PacBio_HiFi/15kb/m54328_180928_230446.Q20.fastq
"

printf "%s\n" $URLS | xargs -n 1 -P "$JOBS" sh -c '
set -eu
DATASETS_DIR="'"$DATASETS_DIR"'"
CHECKSUM="'"$CHECKSUM"'"

url="$1"
file="$(basename "$url")"
out="$DATASETS_DIR/$file"
checksum="$out.sha256"

# Download
if [ ! -f "$out" ]; then
    echo "[DOWNLOADING] $file"
    curl -sS -L --fail --retry 3 -o "$out" "$url"
else
    echo "[FOUND] $file"
fi

# Checksum
if [ "$CHECKSUM" -eq 1 ]; then
    if [ -f "$checksum" ]; then
        echo "[VERIFY] $file"
        sha256sum -c "$checksum" || {
            echo "[REDOWNLOAD] $file due to checksum mismatch"
            rm -f "$out"
            curl -sS -L --fail --retry 3 -o "$out" "$url"
            sha256sum "$out" > "$checksum"
        }
    else
        echo "[CHECKSUM CREATE] $file"
        sha256sum "$out" > "$checksum"
    fi
fi

# Unpack
case "$file" in
    *.gz)
        unpacked="$DATASETS_DIR/${file%.gz}"
        if [ ! -f "$unpacked" ]; then
            echo "[UNPACKING] $file"
            gunzip -k "$out"
        else
            echo "[UNPACKED] ${file%.gz}"
        fi
        ;;
esac
' _

# ---------------------------
# Synthetic FASTA generation
# ---------------------------
SYNTH_SCRIPT="generate_synth.py"
RECORDS=10000
LINES_PER_RECORD=100
SEED=42
POWERS="1 2 4 8 16 32 64 128 256"

printf "%s\n" $POWERS | xargs -n 1 -P "$JOBS" sh -c '
set -eu
DATASETS_DIR="'"$DATASETS_DIR"'"
SYNTH_SCRIPT="'"$SYNTH_SCRIPT"'"
RECORDS="'"$RECORDS"'"
LINES_PER_RECORD="'"$LINES_PER_RECORD"'"
SEED="'"$SEED"'"

pow="$1"
out="$DATASETS_DIR/synt${pow}.fa"
sha256file="$out.sha256"

generate() {
    python3 "$SYNTH_SCRIPT" "$SEED" \
        --records "$RECORDS" \
        --lines-per-record "$LINES_PER_RECORD" \
        --line-length "$pow" \
        > "$out"
}

# Verify existing SHA256
if [ -f "$out" ] && [ -f "$sha256file" ]; then
    actual_sha256=$(sha256sum "$out" | cut -d " " -f1)
    expected_sha256=$(cut -d " " -f1 "$sha256file")
    if [ "$actual_sha256" = "$expected_sha256" ]; then
        echo "[SKIP] $out already exists and SHA256 matches"
        exit 0
    else
        echo "[REGENERATE] $out SHA256 mismatch"
        generate
        sha256sum "$out" > "$sha256file"
    fi
else
    echo "[GENERATE] $out"
    generate
    sha256sum "$out" > "$sha256file"
fi
' _
