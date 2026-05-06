#!/usr/bin/env bash
# benchmark_pipeline_clink.sh — Clink-only pipeline speed + RAM benchmark
#
# Usage:
#   ./benchmark_pipeline_clink.sh -i <input.bam> -o <output_dir> [options]
#
# Options:
#   -i   Input BAM (sorted, indexed)
#   -o   Output directory
#   -n   Number of reads to extract (default: 1,000,000)
#   -t   Threads for pileup.py chromosome-level parallelism (default: 1)
#   -u   UMI length for Clink collapse (default: auto-detect)
#   --fdr  FDR threshold for CIMS/CITS (default: 0.05)

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
N_READS=1000000
INPUT_BAM=""
OUTDIR=""
THREADS=1
UMI_LEN="-1"
FDR=0.05
MIN_SIGNAL=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLINK_DIR="${SCRIPT_DIR}/lib/clink"
export PYTHONPATH="${CLINK_DIR}:${PYTHONPATH:-}"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    echo ""
    echo "Usage: $0 -i <input.bam> -o <output_dir> [options]"
    echo ""
    echo "  -i            Sorted, indexed input BAM"
    echo "  -o            Output directory"
    echo "  -n            Reads to extract (default: 1,000,000)"
    echo "  -t            Threads for pileup.py (default: 1)"
    echo "  -u            UMI length (default: auto-detect from read names)"
    echo "  --fdr         FDR/p-value threshold (default: 0.05)"
    echo "  --min-signal  Minimum mutation/truncation read count per site (default: 1)"
    echo ""
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) INPUT_BAM="$2";        shift 2 ;;
        -o) OUTDIR="$2";           shift 2 ;;
        -n) N_READS="$2";          shift 2 ;;
        -t) THREADS="$2";          shift 2 ;;
        -u) UMI_LEN="$2";          shift 2 ;;
        --fdr) FDR="$2";           shift 2 ;;
        --min-signal) MIN_SIGNAL="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
done

[[ -z "$INPUT_BAM" ]] && { echo "[ERROR] -i <input.bam> required"; usage; }
[[ -z "$OUTDIR"    ]] && { echo "[ERROR] -o <output_dir> required"; usage; }
[[ ! -f "$INPUT_BAM" ]] && { echo "[ERROR] BAM not found: $INPUT_BAM"; exit 1; }
command -v gtime    &>/dev/null || { echo "[ERROR] gtime not found (brew install gnu-time)"; exit 1; }
command -v samtools &>/dev/null || { echo "[ERROR] samtools not in PATH"; exit 1; }
[[ -f "${CLINK_DIR}/pileup.py"   ]] || { echo "[ERROR] pileup.py not found: ${CLINK_DIR}";   exit 1; }
[[ -f "${CLINK_DIR}/collapse.py" ]] || { echo "[ERROR] collapse.py not found: ${CLINK_DIR}"; exit 1; }
[[ -f "${CLINK_DIR}/cits.py"     ]] || { echo "[ERROR] cits.py not found: ${CLINK_DIR}";     exit 1; }
[[ -f "${CLINK_DIR}/cims.py"     ]] || { echo "[ERROR] cims.py not found: ${CLINK_DIR}";     exit 1; }

CLK_DIR="${OUTDIR}/Clink"
mkdir -p "$CLK_DIR"

# Remove stale Clink output so old files don't pollute counts if a step fails
rm -f "${CLK_DIR}/sample"_*.bed

REPORT="${OUTDIR}/benchmark_clink_report.txt"

# ── Helpers ───────────────────────────────────────────────────────────────────

fmt_num() { python3 -c "print(f'{int($1):,}')"; }

count_data() {
    local f="$1"
    [[ ! -s "$f" ]] && echo 0 && return
    local n; n=$(wc -l < "$f")
    echo $(( n > 0 ? n - 1 : 0 ))
}

bam_count() {
    samtools view -c -F 4 "$1" 2>/dev/null || echo 0
}

run_timed() {
    local label="$1"
    local gt_file="$2"
    local log_file="${gt_file%.txt}.log"
    shift 2
    printf "    %-35s " "$label"

    gtime -v "$@" >"$log_file" 2>&1 || true

    grep -E "Elapsed \(wall clock\)|Maximum resident|User time|System time" \
        "$log_file" > "$gt_file" 2>/dev/null || true

    local raw_wall
    raw_wall=$(grep "Elapsed (wall clock)" "$gt_file" | awk '{print $NF}')
    local wall_s
    wall_s=$(echo "$raw_wall" | awk -F: '{
        if (NF==2) printf "%.2f", $1*60+$2
        else       printf "%.2f", $1*3600+$2*60+$3
    }')
    local rss_kb
    rss_kb=$(grep "Maximum resident" "$gt_file" | awk '{print $NF}')
    local rss_mb
    rss_mb=$(awk "BEGIN {printf \"%.1f\", $rss_kb/1024}")
    printf "%8ss  %8s MB\n" "$wall_s" "$rss_mb"
    LAST_WALL="$wall_s"
    LAST_RSS="$rss_mb"
}

accumulate() {
    TOTAL_WALL=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_WALL:-0} + $1}")
    local rss="$2"
    PEAK_RSS=$(awk "BEGIN { a=${PEAK_RSS:-0}; b=${rss}; print (a>b?a:b) }")
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  CLIPittyClip — Clink Pipeline Benchmark"
echo "════════════════════════════════════════════════════════════════"
echo "  Input BAM   : $INPUT_BAM"
echo "  Reads       : $(fmt_num $N_READS)"
echo "  Threads     : $THREADS  (pileup.py)"
echo "  FDR         : $FDR"
echo "  Output dir  : $OUTDIR"
echo "════════════════════════════════════════════════════════════════"
echo ""

# ── Step 0: Subset BAM ────────────────────────────────────────────────────────
echo "[0] Extracting $(fmt_num $N_READS) reads from BAM..."

SUBSET_BAM="${OUTDIR}/subset_${N_READS}.bam"
SUBSET_SORTED="${OUTDIR}/subset_sorted.bam"

set +o pipefail
{ samtools view -H "$INPUT_BAM"
  samtools view -F 4 "$INPUT_BAM" | head -n "$N_READS"
} | samtools view -bS -o "${SUBSET_BAM}.tmp"
set -o pipefail

samtools sort -o "$SUBSET_SORTED" "${SUBSET_BAM}.tmp"
samtools index "$SUBSET_SORTED"
rm -f "${SUBSET_BAM}.tmp"

INPUT_READS=$(bam_count "$SUBSET_SORTED")
echo "    Extracted: $(fmt_num $INPUT_READS) mapped reads"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
echo "── Clink Pipeline ───────────────────────────────────────────────────────"
echo "    Step                                Wall       Peak RAM"
echo "    ─────────────────────────────────────────────────────────"

TOTAL_WALL=0; PEAK_RSS=0

# ── Clink 1: collapse.py (umi_tools dedup) ───────────────────────────────────
DEDUP_BAM="${CLK_DIR}/dedup.bam"
DEDUP_LOG="${CLK_DIR}/dedup.log"

CLINK_COLLAPSE_ARGS="--bam '$SUBSET_SORTED' --out '$DEDUP_BAM'"
[[ "$UMI_LEN" != "-1" ]] && CLINK_COLLAPSE_ARGS="$CLINK_COLLAPSE_ARGS --umi-len $UMI_LEN"
CLINK_COLLAPSE_ARGS="$CLINK_COLLAPSE_ARGS --log '$DEDUP_LOG'"

run_timed "collapse.py (umi_tools dedup)" "${CLK_DIR}/gt_dedup.txt" \
    bash -c "python3 '${CLINK_DIR}/collapse.py' $CLINK_COLLAPSE_ARGS"
accumulate "$LAST_WALL" "$LAST_RSS"

CLK_DEDUP=$(bam_count "$DEDUP_BAM")

# ── Clink 2: pileup.py ───────────────────────────────────────────────────────
NPZ="${CLK_DIR}/pileup.npz"
run_timed "pileup.py (t=${THREADS})" "${CLK_DIR}/gt_pileup.txt" \
    python3 "${CLINK_DIR}/pileup.py" "$DEDUP_BAM" --out "$NPZ" --threads "$THREADS"
accumulate "$LAST_WALL" "$LAST_RSS"

# ── Clink 3: cits.py ─────────────────────────────────────────────────────────
CLK_PREFIX="${CLK_DIR}/sample"
run_timed "cits.py" "${CLK_DIR}/gt_cits.txt" \
    python3 "${CLINK_DIR}/cits.py" \
        --pileup "$NPZ" --prefix "$CLK_PREFIX" \
        --min-cov 5 --min-frac 0.05 --min-signal "$MIN_SIGNAL" --fdr "$FDR"
accumulate "$LAST_WALL" "$LAST_RSS"

CLK_TRUNC=$(count_data "${CLK_PREFIX}_truncations.bed")

# ── Clink 4: cims.py ─────────────────────────────────────────────────────────
run_timed "cims.py" "${CLK_DIR}/gt_cims.txt" \
    python3 "${CLINK_DIR}/cims.py" \
        --pileup "$NPZ" --prefix "$CLK_PREFIX" \
        --min-cov 5 --min-frac 0.05 --min-signal "$MIN_SIGNAL" --fdr "$FDR"
accumulate "$LAST_WALL" "$LAST_RSS"

CLK_DELS=$(count_data "${CLK_PREFIX}_deletions.bed")

CLK_SUBS_TOTAL=0
for sub_bed in "${CLK_PREFIX}"_*to*.bed; do
    [[ -f "$sub_bed" ]] || continue
    n=$(count_data "$sub_bed")
    CLK_SUBS_TOTAL=$(( CLK_SUBS_TOTAL + n ))
done
CLK_TtoC=$(count_data "${CLK_PREFIX}_TtoC.bed")

echo "    ─────────────────────────────────────────────────────────"
printf "    %-35s %8ss  %8s MB\n" "TOTAL" "$TOTAL_WALL" "$PEAK_RSS"
echo ""

# ── Report ────────────────────────────────────────────────────────────────────
{
printf "\n════════════════════════════════════════════════════════════════\n"
printf "  Clink Benchmark Report\n"
printf "════════════════════════════════════════════════════════════════\n"
printf "  Input BAM     : %s\n" "$INPUT_BAM"
printf "  Reads tested  : %s\n" "$(fmt_num $INPUT_READS)"
printf "  Threads       : %s  (pileup.py)\n" "$THREADS"
printf "  FDR           : %s\n" "$FDR"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s\n" "Total wall time (s)"      "$TOTAL_WALL"
printf "  %-38s  %10s\n" "Peak RAM (MB)"             "$PEAK_RSS"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s\n" "OUTPUT COUNTS"             ""
printf "  %-38s  %10s\n" "──────────────────────"    "──────────"
printf "  %-38s  %10s\n" "Input reads"               "$(fmt_num $INPUT_READS)"
printf "  %-38s  %10s\n" "Dedup reads"               "$(fmt_num $CLK_DEDUP)"
printf "  %-38s  %10s\n" "CIMS deletion sites"       "$(fmt_num $CLK_DELS)"
printf "  %-38s  %10s\n" "CIMS substitution sites"   "$(fmt_num $CLK_SUBS_TOTAL)"
printf "  %-38s  %10s\n" "  of which T>C"            "$(fmt_num $CLK_TtoC)"
printf "  %-38s  %10s\n" "CITS truncation sites"     "$(fmt_num $CLK_TRUNC)"
printf "════════════════════════════════════════════════════════════════\n"
} | tee "$REPORT"

echo ""
echo "  Report saved  : $REPORT"
echo "  Clink outputs : $CLK_DIR/"
echo ""
