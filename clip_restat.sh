#!/usr/bin/env bash
# clip_restat.sh — Clink and/or CTK crosslink site calling from a pre-dedup BAM
#
# Both tools start from the same input BAM (optionally subsampled to N reads).
# Prints a side-by-side timing and site-count summary when done.
#
# Clink path:  BAM → pileup.py (npz) → cits.py + cims.py        [BH FDR, minutes]
# CTK path:    BAM → calmd/SAM → parseAlignment.pl → CIMS.pl + CITS.pl  [permutation FDR, hours]
#
# Usage:
#   clip_restat.sh --bam <dedup.bam> --out <dir> --name <str> [options]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_DIR="$SCRIPT_DIR/lib"
CLINK_SCRIPTS="$LIB_DIR/clink"

# ── Defaults ──────────────────────────────────────────────────────────────────
BAM=""
OUT_DIR=""
NAME=""
GENOME=""
SUBSAMPLE=0          # 0 = no subsampling; >0 = target read count
FDR="0.05"
CITS_PVAL=""         # defaults to FDR if not set
CIMS_ITER=10
MIN_COV=5
MIN_FRAC=0.05
THREADS=4
NH=1                 # Clink NH tag threshold: 1=unique, 100=multi-mappers
MAPQ=255             # Clink MAPQ filter: 255=STAR unique, 1=multi-mappers
RUN_CLINK="true"
RUN_CTK="true"
RUN_CITS="true"
RUN_CIMS="true"
VERBOSE="false"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF

Usage: $(basename "$0") --bam <dedup.bam> --out <dir> --name <str> [options]

Required:
  --bam   <path>       Pre-deduplicated BAM file (sorted + indexed)
  --out   <dir>        Output root directory
  --name  <str>        Sample name prefix for output files

Subsampling:
  --subsample <int>    Subsample to N reads before running either tool
                       Both tools receive the exact same subsampled BAM

Thresholds (shared):
  --fdr   <float>      FDR threshold — Clink BH FDR + CTK CIMS FDR   (default: $FDR)
  --cits-pval <float>  CTK CITS p-value threshold  (default: same as --fdr)
  --cims-iter <int>    CTK permutation iterations                      (default: $CIMS_ITER)
  --min-cov   <int>    Minimum read coverage for Clink                 (default: $MIN_COV)
  --min-frac  <float>  Minimum mutation fraction for Clink             (default: $MIN_FRAC)

Alignment filters (Clink pileup):
  --threads <int>      Pileup threads                                  (default: $THREADS)
  --nh <int>           NH tag threshold: 1=unique, 100=multi-mappers   (default: $NH)
  --mapq <int>         MAPQ filter: 255=STAR unique, 1=multi-mappers   (default: $MAPQ)

Genome:
  --genome <path>      Reference FASTA — used by CTK for samtools calmd
                       (strongly recommended for CTK; skipped if not provided)

Tool selection:
  --no-clink           Skip Clink pipeline
  --no-ctk             Skip CTK pipeline
  --no-cits            Skip CITS (both tools)
  --no-cims            Skip CIMS (both tools)

Other:
  -v, --verbose        Print tool output to terminal (also always logged)
  -h, --help           Show this help

Output layout:
  <out>/
    input/<name>[_subN].bam   subsampled BAM (if --subsample used)
    Clink/<name>_pileup.npz   pileup arrays
    Clink/<name>_truncations.bed
    Clink/<name>_deletions.bed + substitution beds
    CTK/preprocessing/        parseAlignment + mutation files
    CTK/CIMS/                 CIMS_del.txt + CIMS_sub.txt
    CTK/CITS/                 CITS.bed
    clip_restat.log
    summary.txt

EOF
    exit 0
}

# ── Parse args ────────────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)        BAM="$2";        shift 2 ;;
        --out)        OUT_DIR="$2";    shift 2 ;;
        --name)       NAME="$2";       shift 2 ;;
        --genome)     GENOME="$2";     shift 2 ;;
        --subsample)  SUBSAMPLE="$2";  shift 2 ;;
        --fdr)        FDR="$2";        shift 2 ;;
        --cits-pval)  CITS_PVAL="$2";  shift 2 ;;
        --cims-iter)  CIMS_ITER="$2";  shift 2 ;;
        --min-cov)    MIN_COV="$2";    shift 2 ;;
        --min-frac)   MIN_FRAC="$2";   shift 2 ;;
        --threads)    THREADS="$2";    shift 2 ;;
        --nh)         NH="$2";         shift 2 ;;
        --mapq)       MAPQ="$2";       shift 2 ;;
        --no-clink)   RUN_CLINK="false"; shift ;;
        --no-ctk)     RUN_CTK="false";   shift ;;
        --no-cits)    RUN_CITS="false";  shift ;;
        --no-cims)    RUN_CIMS="false";  shift ;;
        -v|--verbose) VERBOSE="true";    shift ;;
        -h|--help)    usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$BAM" ]]     && { echo "ERROR: --bam is required" >&2;  exit 1; }
[[ -z "$OUT_DIR" ]] && { echo "ERROR: --out is required" >&2;  exit 1; }
[[ -z "$NAME" ]]    && { echo "ERROR: --name is required" >&2; exit 1; }
[[ ! -f "$BAM" ]]   && { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }

[[ "$RUN_CLINK" == "false" && "$RUN_CTK" == "false" ]] && {
    echo "ERROR: --no-clink and --no-ctk together — nothing to run" >&2; exit 1
}

if [[ "$RUN_CTK" == "true" && -z "$GENOME" ]]; then
    echo "WARNING: --genome not provided. CTK calmd will be skipped." \
         "MD tag quality may be reduced." >&2
fi

# ── Setup ─────────────────────────────────────────────────────────────────────
[[ -z "$CITS_PVAL" ]] && CITS_PVAL="$FDR"

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/clip_restat.log"

# Globals consumed by modules.sh logging/exec functions
export LOG_FILE="$LOG"
export VERBOSE

# Source CLIPittyClip library (CTK functions, logging, utils)
source "$LIB_DIR/utils.sh"
source "$LIB_DIR/modules.sh"

# ── Helpers ───────────────────────────────────────────────────────────────────

# Format seconds → H:MM:SS
fmt_time() {
    local s=$1
    printf "%d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# Count significant lines in a BED/TXT output (subtract 1 for header)
count_sites() {
    local f="$1"
    [[ -f "$f" ]] || { echo "—"; return; }
    local n
    n=$(grep -c -v '^#' "$f" 2>/dev/null || echo 0)
    echo "$n"
}

# Log and optionally echo to terminal
log_step() {
    local msg="$1"
    echo "$msg" | tee -a "$LOG"
}

# Run a command, always logging; also echo stdout if VERBOSE
run_cmd() {
    if [[ "$VERBOSE" == "true" ]]; then
        "$@" 2>&1 | tee -a "$LOG"
        return "${PIPESTATUS[0]}"
    else
        "$@" >> "$LOG" 2>&1
    fi
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════" | tee "$LOG"
echo "  clip_restat.sh"                                         | tee -a "$LOG"
echo "  BAM:     $BAM"                                          | tee -a "$LOG"
echo "  Name:    $NAME"                                         | tee -a "$LOG"
echo "  Output:  $OUT_DIR"                                      | tee -a "$LOG"
[[ "$SUBSAMPLE" -gt 0 ]] && \
echo "  Subset:  $SUBSAMPLE reads"                              | tee -a "$LOG"
echo "  FDR:     $FDR  |  CITS p: $CITS_PVAL  |  iter: $CIMS_ITER (CTK)" | tee -a "$LOG"
echo "  Clink:   $RUN_CLINK  |  CTK: $RUN_CTK"                 | tee -a "$LOG"
echo "  CITS:    $RUN_CITS   |  CIMS: $RUN_CIMS"               | tee -a "$LOG"
echo "════════════════════════════════════════════════════════" | tee -a "$LOG"

# ── Subsample ─────────────────────────────────────────────────────────────────
INPUT_BAM="$BAM"
ACTUAL_READS=""

if [[ "$SUBSAMPLE" -gt 0 ]]; then
    mkdir -p "$OUT_DIR/input"
    SUB_BAM="$OUT_DIR/input/${NAME}_sub${SUBSAMPLE}.bam"

    log_step "[subsample] Counting reads in input BAM..."
    TOTAL=$(samtools view -c "$BAM")
    log_step "[subsample] Input: $TOTAL reads → target: $SUBSAMPLE"

    if [[ "$TOTAL" -le "$SUBSAMPLE" ]]; then
        log_step "[subsample] Target ≥ input; using full BAM (no subsampling)"
        ln -sf "$(realpath "$BAM")" "$SUB_BAM"
        [[ -f "${BAM}.bai" ]] && ln -sf "$(realpath "${BAM}.bai")" "${SUB_BAM}.bai"
        ACTUAL_READS="$TOTAL"
    else
        # samtools -s: integer = random seed, decimal = sampling fraction
        SVAL=$(python3 -c "print(f'{42 + $SUBSAMPLE/$TOTAL:.6f}')")
        log_step "[subsample] Sampling fraction: $SVAL (seed 42)"
        samtools view -s "$SVAL" -b "$BAM" > "$SUB_BAM"
        samtools index "$SUB_BAM"
        ACTUAL_READS=$(samtools view -c "$SUB_BAM")
        log_step "[subsample] Subsampled BAM: $ACTUAL_READS reads → $SUB_BAM"
    fi

    INPUT_BAM="$SUB_BAM"
else
    ACTUAL_READS=$(samtools view -c "$BAM")
fi

# ── Clink ─────────────────────────────────────────────────────────────────────
CLINK_OUT="$OUT_DIR/Clink"
CLINK_PREFIX="$CLINK_OUT/$NAME"
CLINK_NPZ="${CLINK_PREFIX}_pileup.npz"
mkdir -p "$CLINK_OUT"

T_CLINK_PILEUP=0; T_CLINK_CITS=0; T_CLINK_CIMS=0; T_CLINK_TOTAL=0

if [[ "$RUN_CLINK" == "true" ]]; then
    log_step ""
    log_step "[Clink] ── pileup ──────────────────────────────────────────"
    t0=$SECONDS
    run_cmd python3 "$CLINK_SCRIPTS/pileup.py" \
        "$INPUT_BAM"     \
        --out     "$CLINK_NPZ" \
        --threads "$THREADS"   \
        --nh      "$NH"        \
        --mapq    "$MAPQ"
    T_CLINK_PILEUP=$((SECONDS - t0))
    log_step "[Clink] pileup done ($(fmt_time $T_CLINK_PILEUP)) → $CLINK_NPZ"

    if [[ "$RUN_CITS" == "true" ]]; then
        log_step "[Clink] ── CITS ─────────────────────────────────────────────"
        t0=$SECONDS
        run_cmd python3 "$CLINK_SCRIPTS/cits.py" \
            --pileup   "$CLINK_NPZ"   \
            --prefix   "$CLINK_PREFIX" \
            --min-cov  "$MIN_COV"     \
            --min-frac "$MIN_FRAC"    \
            --fdr      "$FDR"
        T_CLINK_CITS=$((SECONDS - t0))
        log_step "[Clink] CITS done ($(fmt_time $T_CLINK_CITS))"
    fi

    if [[ "$RUN_CIMS" == "true" ]]; then
        log_step "[Clink] ── CIMS ─────────────────────────────────────────────"
        t0=$SECONDS
        run_cmd python3 "$CLINK_SCRIPTS/cims.py" \
            --pileup   "$CLINK_NPZ"    \
            --prefix   "$CLINK_PREFIX" \
            --min-cov  "$MIN_COV"      \
            --min-frac "$MIN_FRAC"     \
            --fdr      "$FDR"
        T_CLINK_CIMS=$((SECONDS - t0))
        log_step "[Clink] CIMS done ($(fmt_time $T_CLINK_CIMS))"
    fi

    T_CLINK_TOTAL=$(( T_CLINK_PILEUP + T_CLINK_CITS + T_CLINK_CIMS ))
fi

# ── CTK ───────────────────────────────────────────────────────────────────────
CTK_OUT="$OUT_DIR/CTK"
mkdir -p "$CTK_OUT"
T_CTK_TOTAL=0

if [[ "$RUN_CTK" == "true" ]]; then
    log_step ""
    log_step "[CTK] ── starting ──────────────────────────────────────────"
    t0=$SECONDS

    # run_ctk_full_analysis sources the full CTK pipeline:
    #   calmd + parseAlignment.pl → (skip tag2collapse) →
    #   selectRow + getMutationType → CIMS.pl (parallel by chr) → CITS.pl
    # skip_collapse="true" because input BAM is already UMI-deduplicated
    run_ctk_full_analysis \
        "$INPUT_BAM"    \
        "$CTK_OUT"      \
        "${GENOME:-}"   \
        "$CIMS_ITER"    \
        "$FDR"          \
        "$CITS_PVAL"    \
        "25"            \
        "10"            \
        "no"            \
        "$RUN_CIMS"     \
        "$RUN_CITS"     \
        "true"

    T_CTK_TOTAL=$((SECONDS - t0))
    log_step "[CTK] done ($(fmt_time $T_CTK_TOTAL))"
fi

# ── Collect counts ────────────────────────────────────────────────────────────
# Clink outputs use the PREFIX directly
CLINK_CITS_N=$(  count_sites "${CLINK_PREFIX}_truncations.bed")
CLINK_DEL_N=$(   count_sites "${CLINK_PREFIX}_deletions.bed")

# CTK derives its sample name from the BAM filename (same logic as run_ctk_full_analysis)
CTK_SNAME=$(basename "${INPUT_BAM%.bam}" | sed 's/\.Aligned\.sortedByCoord\.out//')
CTK_CITS_N=$(count_sites "$CTK_OUT/CITS/${CTK_SNAME}_CITS.bed")
CTK_DEL_N=$( count_sites "$CTK_OUT/CIMS/${CTK_SNAME}_CIMS_del.txt")
CTK_SUB_N=$( count_sites "$CTK_OUT/CIMS/${CTK_SNAME}_CIMS_sub.txt")

# ── Summary ───────────────────────────────────────────────────────────────────
SUMMARY=$(cat <<EOF

════════════════════════════════════════════════════════════════
  SUMMARY  |  $NAME
  Reads:   $ACTUAL_READS
  FDR:     $FDR  |  CITS p: $CITS_PVAL  |  CTK iter: $CIMS_ITER
════════════════════════════════════════════════════════════════
                         Clink           CTK
  ──────────────────────────────────────────────────────────────
  CITS (truncations):  $(printf '%-15s' "$CLINK_CITS_N")  $CTK_CITS_N
  CIMS deletions:      $(printf '%-15s' "$CLINK_DEL_N")  $CTK_DEL_N
  CIMS substitutions:  $(printf '%-15s' "(see *toC.bed)")  $CTK_SUB_N
  ──────────────────────────────────────────────────────────────
  Pileup/parse time:   $(printf '%-15s' "$(fmt_time $T_CLINK_PILEUP)")  (included in total)
  CIMS time:           $(printf '%-15s' "$(fmt_time $T_CLINK_CIMS)")  (included in total)
  CITS time:           $(printf '%-15s' "$(fmt_time $T_CLINK_CITS)")  (included in total)
  Total time:          $(printf '%-15s' "$(fmt_time $T_CLINK_TOTAL)")  $(fmt_time $T_CTK_TOTAL)
════════════════════════════════════════════════════════════════

Output directories:
  Clink  →  $CLINK_OUT
  CTK    →  $CTK_OUT
  Log    →  $LOG
EOF
)

echo "$SUMMARY" | tee -a "$LOG"
echo "$SUMMARY" > "$OUT_DIR/summary.txt"
