#!/usr/bin/env bash
# clink_restat.sh — Run Clink CITS/CIMS on an existing pileup.npz OR a BAM file
#
# Two modes:
#   --pileup mode  (original): re-run CITS/CIMS from a pre-built pileup.npz.
#                  Skips BAM reading entirely; useful for re-thresholding.
#
#   --bam mode     (new): start from a pre-deduplicated BAM (e.g. a subsampled
#                  benchmark BAM), run pileup.py to build the npz, then CITS/CIMS.
#                  Mirrors ctk_restat.sh interface for apples-to-apples benchmarks.
#                  Does NOT run umi_tools dedup — input BAM must already be deduped.
#
# Clink pipeline (--bam mode):
#   BAM  →  pileup.py (per-position numpy arrays)  →  pileup.npz
#        →  cits.py  (BH FDR on truncation counts)  →  truncations.bed
#        →  cims.py  (BH FDR on deletion/sub counts) →  deletions.bed + sub beds
#
# Usage:
#   # Re-threshold from existing pileup
#   clink_restat.sh --pileup HNRNPC_K562_pileup.npz --fdr 0.001
#
#   # Benchmark: start from a (pre-dedup) BAM
#   clink_restat.sh \
#       --bam  HNRNPC_K562_sub3M.bam \
#       --out  /path/to/clink_sub3M \
#       --name HNRNPC_K562_sub3M \
#       --fdr 0.001 --threads 8

set -euo pipefail

# ── Locate clink scripts relative to this script ─────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLINK_DIR="$SCRIPT_DIR/lib/clink"

# ── Defaults ──────────────────────────────────────────────────────────────────
# Input mode (one of --bam or --pileup required)
BAM=""
NPZ=""

# Output
OUT_DIR=""       # required in --bam mode; derived from npz path in --pileup mode
SAMPLE_NAME=""   # required in --bam mode; derived from npz filename in --pileup mode
PREFIX=""        # auto-derived; can be overridden with --prefix

# Pileup options (--bam mode only)
THREADS="4"
NH="1"           # NH tag threshold: 1 = unique mappers only; 100 = include multi-mappers
MAPQ="255"       # MAPQ filter: 255 = unique (STAR); use 1 or 0 to include multi-mappers

# CITS/CIMS thresholds
FDR="0.05"
MIN_COV="5"
MIN_FRAC="0.05"

# Mode
RUN_CITS="true"
RUN_CIMS="true"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF

Usage: $(basename "$0") (--bam <dedup.bam> --out <dir> --name <str> | --pileup <path.npz>) [options]

Input (one required):
  --bam   <path>       Pre-deduplicated BAM → run pileup.py then CITS/CIMS
  --pileup <path>      Existing pileup.npz  → re-run CITS/CIMS only (no BAM read)

BAM mode options:
  --out   <dir>        Output directory (created if needed)
  --name  <str>        Sample name prefix for output files
  --threads <int>      Threads for pileup.py              (default: $THREADS)
  --nh  <int>          NH tag max for multi-mappers        (default: $NH)
                         1   = unique mappers only
                         100 = include multi-mappers (use if BAM was built with --clink-multi-map)
  --mapq <int>         Minimum MAPQ filter                 (default: $MAPQ)
                         255 = STAR unique mappers
                         1   = include multi-mappers

Pileup mode options:
  --prefix <str>       Output file prefix  (default: auto-derived from --pileup path)

Thresholds:
  --fdr  <float>       FDR threshold for CITS and CIMS     (default: $FDR)
  --min-cov <int>      Minimum read coverage               (default: $MIN_COV)
  --min-frac <float>   Minimum mutation fraction           (default: $MIN_FRAC)

Mode:
  --cits-only          Run CITS (truncation sites) only
  --cims-only          Run CIMS (deletion/substitution sites) only

Other:
  --clink-dir <path>   Path to clink scripts dir           (default: $CLINK_DIR)
  -h, --help           Show this help

Examples:
  # Re-threshold from existing pileup
  $(basename "$0") --pileup HNRNPC_K562_pileup.npz --fdr 0.001

  # Benchmark vs CTK: both tools, same subsampled BAM, same thresholds
  $(basename "$0") \\
      --bam  HNRNPC_K562_sub3M.bam \\
      --out  /path/to/clink_sub3M \\
      --name HNRNPC_K562_sub3M \\
      --fdr 0.001 --threads 8 --nh 1

EOF
    exit 0
}

# ── Parse args ────────────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)       BAM="$2";        shift 2 ;;
        --pileup)    NPZ="$2";        shift 2 ;;
        --out)       OUT_DIR="$2";    shift 2 ;;
        --name)      SAMPLE_NAME="$2"; shift 2 ;;
        --prefix)    PREFIX="$2";     shift 2 ;;
        --threads)   THREADS="$2";    shift 2 ;;
        --nh)        NH="$2";         shift 2 ;;
        --mapq)      MAPQ="$2";       shift 2 ;;
        --fdr)       FDR="$2";        shift 2 ;;
        --min-cov)   MIN_COV="$2";    shift 2 ;;
        --min-frac)  MIN_FRAC="$2";   shift 2 ;;
        --cits-only) RUN_CIMS="false"; shift ;;
        --cims-only) RUN_CITS="false"; shift ;;
        --clink-dir) CLINK_DIR="$2";  shift 2 ;;
        -h|--help)   usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ ! -d "$CLINK_DIR" ]] && { echo "ERROR: clink dir not found: $CLINK_DIR" >&2; exit 1; }

if [[ -n "$BAM" ]] && [[ -n "$NPZ" ]]; then
    echo "ERROR: specify --bam OR --pileup, not both" >&2; exit 1
fi

if [[ -z "$BAM" ]] && [[ -z "$NPZ" ]]; then
    echo "ERROR: one of --bam or --pileup is required" >&2; exit 1
fi

# ── BAM mode: run pileup.py first ─────────────────────────────────────────────
if [[ -n "$BAM" ]]; then
    [[ ! -f "$BAM" ]]      && { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
    [[ -z "$OUT_DIR" ]]    && { echo "ERROR: --out is required with --bam" >&2; exit 1; }
    [[ -z "$SAMPLE_NAME" ]] && { echo "ERROR: --name is required with --bam" >&2; exit 1; }

    mkdir -p "$OUT_DIR"
    NPZ="${OUT_DIR}/${SAMPLE_NAME}_pileup.npz"
    PREFIX="${OUT_DIR}/${SAMPLE_NAME}"

    echo "════════════════════════════════════════════════════════"
    echo "  Clink re-stat  [BAM mode]"
    echo "  BAM:      $BAM"
    echo "  Name:     $SAMPLE_NAME"
    echo "  Output:   $OUT_DIR"
    echo "  Pileup:   $NPZ"
    echo "  Threads:  $THREADS  |  NH≤$NH  |  MAPQ≥$MAPQ"
    echo "  FDR:      $FDR  |  min-cov: $MIN_COV  |  min-frac: $MIN_FRAC"
    echo "  CITS:     $RUN_CITS  |  CIMS: $RUN_CIMS"
    echo "════════════════════════════════════════════════════════"

    echo "[pileup] Building pileup from BAM..."
    python3 "$CLINK_DIR/pileup.py" \
        "$BAM" \
        --out     "$NPZ" \
        --threads "$THREADS" \
        --nh      "$NH" \
        --mapq    "$MAPQ"

    if [[ ! -s "$NPZ" ]]; then
        echo "ERROR: pileup.py failed or produced empty NPZ: $NPZ" >&2; exit 1
    fi
    echo "[pileup] Done → $NPZ"

# ── Pileup mode: use existing npz ─────────────────────────────────────────────
else
    [[ ! -f "$NPZ" ]] && { echo "ERROR: pileup file not found: $NPZ" >&2; exit 1; }

    # Auto-derive prefix and out dir from npz path if not provided
    if [[ -z "$PREFIX" ]]; then
        PREFIX="${NPZ%_pileup.npz}"
        [[ "$PREFIX" == "$NPZ" ]] && PREFIX="${NPZ%.npz}"
    fi
    [[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$PREFIX")"

    echo "════════════════════════════════════════════════════════"
    echo "  Clink re-stat  [pileup mode]"
    echo "  NPZ:      $NPZ"
    echo "  Prefix:   $PREFIX"
    echo "  FDR:      $FDR  |  min-cov: $MIN_COV  |  min-frac: $MIN_FRAC"
    echo "  CITS:     $RUN_CITS  |  CIMS: $RUN_CIMS"
    echo "════════════════════════════════════════════════════════"
fi

# ── CITS ──────────────────────────────────────────────────────────────────────
if [[ "$RUN_CITS" == "true" ]]; then
    echo "[CITS] Calling truncation sites..."
    python3 "$CLINK_DIR/cits.py" \
        --pileup   "$NPZ"      \
        --prefix   "$PREFIX"   \
        --min-cov  "$MIN_COV"  \
        --min-frac "$MIN_FRAC" \
        --fdr      "$FDR"

    OUT_BED="${PREFIX}_truncations.bed"
    if [[ -f "$OUT_BED" ]]; then
        N=$(( $(grep -c '' "$OUT_BED" 2>/dev/null || echo 0) - 1 ))
        [[ $N -lt 0 ]] && N=0
        echo "[CITS] Done — $N significant truncation sites → $OUT_BED"
    fi
fi

# ── CIMS ──────────────────────────────────────────────────────────────────────
if [[ "$RUN_CIMS" == "true" ]]; then
    echo "[CIMS] Calling deletion/substitution sites..."
    python3 "$CLINK_DIR/cims.py" \
        --pileup   "$NPZ"      \
        --prefix   "$PREFIX"   \
        --min-cov  "$MIN_COV"  \
        --min-frac "$MIN_FRAC" \
        --fdr      "$FDR"

    OUT_DEL="${PREFIX}_deletions.bed"
    if [[ -f "$OUT_DEL" ]]; then
        N=$(( $(grep -c '' "$OUT_DEL" 2>/dev/null || echo 0) - 1 ))
        [[ $N -lt 0 ]] && N=0
        echo "[CIMS] Done — $N significant deletion sites → $OUT_DEL"
    fi
fi

echo "════════════════════════════════════════════════════════"
echo "  Complete. Output files:"
ls -lh "${PREFIX}"_*.bed "${PREFIX}"_*.npz 2>/dev/null \
    | awk '{print "  " $NF, $5}' || true
echo "════════════════════════════════════════════════════════"
