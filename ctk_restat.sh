#!/usr/bin/env bash
# ctk_restat.sh — Run CTK CITS/CIMS on an existing dedup BAM
#
# Skips alignment and collapse. Sources CLIPittyClip modules and calls
# run_ctk_full_analysis() directly on a dedup'd BAM file, producing output
# that matches the structure expected by PEAKittyPeak.sh --ctk-dir.
#
# Usage:
#   ctk_restat.sh --bam <dedup.bam> --out <output_dir> [options]
#
# Examples:
#   # Single sample
#   ctk_restat.sh \
#       --bam  5_Clink/ENCFF040BGS/ENCFF040BGS_dedup.bam \
#       --out  6_CTK/ENCFF040BGS \
#       --name ENCFF040BGS \
#       --cits-pval 0.01 --cims-fdr 0.01
#
#   # Group merged BAM
#   ctk_restat.sh \
#       --bam  5_Clink/GROUP_ELAVL1_HepG2/HNRNPC_HepG2_merged_dedup.bam \
#       --out  6_CTK/GROUP_HepG2 \
#       --name HNRNPC_HepG2 \
#       --cits-pval 0.01 --cims-fdr 0.01

set -uo pipefail   # note: NO -e — modules.sh functions return non-zero on warnings;
                   # we handle errors explicitly below instead

# ── Locate lib dir relative to this script ────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_DIR="$SCRIPT_DIR/lib"

# ── Defaults ──────────────────────────────────────────────────────────────────
BAM=""
OUT_DIR=""
SAMPLE_NAME=""
GENOME_FASTA=""
CIMS_FDR="0.05"
CITS_PVAL="0.05"
CIMS_ITER="10"
CITS_GAP="25"
THREADS="4"
RUN_CITS="true"
RUN_CIMS="true"
RUN_MOTIF="no"     # off by default — slow and needs genome FASTA
SKIP_COLLAPSE="true"  # default true: input is always a pre-deduplicated BAM
VERBOSE="true"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF

Usage: $(basename "$0") --bam <dedup.bam> --out <dir> [options]

Required:
  --bam  <path>        Deduplicated BAM file (sorted + indexed)
  --out  <dir>         Output directory (created if needed)

Naming:
  --name <str>         Sample name for output file prefixes
                       (default: BAM basename without .bam)

Thresholds:
  --cits-pval <float>  CITS p-value threshold        (default: $CITS_PVAL)
  --cims-fdr  <float>  CIMS FDR threshold            (default: $CIMS_FDR)
  --cims-iter <int>    CIMS permutation iterations   (default: $CIMS_ITER)
  --cits-gap  <int>    CITS gap merging window       (default: $CITS_GAP)

Mode:
  --no-cits            Skip CITS (truncation sites)
  --no-cims            Skip CIMS (deletion/substitution sites)
  --run-motif          Enable motif analysis (requires --genome)
  --with-collapse      Run position-based tag2collapse before CIMS/CITS
                       (default: skipped — input BAM is pre-deduplicated;
                       running tag2collapse on UMI-dedup data collapses
                       genuine multi-cell crosslinks to depth=1, breaking CITS)

Other:
  --genome <path>      Reference FASTA (required only with --run-motif)
  --threads <int>      Parallel threads                (default: $THREADS)
  -h, --help           Show this help

Output structure (matches PEAKittyPeak --ctk-dir expectation):
  <out>/CITS/<name>_CITS.txt
  <out>/CIMS/<name>_CIMS_del.txt
  <out>/CIMS/<name>_CIMS_sub.txt

EOF
    exit 0
}

# ── Parse args ────────────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)        BAM="$2";          shift 2 ;;
        --out)        OUT_DIR="$2";      shift 2 ;;
        --name)       SAMPLE_NAME="$2";  shift 2 ;;
        --genome)     GENOME_FASTA="$2"; shift 2 ;;
        --cims-fdr)   CIMS_FDR="$2";     shift 2 ;;
        --cits-pval)  CITS_PVAL="$2";    shift 2 ;;
        --cims-iter)  CIMS_ITER="$2";    shift 2 ;;
        --cits-gap)   CITS_GAP="$2";     shift 2 ;;
        --threads)    THREADS="$2";      shift 2 ;;
        --no-cits)       RUN_CITS="false";      shift ;;
        --no-cims)       RUN_CIMS="false";      shift ;;
        --run-motif)     RUN_MOTIF="yes";       shift ;;
        --with-collapse) SKIP_COLLAPSE="false"; shift ;;
        -h|--help)       usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$BAM" ]]           && { echo "ERROR: --bam is required" >&2; exit 1; }
[[ ! -f "$BAM" ]]         && { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -z "$OUT_DIR" ]]       && { echo "ERROR: --out is required" >&2; exit 1; }
[[ ! -d "$LIB_DIR" ]]     && { echo "ERROR: lib dir not found: $LIB_DIR" >&2; exit 1; }

if [[ "$RUN_MOTIF" == "yes" && -z "$GENOME_FASTA" ]]; then
    echo "ERROR: --run-motif requires --genome <fasta>" >&2; exit 1
fi

# Auto-derive sample name from BAM if not given
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME="$(basename "${BAM%.bam}")"
fi

# ── Source CLIPittyClip modules ───────────────────────────────────────────────
mkdir -p "$OUT_DIR"
export LOG_FILE="${OUT_DIR}/ctk_restat.log"
export VERBOSE
export RUN_CIMS RUN_CITS   # checked as env vars inside run_parse_alignment

# modules.sh writes the SAM temp file next to the input BAM by default, which
# may not be writable. Symlink the BAM into the output dir so the SAM lands there.
BAM_LINK="${OUT_DIR}/$(basename "$BAM")"
ln -sf "$BAM" "$BAM_LINK"
ln -sf "${BAM}.bai" "${BAM_LINK}.bai" 2>/dev/null || true
BAM="$BAM_LINK"

# shellcheck source=lib/utils.sh
source "$LIB_DIR/utils.sh"
# shellcheck source=lib/modules.sh
source "$LIB_DIR/modules.sh"

# ── Summary ───────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════"
echo "  CTK re-stat"
echo "  BAM:      $BAM"
echo "  Sample:   $SAMPLE_NAME"
echo "  Output:   $OUT_DIR"
echo "  CITS p:   $CITS_PVAL  |  CIMS FDR: $CIMS_FDR"
echo "  CITS:     $RUN_CITS   |  CIMS:     $RUN_CIMS"
echo "  Collapse: $([ "$SKIP_COLLAPSE" == "true" ] && echo "skipped (pre-dedup BAM)" || echo "enabled")"
echo "  Log:      $LOG_FILE"
echo "════════════════════════════════════════════════════════"

# ── Run ───────────────────────────────────────────────────────────────────────
run_ctk_full_analysis \
    "$BAM"            \
    "$OUT_DIR"        \
    "$GENOME_FASTA"   \
    "$CIMS_ITER"      \
    "$CIMS_FDR"       \
    "$CITS_PVAL"      \
    "$CITS_GAP"       \
    "10"              \
    "$RUN_MOTIF"      \
    "$RUN_CIMS"       \
    "$RUN_CITS"       \
    "$SKIP_COLLAPSE"
CTK_EXIT=$?

# Clean up symlink
rm -f "$BAM_LINK" "${BAM_LINK}.bai" 2>/dev/null || true

if [[ $CTK_EXIT -ne 0 ]]; then
    echo "════════════════════════════════════════════════════════" >&2
    echo "  ERROR: CTK analysis failed (exit $CTK_EXIT)" >&2
    echo "  Last 20 lines of log:" >&2
    tail -20 "$LOG_FILE" >&2
    echo "════════════════════════════════════════════════════════" >&2
    exit $CTK_EXIT
fi

echo "════════════════════════════════════════════════════════"
echo "  Complete. Output files:"
find "$OUT_DIR" \( -name "*.txt" -o -name "*.bed" \) 2>/dev/null \
    | sort | awk '{print "  "$0}' || true
echo "════════════════════════════════════════════════════════"
