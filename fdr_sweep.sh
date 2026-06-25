#!/usr/bin/env bash
# fdr_sweep.sh — sweep significance thresholds and compute peak overlap fraction
#
# Works with:
#   Clink truncations.bed / deletions.bed  — qvalue in col 11
#   CTK CITS.txt                           — p-value in col 5
#   CTK CIMS_del/sub.txt                   — p-value col 9, FDR col 10
#
# For each threshold, reports: how many sites pass, how many fall inside a peak,
# and what fraction that represents. Output is a TSV suitable for plotting.
#
# Usage:
#   fdr_sweep.sh --sites <file> --peaks <bed> [options]
#
# Examples:
#   # Clink CITS — auto-detects col 11
#   fdr_sweep.sh \
#       --sites  5_Clink/GROUP_HNRNPC_HepG2/HNRNPC_HepG2_truncations.bed \
#       --peaks  4_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed \
#       --label  HNRNPC_HepG2_Clink_CITS
#
#   # CTK CITS — p-value in col 5
#   fdr_sweep.sh \
#       --sites  6_CTK/GROUP_HNRNPC_HepG2/CITS/HNRNPC_HepG2_CITS.txt \
#       --peaks  4_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed \
#       --val-col 5 --label HNRNPC_HepG2_CTK_CITS
#
#   # CTK CIMS (FDR col 10)
#   fdr_sweep.sh \
#       --sites  6_CTK/GROUP_HNRNPC_HepG2/CIMS/HNRNPC_HepG2_CIMS_del.txt \
#       --peaks  4_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed \
#       --val-col 10 --label HNRNPC_HepG2_CTK_CIMS_del

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
SITES=""
PEAKS=""
LABEL=""
VAL_COL=""           # auto-detect if not set
STRAND="-s"          # strand-specific overlap; use --no-strand to disable
SLOP_5=0             # nt to extend peaks upstream   (5' / left  on + strand)
SLOP_3=0             # nt to extend peaks downstream (3' / right on + strand)
# Log-scale sweep from loose → strict
THRESHOLDS="0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00001"
OUT=""               # default: stdout

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF

Usage: $(basename "$0") --sites <file> --peaks <bed> [options]

Required:
  --sites  <file>      CITS/CIMS site file (Clink BED or CTK txt)
  --peaks  <bed>       Peak BED file to intersect against

Value column:
  --val-col <int>      Column (1-based) containing p/q-value to sweep
                       Auto-detected if not set:
                         11 cols → col 11 (Clink qvalue)
                         filename contains CIMS → col 10 (CTK CIMS FDR)
                         otherwise            → col 5  (CTK CITS p-value)

Thresholds:
  --thresholds <csv>   Comma-separated threshold list (default: log scale
                       0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,
                       0.0005,0.0002,0.0001,0.00001)

Peak extension (strand-aware via bedtools slop -s):
  --slop-5 <int>       Extend peaks upstream   (5' end) by N nt (default: 0)
  --slop-3 <int>       Extend peaks downstream (3' end) by N nt (default: 0)
                       Example: --slop-5 25 extends each peak 25 nt toward the
                       5' end (capturing CITS that sit upstream of peak center)

Overlap:
  --no-strand          Ignore strand for peak overlap (default: strand-specific)

Output:
  --label  <str>       Label column in output (default: sites filename basename)
  --out    <file>      Write TSV to file instead of stdout

  -h, --help           Show this help

Output columns:
  label  threshold  n_sites  n_in_peak  n_out_peak  pct_in_peak

EOF
    exit 0
}

# ── Parse args ────────────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites)       SITES="$2";       shift 2 ;;
        --peaks)       PEAKS="$2";       shift 2 ;;
        --label)       LABEL="$2";       shift 2 ;;
        --val-col)     VAL_COL="$2";     shift 2 ;;
        --thresholds)  THRESHOLDS="$2";  shift 2 ;;
        --slop-5)      SLOP_5="$2";      shift 2 ;;
        --slop-3)      SLOP_3="$2";      shift 2 ;;
        --no-strand)   STRAND="";        shift ;;
        --out)         OUT="$2";         shift 2 ;;
        -h|--help)     usage ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$SITES" ]]   && { echo "ERROR: --sites required" >&2; exit 1; }
[[ ! -f "$SITES" ]] && { echo "ERROR: sites file not found: $SITES" >&2; exit 1; }
[[ -z "$PEAKS" ]]   && { echo "ERROR: --peaks required" >&2; exit 1; }
[[ ! -f "$PEAKS" ]] && { echo "ERROR: peaks file not found: $PEAKS" >&2; exit 1; }
command -v bedtools &>/dev/null || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }

# ── Auto-detect value column ──────────────────────────────────────────────────
if [[ -z "$VAL_COL" ]]; then
    # Count columns in first non-comment line
    NCOLS=$(grep -v "^#" "$SITES" | head -1 | awk '{print NF}')
    BASENAME=$(basename "$SITES")
    if [[ "$NCOLS" -ge 11 ]]; then
        VAL_COL=11
        echo "[INFO] Auto-detected Clink format ($NCOLS cols) → using col 11 (qvalue)" >&2
    elif echo "$BASENAME" | grep -qi "CIMS"; then
        VAL_COL=9
        echo "[INFO] Auto-detected CTK CIMS format → using col 9 (FDR)" >&2
    else
        VAL_COL=5
        echo "[INFO] Auto-detected CTK CITS format → using col 5 (p-value)" >&2
        echo "[INFO] NOTE: if sites file is a CITS.bed, p-value is in col 4 name field;" >&2
        echo "[INFO]       use --val-col name_p or use BiTS/fdr_sweep.sh for auto-handling." >&2
    fi
fi

# ── Default label ─────────────────────────────────────────────────────────────
[[ -z "$LABEL" ]] && LABEL=$(basename "$SITES" | sed 's/\.\(bed\|txt\)$//')

# ── Setup temp dir ────────────────────────────────────────────────────────────
TMP=$(mktemp -d /tmp/fdr_sweep.XXXXXX)
trap "rm -rf $TMP" EXIT

# Need genome chrom sizes for bedtools slop
CHROMSIZES="$TMP/chrom.sizes"
# Derive from the peaks file itself (max end per chrom — safe lower bound)
awk 'BEGIN{OFS="\t"} {if($3>s[$1]) s[$1]=$3} END{for(c in s) print c, s[c]+1000000}' \
    "$PEAKS" > "$CHROMSIZES"

# Pre-sort and optionally slop peaks
PEAKS_SORTED="$TMP/peaks.bed"
if [[ "$SLOP_5" -gt 0 || "$SLOP_3" -gt 0 ]]; then
    sort -k1,1 -k2,2n "$PEAKS" \
        | bedtools slop -i - -g "$CHROMSIZES" -s -l "$SLOP_5" -r "$SLOP_3" \
        | sort -k1,1 -k2,2n > "$PEAKS_SORTED"
    echo "[INFO] Peak extension: upstream(5') +${SLOP_5} nt, downstream(3') +${SLOP_3} nt" >&2
else
    sort -k1,1 -k2,2n "$PEAKS" > "$PEAKS_SORTED"
fi

# Strip header from sites, sort once
SITES_CLEAN="$TMP/sites_clean.bed"
grep -v "^#" "$SITES" | sort -k1,1 -k2,2n > "$SITES_CLEAN"
TOTAL_IN_FILE=$(wc -l < "$SITES_CLEAN")

echo "[INFO] Sites file: $SITES" >&2
echo "[INFO] Total sites in file: $TOTAL_IN_FILE" >&2
echo "[INFO] Value column: $VAL_COL" >&2
echo "[INFO] Peak file: $PEAKS ($(wc -l < $PEAKS_SORTED) peaks)" >&2
echo "[INFO] Strand-specific overlap: ${STRAND:-(none)}" >&2
echo "" >&2

# ── Output header ─────────────────────────────────────────────────────────────
HEADER="label\tthreshold\tslop_5\tslop_3\tn_sites\tn_in_peak\tn_out_peak\tpct_in_peak"
if [[ -n "$OUT" ]]; then
    printf "%b\n" "$HEADER" > "$OUT"
else
    printf "%b\n" "$HEADER"
fi

# ── Sweep ─────────────────────────────────────────────────────────────────────
IFS=',' read -ra THRESH_LIST <<< "$THRESHOLDS"

for T in "${THRESH_LIST[@]}"; do
    FILTERED="$TMP/filtered.bed"

    # Filter by threshold — awk handles scientific notation
    awk -v col="$VAL_COL" -v thresh="$T" \
        'BEGIN{OFS="\t"} NR>0 && $col+0 <= thresh+0 {print}' \
        "$SITES_CLEAN" > "$FILTERED"

    N_TOTAL=$(wc -l < "$FILTERED")

    if [[ "$N_TOTAL" -eq 0 ]]; then
        ROW="${LABEL}\t${T}\t${SLOP_5}\t${SLOP_3}\t0\t0\t0\tNA"
    else
        N_IN=$(bedtools intersect -a "$FILTERED" -b "$PEAKS_SORTED" $STRAND -u 2>/dev/null | wc -l)
        N_OUT=$(( N_TOTAL - N_IN ))
        PCT=$(awk "BEGIN{printf \"%.2f\", 100*${N_IN}/${N_TOTAL}}")
        ROW="${LABEL}\t${T}\t${SLOP_5}\t${SLOP_3}\t${N_TOTAL}\t${N_IN}\t${N_OUT}\t${PCT}"
    fi

    printf "%b\n" "$ROW"
    if [[ -n "$OUT" ]]; then
        printf "%b\n" "$ROW" >> "$OUT"
    fi

    echo "[INFO] threshold=$T  sites=$N_TOTAL  in_peak=${N_IN:-0}  pct=${PCT:-NA}%" >&2
done

echo "" >&2
echo "[INFO] Done." >&2
