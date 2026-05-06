#!/usr/bin/env bash
# benchmark_pipeline.sh — CTK vs Clink pipeline speed + RAM benchmark
#
# Usage:
#   ./benchmark_pipeline.sh -i <input.bam> -o <output_dir> [options]
#
# Options:
#   -i   Input BAM (sorted, indexed)
#   -o   Output directory
#   -n   Number of reads to extract (default: 1,000,000)
#   -t   Threads for pileup.py chromosome-level parallelism (default: 1)
#   -u   UMI length for Clink collapse (default: auto-detect)
#   -f   Genome FASTA for samtools calmd (optional but recommended for CTK CIMS)
#   --cims-iter  CIMS permutation iterations (default: 5; use 200+ for reliable FDR)
#   --cims-m     Min mutations per site to report (default: 1; post-filters CIMS output)
#   --cims-k     Min coverage per site to report  (default: 1; post-filters CIMS output)
#   --fdr        FDR threshold for CIMS/CITS (default: 0.05)

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
N_READS=1000000
INPUT_BAM=""
OUTDIR=""
THREADS=1
UMI_LEN="-1"
GENOME_FASTA=""
CIMS_ITER=5
CIMS_M=1
CIMS_K=1
FDR=0.05

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLINK_DIR="${SCRIPT_DIR}/lib/clink"
export PERL5LIB="${CONDA_PREFIX:-}/lib/czplib:${PERL5LIB:-}"
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
    echo "  -u            UMI length for Clink (default: auto-detect from read names)"
    echo "  -f            Genome FASTA for samtools calmd (recommended for CTK CIMS)"
    echo "  --cims-iter   CIMS permutation iterations (default: 5; use 200+ for reliable FDR)"
    echo "  --cims-m      CIMS minimum mutations per site (default: 1)"
    echo "  --cims-k      CIMS minimum coverage per site (default: 1)"
    echo "  --fdr         FDR/p-value threshold (default: 0.05)"
    echo ""
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) INPUT_BAM="$2";   shift 2 ;;
        -o) OUTDIR="$2";      shift 2 ;;
        -n) N_READS="$2";     shift 2 ;;
        -t) THREADS="$2";     shift 2 ;;
        -u) UMI_LEN="$2";     shift 2 ;;
        -f) GENOME_FASTA="$2"; shift 2 ;;
        --cims-iter) CIMS_ITER="$2"; shift 2 ;;
        --cims-m)    CIMS_M="$2";    shift 2 ;;
        --cims-k)    CIMS_K="$2";    shift 2 ;;
        --fdr) FDR="$2";      shift 2 ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
done

[[ -z "$INPUT_BAM"  ]] && { echo "[ERROR] -i <input.bam> required"; usage; }
[[ -z "$OUTDIR"     ]] && { echo "[ERROR] -o <output_dir> required"; usage; }
[[ ! -f "$INPUT_BAM" ]] && { echo "[ERROR] BAM not found: $INPUT_BAM"; exit 1; }
command -v gtime          &>/dev/null || { echo "[ERROR] gtime not found (brew install gnu-time)"; exit 1; }
command -v samtools       &>/dev/null || { echo "[ERROR] samtools not in PATH"; exit 1; }
command -v parseAlignment.pl &>/dev/null || { echo "[ERROR] parseAlignment.pl not in PATH"; exit 1; }
command -v CIMS.pl        &>/dev/null || { echo "[ERROR] CIMS.pl not in PATH"; exit 1; }
command -v CITS.pl        &>/dev/null || { echo "[ERROR] CITS.pl not in PATH"; exit 1; }
[[ -f "${CLINK_DIR}/pileup.py" ]] || { echo "[ERROR] pileup.py not found: ${CLINK_DIR}"; exit 1; }

CTK_DIR="${OUTDIR}/CTK"
CLK_DIR="${OUTDIR}/Clink"
mkdir -p "$CTK_DIR" "$CLK_DIR"

# Remove stale CIMS/CITS output so old files don't pollute counts if a step fails
rm -f "${CTK_DIR}/CIMS_del.txt" "${CTK_DIR}/CIMS_sub.txt" "${CTK_DIR}/CITS.txt"
rm -f "${CLK_DIR}/sample"_*.bed

REPORT="${OUTDIR}/benchmark_report.txt"

# ── Helpers ───────────────────────────────────────────────────────────────────

# fmt_num: add thousands commas (python — avoids BSD/GNU sed differences)
fmt_num() { python3 -c "print(f'{int($1):,}')"; }

# count_data_lines FILE: lines minus 1 header (0 if missing/empty)
count_data() {
    local f="$1"
    [[ ! -s "$f" ]] && echo 0 && return
    local n
    n=$(wc -l < "$f")
    echo $(( n > 0 ? n - 1 : 0 ))
}

# plain_count FILE: raw line count
plain_count() {
    [[ ! -s "$1" ]] && echo 0 && return
    wc -l < "$1"
}

# bam_count FILE: mapped read count
bam_count() {
    samtools view -c -F 4 "$1" 2>/dev/null || echo 0
}

# run_timed LABEL GTIME_FILE CMD...: run CMD under gtime, silence tool output
# All stdout+stderr from the command goes to a .log file alongside the gtime file.
# gtime stats are written to gtime's own stderr → mixed into the log, then extracted.
run_timed() {
    local label="$1"
    local gt_file="$2"
    local log_file="${gt_file%.txt}.log"
    shift 2
    printf "    %-35s " "$label"

    # Run: redirect cmd stdout+stderr to log; gtime writes its stats to its own
    # stderr which also lands in the log (fd2 → log via 2>&1).
    gtime -v "$@" >"$log_file" 2>&1 || true

    # Extract gtime stats lines (distinctive format, safe to grep)
    grep -E "Elapsed \(wall clock\)|Maximum resident|User time|System time" \
        "$log_file" > "$gt_file" 2>/dev/null || true

    # Parse
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

# accumulate_time WALL_S: add to TOTAL_WALL
accumulate() {
    TOTAL_WALL=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_WALL:-0} + $1}")
    # track peak RAM
    local rss="$2"
    PEAK_RSS=$(awk "BEGIN { a=${PEAK_RSS:-0}; b=${rss}; print (a>b?a:b) }")
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  CLIPittyClip — CTK vs Clink Pipeline Benchmark"
echo "════════════════════════════════════════════════════════════════"
echo "  Input BAM   : $INPUT_BAM"
echo "  Reads       : $(fmt_num $N_READS)"
echo "  Threads     : $THREADS  (pileup.py)"
echo "  CIMS iter   : $CIMS_ITER   min-mut: $CIMS_M   min-cov: $CIMS_K   FDR: $FDR"
echo "  Genome FASTA: ${GENOME_FASTA:-none (calmd skipped)}"
echo "  Output dir  : $OUTDIR"
echo "════════════════════════════════════════════════════════════════"
echo ""

# ── Step 0: Subset BAM ────────────────────────────────────────────────────────
echo "[0] Extracting $(fmt_num $N_READS) reads from BAM..."

SUBSET_BAM="${OUTDIR}/subset_${N_READS}.bam"
SUBSET_SORTED="${OUTDIR}/subset_sorted.bam"

# Disable pipefail here: head closes the pipe early once it has N lines,
# causing samtools view to receive SIGPIPE (exit 141). That is expected behaviour.
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
if [[ "$UMI_LEN" != "-1" ]]; then
    CLINK_COLLAPSE_ARGS="$CLINK_COLLAPSE_ARGS --umi-len $UMI_LEN"
fi
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

# ── Clink 2b: raw pileup signal counts (positions with any signal, min-cov=5) ─
CLK_RAW_DEL_POS=0; CLK_RAW_SUB_POS=0; CLK_RAW_TRUNC_POS=0
if [[ -f "$NPZ" ]]; then
    read -r CLK_RAW_DEL_POS CLK_RAW_SUB_POS CLK_RAW_TRUNC_POS < <(python3 - "$NPZ" <<'PYEOF'
import sys, numpy as np
data = np.load(sys.argv[1])
# New stranded format: keys are {chrom}__{fwd|rev}__{field}
dels = subs = trunc = 0
for k in data.files:
    parts = k.split('__')
    if len(parts) < 3 or parts[1] not in ('fwd', 'rev'):
        continue
    field = parts[2]
    if field == 'coverage':
        pass  # needed for mask below
    pfx = f'{parts[0]}__{parts[1]}'
    if field != 'coverage':
        continue
    cov = data[k]
    mask = cov >= 5
    d_key = f'{pfx}__deletions'
    t_key = f'{pfx}__truncations'
    if d_key in data.files:
        d = data[d_key]
        dels  += int(np.sum(d[mask] > 0))
    if t_key in data.files:
        t = data[t_key]
        trunc += int(np.sum(t[mask] > 0))
    sub_keys = [sk for sk in data.files if sk.startswith(f'{pfx}__sub__')]
    if sub_keys:
        sub_any = np.zeros(int(mask.sum()), dtype=bool)
        for sk in sub_keys:
            sub_any |= (data[sk][mask] > 0)
        subs += int(np.sum(sub_any))
print(dels, subs, trunc)
PYEOF
)
fi

# ── Clink 3: cits.py ─────────────────────────────────────────────────────────
CLK_PREFIX="${CLK_DIR}/sample"
run_timed "cits.py" "${CLK_DIR}/gt_cits.txt" \
    python3 "${CLINK_DIR}/cits.py" \
        --pileup "$NPZ" --prefix "$CLK_PREFIX" \
        --min-cov 5 --min-frac 0.05 --min-signal "$CIMS_M" --fdr "$FDR"
accumulate "$LAST_WALL" "$LAST_RSS"

CLK_TRUNC=$(count_data "${CLK_PREFIX}_truncations.bed")

# ── Clink 4: cims.py ─────────────────────────────────────────────────────────
run_timed "cims.py" "${CLK_DIR}/gt_cims.txt" \
    python3 "${CLINK_DIR}/cims.py" \
        --pileup "$NPZ" --prefix "$CLK_PREFIX" \
        --min-cov 5 --min-frac 0.05 --min-signal "$CIMS_M" --fdr "$FDR"
accumulate "$LAST_WALL" "$LAST_RSS"

CLK_DELS=$(count_data "${CLK_PREFIX}_deletions.bed")

# Count all substitution beds and find T>C specifically
CLK_SUBS_TOTAL=0
for sub_bed in "${CLK_PREFIX}"_*to*.bed; do
    [[ -f "$sub_bed" ]] || continue
    n=$(count_data "$sub_bed")
    CLK_SUBS_TOTAL=$(( CLK_SUBS_TOTAL + n ))
done
CLK_TtoC=$(count_data "${CLK_PREFIX}_TtoC.bed")

CLK_TOTAL_WALL="$TOTAL_WALL"
CLK_PEAK_RSS="$PEAK_RSS"

echo "    ─────────────────────────────────────────────────────────"
printf "    %-35s %8ss  %8s MB\n" "TOTAL" "$CLK_TOTAL_WALL" "$CLK_PEAK_RSS"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
echo "── CTK Pipeline ─────────────────────────────────────────────────────────"
echo "    Step                                Wall       Peak RAM"
echo "    ─────────────────────────────────────────────────────────"

TOTAL_WALL=0; PEAK_RSS=0

# ── CTK 1: BAM → SAM (+ optional calmd) ──────────────────────────────────────
SAM_FILE="${CTK_DIR}/subset.sam"
if [[ -n "$GENOME_FASTA" && -f "$GENOME_FASTA" ]]; then
    run_timed "samtools calmd → SAM" "${CTK_DIR}/gt_calmd.txt" \
        bash -c "samtools calmd -S '$SUBSET_SORTED' '$GENOME_FASTA' > '$SAM_FILE' 2>/dev/null"
else
    run_timed "samtools view → SAM" "${CTK_DIR}/gt_sam.txt" \
        bash -c "samtools view -h '$SUBSET_SORTED' > '$SAM_FILE'"
fi
accumulate "$LAST_WALL" "$LAST_RSS"

# ── CTK 2: parseAlignment.pl ──────────────────────────────────────────────────
TAGS_BED="${CTK_DIR}/tags.bed"
MUTATIONS_TXT="${CTK_DIR}/mutations.txt"
run_timed "parseAlignment.pl" "${CTK_DIR}/gt_parse.txt" \
    parseAlignment.pl --map-qual 1 --min-len 16 \
        --mutation-file "$MUTATIONS_TXT" "$SAM_FILE" "$TAGS_BED"
accumulate "$LAST_WALL" "$LAST_RSS"
rm -f "$SAM_FILE"

CTK_TAGS=$(plain_count "$TAGS_BED")

# ── CTK 3: tag2collapse.pl ────────────────────────────────────────────────────
COLLAPSED_BED="${CTK_DIR}/collapsed.bed"
run_timed "tag2collapse.pl" "${CTK_DIR}/gt_collapse.txt" \
    tag2collapse.pl --keep-tag-name "$TAGS_BED" "$COLLAPSED_BED"
accumulate "$LAST_WALL" "$LAST_RSS"

CTK_COLLAPSED=$(plain_count "$COLLAPSED_BED")

# ── CTK 4: selectRow + getMutationType ───────────────────────────────────────
MATCHED_MUT="${CTK_DIR}/mutations_matched.txt"
DEL_BED="${CTK_DIR}/deletions.bed"
SUB_BED="${CTK_DIR}/substitutions.bed"
run_timed "selectRow + getMutationType" "${CTK_DIR}/gt_muttype.txt" \
    bash -c "selectRow.pl -q 3 -f 3 '$MUTATIONS_TXT' '$COLLAPSED_BED' > '$MATCHED_MUT' && \
             getMutationType.pl -t del '$MATCHED_MUT' '$DEL_BED' && \
             getMutationType.pl -t sub '$MATCHED_MUT' '$SUB_BED'"
accumulate "$LAST_WALL" "$LAST_RSS"

CTK_DELS=$(plain_count "$DEL_BED")
CTK_SUBS=$(plain_count "$SUB_BED")

# ── CTK 5: CIMS.pl — deletions ───────────────────────────────────────────────
# CIMS.pl output cols: chrom chromStart chromEnd name score strand tagNumber(k) mutationFreq(m) FDR count
# col7=coverage(k)  col8=mutations(m)
CIMS_DEL="${CTK_DIR}/CIMS_del.txt"
if [[ -s "$DEL_BED" ]]; then
    CACHE_CIMS_DEL=$(mktemp -u "${TMPDIR:-/tmp}/cims_del.XXXXXX")
    run_timed "CIMS.pl (deletions)" "${CTK_DIR}/gt_cims_del.txt" \
        CIMS.pl -big -c "$CACHE_CIMS_DEL" -n "$CIMS_ITER" -FDR "$FDR" \
            "$COLLAPSED_BED" "$DEL_BED" "$CIMS_DEL"
    accumulate "$LAST_WALL" "$LAST_RSS"
    rm -rf "$CACHE_CIMS_DEL"
    # Post-filter by min mutations (col8) and min coverage (col7)
    if [[ -s "$CIMS_DEL" ]] && { [[ "$CIMS_M" -gt 1 ]] || [[ "$CIMS_K" -gt 1 ]]; }; then
        awk -v m="$CIMS_M" -v k="$CIMS_K" \
            'NR==1 || ($8 >= m && $7 >= k)' "$CIMS_DEL" > "${CIMS_DEL}.tmp" \
            && mv "${CIMS_DEL}.tmp" "$CIMS_DEL"
    fi
else
    printf "    %-35s %s\n" "CIMS.pl (deletions)" "SKIPPED (no deletions)"
fi

CTK_CIMS_DEL=$(count_data "$CIMS_DEL")

# ── CTK 6: CIMS.pl — substitutions ───────────────────────────────────────────
CIMS_SUB="${CTK_DIR}/CIMS_sub.txt"
if [[ -s "$SUB_BED" ]]; then
    CACHE_CIMS_SUB=$(mktemp -u "${TMPDIR:-/tmp}/cims_sub.XXXXXX")
    run_timed "CIMS.pl (substitutions)" "${CTK_DIR}/gt_cims_sub.txt" \
        CIMS.pl -big -c "$CACHE_CIMS_SUB" -n "$CIMS_ITER" -FDR "$FDR" \
            "$COLLAPSED_BED" "$SUB_BED" "$CIMS_SUB"
    accumulate "$LAST_WALL" "$LAST_RSS"
    rm -rf "$CACHE_CIMS_SUB"
    # Post-filter by min mutations (col8) and min coverage (col7)
    if [[ -s "$CIMS_SUB" ]] && { [[ "$CIMS_M" -gt 1 ]] || [[ "$CIMS_K" -gt 1 ]]; }; then
        awk -v m="$CIMS_M" -v k="$CIMS_K" \
            'NR==1 || ($8 >= m && $7 >= k)' "$CIMS_SUB" > "${CIMS_SUB}.tmp" \
            && mv "${CIMS_SUB}.tmp" "$CIMS_SUB"
    fi
else
    printf "    %-35s %s\n" "CIMS.pl (substitutions)" "SKIPPED (no substitutions)"
fi

CTK_CIMS_SUB=$(count_data "$CIMS_SUB")

# ── CTK 7: CITS.pl ───────────────────────────────────────────────────────────
CITS_TXT="${CTK_DIR}/CITS.txt"
CACHE_CITS=$(mktemp -u "${TMPDIR:-/tmp}/cits.XXXXXX")
run_timed "CITS.pl" "${CTK_DIR}/gt_cits.txt" \
    CITS.pl -big -c "$CACHE_CITS" -p "$FDR" --gap 25 \
        "$COLLAPSED_BED" "$DEL_BED" "$CITS_TXT"
accumulate "$LAST_WALL" "$LAST_RSS"
rm -rf "$CACHE_CITS"

CTK_CITS=$(count_data "$CITS_TXT")

CTK_TOTAL_WALL="$TOTAL_WALL"
CTK_PEAK_RSS="$PEAK_RSS"

echo "    ─────────────────────────────────────────────────────────"
printf "    %-35s %8ss  %8s MB\n" "TOTAL" "$CTK_TOTAL_WALL" "$CTK_PEAK_RSS"
echo ""

# ── Speedup ───────────────────────────────────────────────────────────────────
SPEEDUP=$(awk "BEGIN {printf \"%.1fx\", $CTK_TOTAL_WALL / $CLK_TOTAL_WALL}")

# ══════════════════════════════════════════════════════════════════════════════
# Report
# ══════════════════════════════════════════════════════════════════════════════
{
printf "\n════════════════════════════════════════════════════════════════\n"
printf "  CTK vs Clink Benchmark Report\n"
printf "════════════════════════════════════════════════════════════════\n"
printf "  Input BAM     : %s\n" "$INPUT_BAM"
printf "  Reads tested  : %s\n" "$(fmt_num $INPUT_READS)"
printf "  CIMS iter     : %s    min-mut: %s    min-cov: %s    FDR: %s\n" "$CIMS_ITER" "$CIMS_M" "$CIMS_K" "$FDR"
printf "  Genome FASTA  : %s\n" "${GENOME_FASTA:-none}"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s  %10s\n" "PERFORMANCE"             "Clink"      "CTK"
printf "  %-38s  %10s  %10s\n" "──────────────────────" "──────────" "──────────"
printf "  %-38s  %10s  %10s\n" "Total wall time (s)"     "$CLK_TOTAL_WALL"  "$CTK_TOTAL_WALL"
printf "  %-38s  %10s  %10s\n" "Peak RAM (MB)"            "$CLK_PEAK_RSS"    "$CTK_PEAK_RSS"
printf "  %-38s  %10s\n"       "Clink speedup"            "$SPEEDUP"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s  %10s\n" "INPUT / DEDUPLICATION"  "Clink"      "CTK"
printf "  %-38s  %10s  %10s\n" "──────────────────────" "──────────" "──────────"
printf "  %-38s  %10s  %10s\n" "Input reads"              "$(fmt_num $INPUT_READS)"  "$(fmt_num $INPUT_READS)"
printf "  %-38s  %10s  %10s\n" "Tags (parsed from SAM)"   "—"                        "$(fmt_num $CTK_TAGS)"
printf "  %-38s  %10s  %10s\n" "After dedup"              "$(fmt_num $CLK_DEDUP)"    "$(fmt_num $CTK_COLLAPSED)"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s  %10s\n" "RAW SIGNAL (min-cov≥5)"  "Clink"      "CTK"
printf "  NOTE: Clink=unique genomic positions; CTK=reads carrying event\n"
printf "  (same reads, different aggregation — not directly comparable)\n"
printf "  %-38s  %10s  %10s\n" "──────────────────────" "──────────" "──────────"
printf "  %-38s  %10s  %10s\n" "Positions/reads with deletion"     "$(fmt_num $CLK_RAW_DEL_POS)"  "$(fmt_num $CTK_DELS)"
printf "  %-38s  %10s  %10s\n" "Positions/reads with substitution" "$(fmt_num $CLK_RAW_SUB_POS)"  "$(fmt_num $CTK_SUBS)"
printf "  %-38s  %10s  %10s\n" "Positions/reads with truncation"   "$(fmt_num $CLK_RAW_TRUNC_POS)" "—"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s  %10s\n" "CIMS SIGNIFICANT SITES"  "Clink"      "CTK"
printf "  (FDR≤${FDR}, min-mut≥${CIMS_M}, min-cov≥${CIMS_K}; CTK: ${CIMS_ITER} permutations)\n"
printf "  %-38s  %10s  %10s\n" "──────────────────────" "──────────" "──────────"
printf "  %-38s  %10s  %10s\n" "Deletion sites"           "$(fmt_num $CLK_DELS)"         "$(fmt_num $CTK_CIMS_DEL)"
printf "  %-38s  %10s  %10s\n" "Substitution sites"       "$(fmt_num $CLK_SUBS_TOTAL)"   "$(fmt_num $CTK_CIMS_SUB)"
printf "  %-38s  %10s  %10s\n" "  of which T>C"           "$(fmt_num $CLK_TtoC)"         "—"
printf "────────────────────────────────────────────────────────────────\n"
printf "  %-38s  %10s  %10s\n" "CITS TRUNCATION SITES"   "Clink"      "CTK"
printf "  (FDR≤${FDR}, min-cov≥5, min-frac≥0.05; CTK requires read clusters)\n"
printf "  %-38s  %10s  %10s\n" "──────────────────────" "──────────" "──────────"
printf "  %-38s  %10s  %10s\n" "Truncation sites"         "$(fmt_num $CLK_TRUNC)"        "$(fmt_num $CTK_CITS)"
printf "════════════════════════════════════════════════════════════════\n"
} | tee "$REPORT"

echo ""
echo "  Report saved : $REPORT"
echo "  CTK outputs  : $CTK_DIR/"
echo "  Clink outputs: $CLK_DIR/"
echo ""
