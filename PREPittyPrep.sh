#!/bin/bash

# PREPittyPrep.sh - Preprocessing Module for CLIPittyClip (v1.0)
#
# Runs the full CLIP-seq preprocessing stack and stops before alignment:
#   dedup → [demux] → fastp → *_prepped.fastq.gz (ready to map)
#
# GEO mode (--geo):
#   Raw barcode splitting only — reads written exactly as received.
#   No dedup, no fastp. For GEO deposit of pooled libraries.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/dedup.sh"
source "${SCRIPT_DIR}/lib/modules.sh"

# ─── Defaults ────────────────────────────────────────────────────────────────
THREADS=1
INPUT_FILE=""
INPUT_DIR=""
BARCODE_FILE=""
EXP_ID=""
UMI_LEN=0
BC_LEN=""
SPACER_LEN="0"
ADAPTER_3="AGATCGGAAGAGC"   # Illumina universal adapter (override with -a)
DEDUP_MODE="true"
GEO_MODE="false"
ECLIP_MODE="false"
FILTER_NCRNA="false"
GENOME_INDEX=""
DEMUX_MISMATCHES="1"
SAMPLE_SIZE=0
KEEP_INTERMEDIATE="no"

# ─── Usage ───────────────────────────────────────────────────────────────────
function show_usage {
    echo ""
    echo "Usage: $0 [-i <input.fastq.gz> | -d <dir>] [options]"
    echo ""
    echo "PREPittyPrep v1.0 - CLIP-seq Preprocessing Module"
    echo "Output: groomed, deduplicated, adapter-trimmed *_prepped.fastq.gz"
    echo ""
    echo "INPUT (choose one):"
    echo "  -i <path>              Input FASTQ file (.fastq or .fastq.gz)"
    echo "  -d <dir>               Directory of FASTQ files (batch mode)"
    echo ""
    echo "PREPROCESSING:"
    echo "  -b <path>              Barcode file; triggers demux mode"
    echo "  -u <int>               UMI length (default: 0)"
    echo "  -a <str>               3' adapter sequence"
    echo "  --bc-len <int>         Barcode length (auto-detected from -b)"
    echo "  --spacer-len <int>     Spacer length to trim after barcode (default: 0)"
    echo "  --no-dedup             Skip deduplication"
    echo "  --eclip                eCLIP mode (UMI in read header)"
    echo "  --demux-mismatches <N> Max barcode mismatches (default: 1)"
    echo "  --filter-ncrna         Enable ncRNA pre-filtering (requires -x)"
    echo "  -x <path>              Genome index (required only with --filter-ncrna)"
    echo "  -s <int>               Test mode: process N reads only"
    echo ""
    echo "OUTPUT:"
    echo "  -o <str>               Output folder name / path"
    echo "  -k                     Keep intermediate files (including 0_DEMUX_FASTQ)"
    echo "  -t <int>               Threads (default: 1)"
    echo "  -h, --help             Show this help"
    echo ""
    echo "MODES:"
    echo "  --geo                  GEO mode: raw demux only (requires -i and -b)"
    echo "                         Reads split by barcode, not modified."
    echo ""
    echo "EXAMPLES:"
    echo "  # Single file preprocessing"
    echo "  $0 -i reads.fastq.gz -u 7 -t 8"
    echo ""
    echo "  # Demux + preprocess pooled library"
    echo "  $0 -i pool.fastq.gz -b barcodes.txt -u 7 -t 8"
    echo ""
    echo "  # Batch preprocess a directory of FASTQs"
    echo "  $0 -d /path/to/samples/ -u 7 -t 8"
    echo ""
    echo "  # GEO deposit: raw demux only"
    echo "  $0 -i pool.fastq.gz -b barcodes.txt -u 7 --geo -o my_GEO"
    echo ""
}

if [[ $# -eq 0 ]]; then show_usage; exit 1; fi

# ─── Argument Parsing ────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        -i)                   INPUT_FILE="$2"; shift 2 ;;
        -d)                   INPUT_DIR="$2"; shift 2 ;;
        -b)                   BARCODE_FILE="$2"; shift 2 ;;
        -o)                   EXP_ID="$2"; shift 2 ;;
        -t)                   THREADS="$2"; shift 2 ;;
        -u)                   UMI_LEN="$2"; shift 2 ;;
        -a)                   ADAPTER_3="$2"; shift 2 ;;
        --bc-len)             BC_LEN="$2"; shift 2 ;;
        --spacer-len)         SPACER_LEN="$2"; shift 2 ;;
        --no-dedup)           DEDUP_MODE="false"; shift ;;
        --eclip)              ECLIP_MODE="true"; shift ;;
        --filter-ncrna)       FILTER_NCRNA="true"; shift ;;
        -x)                   GENOME_INDEX="$2"; shift 2 ;;
        --demux-mismatches)   DEMUX_MISMATCHES="$2"; shift 2 ;;
        -s)                   SAMPLE_SIZE="$2"; shift 2 ;;
        -k)                   KEEP_INTERMEDIATE="yes"; shift ;;
        --geo)                GEO_MODE="true"; shift ;;
        -h|--help)            show_usage; exit 0 ;;
        *)                    log_error "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# ─── Validation ──────────────────────────────────────────────────────────────
if [[ -n "$INPUT_DIR" && -n "$INPUT_FILE" ]]; then
    log_error "Cannot use both -i and -d. Choose one."; exit 1
fi
if [[ -z "$INPUT_FILE" && -z "$INPUT_DIR" ]]; then
    log_error "Missing required input (-i or -d)."; show_usage; exit 1
fi
if [[ "$GEO_MODE" == "true" && -z "$BARCODE_FILE" ]]; then
    log_error "--geo requires a barcode file (-b)."; exit 1
fi
if [[ "$GEO_MODE" == "true" && -z "$INPUT_FILE" ]]; then
    log_error "--geo requires a single pooled input file (-i)."; exit 1
fi
if [[ -n "$BARCODE_FILE" && "$GEO_MODE" == "false" && -z "$INPUT_FILE" ]]; then
    log_error "Demux mode (-b) requires a single pooled input file (-i), not -d."; exit 1
fi
if [[ "$FILTER_NCRNA" == "true" && -z "$GENOME_INDEX" ]]; then
    log_error "--filter-ncrna requires a genome index (-x)."; exit 1
fi

# ─── Dependency checks ───────────────────────────────────────────────────────
check_dependency fastp
if [[ -n "$BARCODE_FILE" ]]; then check_dependency cutadapt; fi
if [[ "$FILTER_NCRNA" == "true" ]]; then check_dependency bowtie2; fi

# ─── Barcode length auto-detection ───────────────────────────────────────────
if [[ -n "$BARCODE_FILE" ]]; then
    [[ -f "$BARCODE_FILE" ]] || { log_error "Barcode file not found: $BARCODE_FILE"; exit 1; }
    bc_file_len=$(awk '!/^#/{print length($2); exit}' "$BARCODE_FILE")
    BC_LEN="${BC_LEN:-$bc_file_len}"
fi
BC_LEN="${BC_LEN:-0}"

# ─── Resolve absolute paths ──────────────────────────────────────────────────
if [[ -n "$INPUT_FILE" ]]; then
    check_file "$INPUT_FILE" || exit 1
    INPUT_FILE="$(cd "$(dirname "$INPUT_FILE")" && pwd)/$(basename "$INPUT_FILE")"
fi
if [[ -n "$INPUT_DIR" ]]; then
    [[ -d "$INPUT_DIR" ]] || { log_error "Input directory not found: $INPUT_DIR"; exit 1; }
    INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
fi
if [[ -n "$GENOME_INDEX" ]]; then
    GENOME_INDEX="$(cd "$GENOME_INDEX" && pwd)"
fi
if [[ -n "$BARCODE_FILE" ]]; then
    BARCODE_FILE="$(cd "$(dirname "$BARCODE_FILE")" && pwd)/$(basename "$BARCODE_FILE")"
fi

# ─── Determine base name and output root ─────────────────────────────────────
if [[ -n "$INPUT_FILE" ]]; then
    WORK_PARENT="$(cd "$(dirname "$INPUT_FILE")" && pwd)"
    BASE=$(basename "$INPUT_FILE"); BASE="${BASE%.fastq.gz}"; BASE="${BASE%.fq.gz}"; BASE="${BASE%.fastq}"; BASE="${BASE%.fq}"
else
    WORK_PARENT="$(cd "$(dirname "$INPUT_DIR")" && pwd)"
    BASE=$(basename "$INPUT_DIR")
fi
[[ -n "$EXP_ID" ]] && BASE="$EXP_ID"

# ─── Banner ──────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════════════════"
echo "  PREPittyPrep v1.0"
echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
if [[ -n "$INPUT_DIR" ]]; then
    echo "  Input:   $INPUT_DIR  [directory mode]"
else
    echo "  Input:   $(basename "$INPUT_FILE")"
fi
if [[ "$GEO_MODE" == "true" ]]; then
    echo "  Mode:    GEO (raw demux — no read modification)"
else
    echo "  Mode:    PREPROCESS  [dedup → fastp → _prepped.fastq.gz]"
    [[ -n "$BARCODE_FILE" ]] && echo "           + DEMUX  [cutadapt → per-sample]"
fi
echo "  UMI:     ${UMI_LEN}bp | BC: ${BC_LEN}bp | Spacer: ${SPACER_LEN}bp | Threads: $THREADS"
echo "════════════════════════════════════════════════════════════════════════════════"

# ─────────────────────────────────────────────────────────────────────────────
# GEO MODE — demux only, exit immediately after
# ─────────────────────────────────────────────────────────────────────────────
if [[ "$GEO_MODE" == "true" ]]; then
    GEO_OUT="${WORK_PARENT}/${BASE}_GEO"
    echo ""
    echo "[GEO DEMUX]"
    echo "  Barcodes: $(basename "$BARCODE_FILE")"
    echo "  UMI offset: ${UMI_LEN}bp"
    echo "  Output: $GEO_OUT/"
    run_geo_demux "$INPUT_FILE" "$BARCODE_FILE" "$GEO_OUT" "$UMI_LEN" "$DEMUX_MISMATCHES"
    echo ""
    echo "[DONE] GEO files saved to: $GEO_OUT/"
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# NORMAL MODE — output directory setup
# ─────────────────────────────────────────────────────────────────────────────
OUTPUT_ROOT="${WORK_PARENT}/${BASE}_prepped"
mkdir -p "$OUTPUT_ROOT/PREPPED_FASTQ" \
         "$OUTPUT_ROOT/REPORTS/FASTP_REPORT"

# Log file inside output dir
LOG_FILE="${OUTPUT_ROOT}/REPORTS/detailed_output.log"
> "$LOG_FILE"
log_info "PREPittyPrep started: $(date '+%Y-%m-%d %H:%M:%S')"

# Scratch dir for intermediate files
SCRATCH="${OUTPUT_ROOT}/.prep_scratch"
mkdir -p "$SCRATCH"

# ─────────────────────────────────────────────────────────────────────────────
# Helper: preprocess one sample → _prepped.fastq.gz in $out_dir
# Args: $1=input_fastq $2=sample_name $3=work_dir
# ─────────────────────────────────────────────────────────────────────────────
preprocess_one() {
    local in_file="$1"
    local sname="$2"
    local wdir="$3"
    local work="$in_file"

    # 1. Dedup
    if [[ "$DEDUP_MODE" == "true" ]]; then
        local dd_out="${wdir}/${sname}_dedup.fastq"
        if run_dedup "$in_file" "$dd_out"; then
            work="$dd_out"
        else
            log_warning "Dedup failed for $sname. Continuing with original."
        fi
    fi

    # 2. Optional ncRNA filter
    if [[ "$FILTER_NCRNA" == "true" && -n "$GENOME_INDEX" ]]; then
        local ncrna_dir="${wdir}/ncRNA_Mapping"
        local ncrna_out="${wdir}/${sname}_ncrna_filtered.fastq"
        local ncrna_idx
        ncrna_idx=$(check_ncrna_index "$GENOME_INDEX")
        if [[ -n "$ncrna_idx" ]]; then
            run_ncrna_filter "$work" "$ncrna_out" "$ncrna_dir" "$ncrna_idx" "$THREADS" "$sname"
            work="$ncrna_out"
        else
            log_warning "ncRNA index not found. Skipping ncRNA filter for $sname."
        fi
    fi

    # 3. fastp (adapter trim, UMI extraction, front trim)
    run_fastp "$work" "${wdir}/${sname}" "$UMI_LEN" "$ADAPTER_3" "$THREADS" "$SAMPLE_SIZE" "$ECLIP_MODE" "$BC_LEN" "$SPACER_LEN"

    # 4. Gzip cleaned output → _prepped.fastq.gz
    local cleaned="${wdir}/${sname}_cleaned.fastq"
    local prepped="${wdir}/${sname}_prepped.fastq.gz"
    if [[ -s "$cleaned" ]]; then
        gzip -c "$cleaned" > "$prepped"
        rm -f "$cleaned"
        log_info "Prepared: $(basename "$prepped")"
        return 0
    else
        log_error "fastp produced empty output for $sname"
        return 1
    fi
}

# Helper: collect outputs from a work dir into the output root
collect_one() {
    local sname="$1"
    local wdir="$2"
    [[ -f "${wdir}/${sname}_prepped.fastq.gz" ]] && \
        cp "${wdir}/${sname}_prepped.fastq.gz" "$OUTPUT_ROOT/PREPPED_FASTQ/"
    [[ -f "${wdir}/${sname}_fastp.html" ]] && \
        cp "${wdir}/${sname}_fastp.html" "$OUTPUT_ROOT/REPORTS/FASTP_REPORT/"
    [[ -f "${wdir}/${sname}_fastp.json" ]] && \
        cp "${wdir}/${sname}_fastp.json" "$OUTPUT_ROOT/REPORTS/FASTP_REPORT/"
}

# ─────────────────────────────────────────────────────────────────────────────
# DIRECTORY MODE
# ─────────────────────────────────────────────────────────────────────────────
if [[ -n "$INPUT_DIR" ]]; then
    mapfile -t SAMPLE_FILES < <(ls "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz \
                                    "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq 2>/dev/null | sort)
    total=${#SAMPLE_FILES[@]}
    echo ""
    echo "[DIRECTORY MODE] Found $total input files"
    echo ""
    echo "[PREPROCESSING]"
    current=0
    for sample in "${SAMPLE_FILES[@]}"; do
        [[ -f "$sample" ]] || continue
        ((current++))
        sname=$(basename "$sample")
        sname="${sname%.fastq.gz}"; sname="${sname%.fq.gz}"
        sname="${sname%.fastq}";    sname="${sname%.fq}"
        sdir="${SCRATCH}/${sname}"
        mkdir -p "$sdir"
        printf "  %3d/%-3d  %-28s" "$current" "$total" "$sname"
        preprocess_one "$sample" "$sname" "$sdir"
        if [[ $? -eq 0 ]]; then
            printf "  ✓\n"
            collect_one "$sname" "$sdir"
        else
            printf "  ✗ FAILED\n"
        fi
        [[ "$KEEP_INTERMEDIATE" != "yes" ]] && rm -rf "$sdir"
    done

# ─────────────────────────────────────────────────────────────────────────────
# POOLED + DEMUX MODE
# ─────────────────────────────────────────────────────────────────────────────
elif [[ -n "$BARCODE_FILE" ]]; then
    # Step 1: Dedup the pool
    echo ""
    echo "[DEDUPLICATING]"
    POOL_DEDUP="${SCRATCH}/pooled_dedup.fastq"
    WORK_INPUT="$INPUT_FILE"
    if [[ "$DEDUP_MODE" == "true" ]]; then
        if run_dedup "$INPUT_FILE" "$POOL_DEDUP"; then
            echo "  > Deduplication complete."
            WORK_INPUT="$POOL_DEDUP"
        else
            log_warning "Deduplication failed. Using original file."
        fi
    else
        echo "  > Skipped (--no-dedup)."
    fi

    # Step 2: Demux (run_demultiplexing uses CWD, so cd into SCRATCH)
    echo ""
    echo "[DEMULTIPLEXING]"
    echo "  > Barcodes: $(basename "$BARCODE_FILE")"
    pushd "$SCRATCH" > /dev/null
    run_demultiplexing "$WORK_INPUT" "$BARCODE_FILE" "$SAMPLE_SIZE" "false"
    popd > /dev/null
    [[ -f "$POOL_DEDUP" ]] && rm -f "$POOL_DEDUP"

    # Demux summary
    echo ""
    echo "[DEMUX SUMMARY]"
    printf "  %-25s %-12s %s\n" "Sample" "Reads" "% of Total"
    echo "  -----------------------------------------------"
    total_reads=0
    for f in "${SCRATCH}/demux_fastq"/*.fastq; do
        [[ -f "$f" ]] || continue
        lines=$(wc -l < "$f"); total_reads=$((total_reads + lines / 4))
    done
    for f in "${SCRATCH}/demux_fastq"/*.fastq; do
        [[ -f "$f" ]] || continue
        sname=$(basename "$f" .fastq)
        lines=$(wc -l < "$f"); count=$((lines / 4))
        pct=$(awk "BEGIN {printf \"%.1f\", ($total_reads>0)?($count/$total_reads)*100:0}")
        printf "  %-25s %-12s %s%%\n" "$sname" "$count" "$pct"
    done
    echo "  -----------------------------------------------"
    echo "  Total demuxed: $total_reads reads"

    # Step 3: Per-sample fastp
    echo ""
    echo "[PREPROCESSING (per sample)]"
    total_samples=$(ls "${SCRATCH}/demux_fastq"/*.fastq 2>/dev/null | grep -vc "unknown.fastq")
    current=0
    for f in "${SCRATCH}/demux_fastq"/*.fastq; do
        [[ -f "$f" ]] || continue
        sname=$(basename "$f" .fastq)
        [[ "$sname" == "unknown" ]] && continue
        ((current++))
        sdir="${SCRATCH}/${sname}"
        mkdir -p "$sdir"
        printf "  %3d/%-3d  %-28s" "$current" "$total_samples" "$sname"
        # For demux samples, dedup already done on pool; skip dedup here
        # Run fastp directly (dedup was already done on the pooled file)
        run_fastp "$f" "${sdir}/${sname}" "$UMI_LEN" "$ADAPTER_3" "$THREADS" "$SAMPLE_SIZE" "$ECLIP_MODE" "$BC_LEN" "$SPACER_LEN"
        cleaned="${sdir}/${sname}_cleaned.fastq"
        if [[ -s "$cleaned" ]]; then
            gzip -c "$cleaned" > "${sdir}/${sname}_prepped.fastq.gz"
            rm -f "$cleaned"
            printf "  ✓\n"
            collect_one "$sname" "$sdir"
        else
            printf "  ✗ FAILED\n"
        fi
        [[ "$KEEP_INTERMEDIATE" != "yes" ]] && rm -rf "$sdir"
    done

    # Handle 0_DEMUX_FASTQ
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then
        mkdir -p "$OUTPUT_ROOT/0_DEMUX_FASTQ"
        mv "${SCRATCH}/demux_fastq"/*.fastq "$OUTPUT_ROOT/0_DEMUX_FASTQ/" 2>/dev/null
        log_info "Demux FASTQs kept in $OUTPUT_ROOT/0_DEMUX_FASTQ/"
    else
        rm -rf "${SCRATCH}/demux_fastq"
    fi

# ─────────────────────────────────────────────────────────────────────────────
# SINGLE FILE MODE
# ─────────────────────────────────────────────────────────────────────────────
else
    echo ""
    echo "[PREPROCESSING]"
    sdir="${SCRATCH}/${BASE}"
    mkdir -p "$sdir"
    printf "  1/1    %-28s" "$BASE"
    preprocess_one "$INPUT_FILE" "$BASE" "$sdir"
    if [[ $? -eq 0 ]]; then
        printf "  ✓\n"
        collect_one "$BASE" "$sdir"
    else
        printf "  ✗ FAILED\n"
    fi
    [[ "$KEEP_INTERMEDIATE" != "yes" ]] && rm -rf "$sdir"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Cleanup scratch and summary
# ─────────────────────────────────────────────────────────────────────────────
[[ "$KEEP_INTERMEDIATE" != "yes" ]] && rm -rf "$SCRATCH"

echo ""
echo "[OUTPUT]"
echo "  All results saved to: $OUTPUT_ROOT/"
echo "  ├── PREPPED_FASTQ/   ← ready-to-map gzipped FASTQs"
if [[ -n "$BARCODE_FILE" && "$KEEP_INTERMEDIATE" == "yes" ]]; then
    echo "  ├── 0_DEMUX_FASTQ/  ← raw demux FASTQs (kept with -k)"
fi
echo "  └── REPORTS/"
echo "        ├── FASTP_REPORT/"
echo "        └── detailed_output.log"
echo ""

n_prepped=$(ls "$OUTPUT_ROOT/PREPPED_FASTQ"/*.fastq.gz 2>/dev/null | wc -l)
echo "[SUCCESS] PREPittyPrep finished. $n_prepped prepped file(s) ready."
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
exit 0
