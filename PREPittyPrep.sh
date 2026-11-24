#!/bin/bash

# PREPittyPrep.sh - Preprocessing Module for CLIPittyClip (v3.5.3)
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
source "${SCRIPT_DIR}/lib/wizard.sh"

# ─── Defaults ────────────────────────────────────────────────────────────────
THREADS=1
INPUT_FILE=""
INPUT_DIR=""
BARCODE_FILE=""
EXP_ID=""
UMI_LEN=0
BC_LEN=""
SPACER_LEN="0"
ADAPTER_3="L32"             # Same default as CLIPittyClip.sh (override with -a)
BC_FIRST="false"
TRIM3_LEN="0"
FASTP_MIN_QUAL="30"
DEDUP_MODE="true"
LOW_MEMORY="false"
GEO_MODE="false"
ECLIP_MODE=""               # "pe" or "se"; empty = off
PARCLIP_MODE="false"
PARCLIP_ADAPTERS=""
FILTER_REPEAT="false"
GENOME_INDEX=""
DEMUX_MISMATCHES="1"
SAMPLE_SIZE=0
KEEP_INTERMEDIATE="no"
WIZARD_MODE="false"

# ─── Usage ───────────────────────────────────────────────────────────────────
function show_usage {
    echo ""
    echo "Usage: $0 [-i <input.fastq.gz> | -d <dir>] [options]"
    echo ""
    echo "PREPittyPrep v3.5.3 - CLIP-seq Preprocessing Module"
    echo "Output: groomed, deduplicated, adapter-trimmed *_prepped.fastq.gz"
    echo ""
    echo "INPUT (choose one):"
    echo "  -i, --input-file <path>  Input FASTQ file (.fastq or .fastq.gz)"
    echo "  -d, --input-dir <dir>    Directory of FASTQ files (batch mode)"
    echo ""
    echo "PREPROCESSING:"
    echo "  -b, --barcodes <path>    Barcode file; triggers demux mode (needs -i)"
    echo "  -u, --umi-length <int>   UMI length (default: 0)"
    echo "  -a, --adapter <str>      3' adapter sequence, or a name: L32, L3, L3MOD, L7 (L19), SRA3"
    echo "                             (default: L32)"
    echo "  --bc-len <int>           Barcode length (auto-detected from -b)"
    echo "  --spacer-len <int>       Spacer length to trim after barcode (default: 0)"
    echo "  --bc-first               Barcode precedes UMI: layout [BC][UMI][spacer][READ]"
    echo "  --trim3 <int>            Bases to trim from the 3' end after adapter removal (default: 0)"
    echo "  --min-qual <int>         fastp average quality threshold (default: 30)"
    echo "  --no-dedup               Skip deduplication"
    echo "  --low-memory             Use the flat-RAM sort-based dedup engine"
    echo "  --eclip <pe|se>          eCLIP mode: 'pe' (post-eclipdemux R2, UMI in header) or"
    echo "                             'se' (raw seCLIP R1, UMI in sequence)"
    echo "  --parclip                PAR-CLIP mode: [UMI][READ][2nt spacer][6mer barcode][adapter]."
    echo "                             Requires -u."
    echo "  --parclip-adapters <fa>  Custom PAR-CLIP adapter FASTA (default: lib/parclip_adapters.fa)"
    echo "  --demux-mismatches <N>   Max barcode mismatches (default: 1)"
    echo "  --filter-repeat          Remove reads that map to the repeat index (needs -x)"
    echo "  -x, --index <path>       Genome index that holds the repeat index (only with --filter-repeat)"
    echo "  -s, --sample <int>       Test mode: process N reads only"
    echo ""
    echo "OUTPUT:"
    echo "  -o, --output <str>       Output folder name / path"
    echo "  -k, --keep               Keep intermediate files (including 0_DEMUX_FASTQ)"
    echo "  -t, --threads <int>      Threads (default: 1)"
    echo "  -w, --wizard           Launch interactive configuration wizard"
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
        -i|--input-file)      INPUT_FILE="$2"; shift 2 ;;
        -d|--input-dir)       INPUT_DIR="$2"; shift 2 ;;
        -b|--barcodes)        BARCODE_FILE="$2"; shift 2 ;;
        -o|--output)          EXP_ID="$2"; shift 2 ;;
        -t|--threads)         THREADS="$2"; shift 2 ;;
        -u|--umi-length)      UMI_LEN="$2"; shift 2 ;;
        -a|--adapter)         ADAPTER_3="$2"; shift 2 ;;
        --bc-len)             BC_LEN="$2"; shift 2 ;;
        --spacer-len)         SPACER_LEN="$2"; shift 2 ;;
        --bc-first)           BC_FIRST="true"; shift ;;
        --trim3)              TRIM3_LEN="$2"; shift 2 ;;
        --min-qual)           FASTP_MIN_QUAL="$2"; shift 2 ;;
        --no-dedup)           DEDUP_MODE="false"; shift ;;
        --low-memory)         LOW_MEMORY="true"; shift ;;
        --eclip)              ECLIP_MODE="$2"; shift 2 ;;
        --parclip)            PARCLIP_MODE="true"; shift ;;
        --parclip-adapters)   PARCLIP_ADAPTERS="$2"; shift 2 ;;
        --filter-repeat)      FILTER_REPEAT="true"; shift ;;
        --filter-ncrna)       echo "[WARNING] --filter-ncrna is deprecated; use --filter-repeat instead." >&2
                              FILTER_REPEAT="true"; shift ;;
        -x|--index)           GENOME_INDEX="$2"; shift 2 ;;
        --demux-mismatches)   DEMUX_MISMATCHES="$2"; shift 2 ;;
        -s|--sample)          SAMPLE_SIZE="$2"; shift 2 ;;
        -k|--keep)            KEEP_INTERMEDIATE="yes"; shift ;;
        --geo)                GEO_MODE="true"; shift ;;
        -w|--wizard|--advanced) WIZARD_MODE="true"; shift ;;
        -h|--help)            show_usage; exit 0 ;;
        *)                    log_error "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# ─── Wizard (if requested) ───────────────────────────────────────────────────
if [[ "$WIZARD_MODE" == "true" ]]; then
    run_wizard_prepittyprep || exit 1

    # Input mode
    if [[ "$WIZ_MODE" == "single" ]]; then
        INPUT_FILE="$WIZ_INPUT_FILE"
    else
        INPUT_DIR="$WIZ_INPUT_DIR"
    fi
    [[ -n "$WIZ_BARCODE_FILE" ]] && BARCODE_FILE="$WIZ_BARCODE_FILE"
    [[ "$WIZ_GEO_MODE" == "true" ]] && GEO_MODE="true"

    THREADS="$WIZ_THREADS"
    UMI_LEN="$WIZ_UMI_LEN"
    ADAPTER_3="$WIZ_ADAPTER"
    [[ -n "$WIZ_BC_LEN" ]] && BC_LEN="$WIZ_BC_LEN"
    SPACER_LEN="$WIZ_SPACER_LEN"
    [[ "$WIZ_BC_FIRST" == "true" ]] && BC_FIRST="true"
    TRIM3_LEN="$WIZ_TRIM3_LEN"
    [[ "$WIZ_NO_DEDUP" == "true" ]] && DEDUP_MODE="false"
    [[ -n "$WIZ_ECLIP_MODE" ]] && ECLIP_MODE="$WIZ_ECLIP_MODE"
    [[ "$WIZ_PARCLIP_MODE" == "true" ]] && PARCLIP_MODE="true"
    DEMUX_MISMATCHES="$WIZ_DEMUX_MISMATCHES"
    [[ "$WIZ_FILTER_REPEAT" == "true" ]] && FILTER_REPEAT="true"
    [[ -n "$WIZ_GENOME_INDEX" ]] && GENOME_INDEX="$WIZ_GENOME_INDEX"
    [[ -n "$WIZ_OUTPUT" ]] && EXP_ID="$WIZ_OUTPUT"
    [[ "$WIZ_KEEP" == "true" ]] && KEEP_INTERMEDIATE="yes"
    SAMPLE_SIZE="$WIZ_SAMPLE_SIZE"
fi

# Convert an adapter name to its sequence
ADAPTER_3=$(resolve_adapter "$ADAPTER_3") || exit 1

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
if [[ "$FILTER_REPEAT" == "true" && -z "$GENOME_INDEX" ]]; then
    echo "[ERROR] --filter-repeat requires a genome index (-x) that holds the repeat index." >&2; exit 1
fi
if [[ -n "$ECLIP_MODE" && "$ECLIP_MODE" != "pe" && "$ECLIP_MODE" != "se" ]]; then
    echo "[ERROR] --eclip requires 'pe' or 'se' (got: '$ECLIP_MODE')." >&2; exit 1
fi
if [[ -n "$ECLIP_MODE" && "$PARCLIP_MODE" == "true" ]]; then
    echo "[ERROR] Use --eclip or --parclip, not both." >&2; exit 1
fi
if [[ "$PARCLIP_MODE" == "true" && "${UMI_LEN:-0}" -le 0 ]]; then
    echo "[ERROR] --parclip requires -u <umi_length>. Specify the UMI length (for example -u 4)." >&2; exit 1
fi
if [[ -n "$PARCLIP_ADAPTERS" && ! -f "$PARCLIP_ADAPTERS" ]]; then
    echo "[ERROR] PAR-CLIP adapter FASTA not found: $PARCLIP_ADAPTERS" >&2; exit 1
fi

# ─── Dependency checks ───────────────────────────────────────────────────────
check_dependency fastp
if [[ -n "$BARCODE_FILE" ]]; then check_dependency cutadapt; fi
if [[ "$FILTER_REPEAT" == "true" ]]; then check_dependency bowtie2; fi

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

# ─── Banner ──────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════════════════"
echo "  PREPittyPrep v3.5.3"
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
if   [[ -n "$ECLIP_MODE" ]];           then echo "  Protocol: eCLIP ${ECLIP_MODE}"
elif [[ "$PARCLIP_MODE" == "true" ]]; then echo "  Protocol: PAR-CLIP"
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
if [[ -n "$EXP_ID" ]]; then
    if [[ "$EXP_ID" == /* ]]; then
        OUTPUT_ROOT="${EXP_ID}"
    else
        OUTPUT_ROOT="${WORK_PARENT}/${EXP_ID}"
    fi
else
    OUTPUT_ROOT="${WORK_PARENT}/${BASE}_prepped"
fi

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
# Same steps and order as the CLIPittyClip.sh per-sample preprocessing.
# Args: $1=input_fastq $2=sample_name $3=work_dir $4=skip_dedup (true after pool dedup)
preprocess_one() {
    local in_file="$1"
    local sname="$2"
    local wdir="$3"
    local skip_dedup="${4:-false}"
    local work="$in_file"
    local prefix="${wdir}/${sname}"

    # 1. Dedup (the eCLIP preprocessors deduplicate inside their own steps)
    if [[ "$DEDUP_MODE" == "true" && -z "$ECLIP_MODE" && "$skip_dedup" != "true" ]]; then
        local dd_out="${prefix}_dedup.fastq"
        if run_dedup "$in_file" "$dd_out"; then
            work="$dd_out"
        else
            log_warning "Dedup failed for $sname. Continuing with original."
        fi
    fi

    # 2. Adapter trimming and UMI extraction, by protocol
    if   [[ "$ECLIP_MODE" == "pe" ]]; then
        run_eclip_pe_preprocessing "$work" "$prefix" "$THREADS" "$SAMPLE_SIZE" "$UMI_LEN"
    elif [[ "$ECLIP_MODE" == "se" ]]; then
        run_eclip_se_preprocessing "$work" "$prefix" "$THREADS" "$SAMPLE_SIZE"
    elif [[ "$PARCLIP_MODE" == "true" ]]; then
        run_parclip_preprocessing "$work" "$prefix" "$UMI_LEN" "$THREADS" "$SAMPLE_SIZE" "$PARCLIP_ADAPTERS"
    else
        run_fastp "$work" "$prefix" "$UMI_LEN" "$ADAPTER_3" "$THREADS" "$SAMPLE_SIZE" \
            "$BC_LEN" "$SPACER_LEN" "$BC_FIRST" "$FASTP_MIN_QUAL" "$TRIM3_LEN"
    fi
    local cleaned="${prefix}_cleaned.fastq"

    # 3. Optional repeat-element filter (after trimming, as in CLIPittyClip.sh)
    if [[ "$FILTER_REPEAT" == "true" && -s "$cleaned" ]]; then
        local repeat_idx
        repeat_idx=$(check_repeat_index "$GENOME_INDEX")
        if [[ -n "$repeat_idx" ]]; then
            local repeat_dir="${wdir}/Repeat_Mapping"
            local repeat_out="${prefix}_repeat_filtered.fastq"
            run_repeat_filter "$cleaned" "$repeat_out" "$repeat_dir" "$repeat_idx" "$THREADS" "$sname"
            run_repeat_quantify "${repeat_dir}/${sname}_repeat.bam" \
                "${repeat_dir}/${sname}_repeat_stats.txt" "$repeat_dir" "$sname"
            if [[ -s "$repeat_out" ]]; then
                mv "$repeat_out" "$cleaned"
            fi
        else
            log_warning "Repeat index not found in $GENOME_INDEX or $GENOME_INDEX/Repeat. Skipping repeat filter for $sname."
        fi
    fi

    # 4. Gzip cleaned output → _prepped.fastq.gz
    local prepped="${prefix}_prepped.fastq.gz"
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
    if [[ -d "${wdir}/Repeat_Mapping" ]]; then
        mkdir -p "$OUTPUT_ROOT/REPORTS/Repeat_Mapping"
        cp "${wdir}/Repeat_Mapping"/*.{txt,tsv} "$OUTPUT_ROOT/REPORTS/Repeat_Mapping/" 2>/dev/null
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# DIRECTORY MODE
# ─────────────────────────────────────────────────────────────────────────────
if [[ -n "$INPUT_DIR" ]]; then
    # Plain array assignment instead of mapfile, for macOS bash 3.2
    SAMPLE_FILES=($(ls "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz \
                       "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq 2>/dev/null | sort))
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
        DEDUP_START=$(date +%s)
        if run_dedup "$INPUT_FILE" "$POOL_DEDUP"; then
            DEDUP_END=$(date +%s)
            DEDUP_DURATION=$((DEDUP_END - DEDUP_START))
            DEDUP_M=$((DEDUP_DURATION/60))
            DEDUP_S=$((DEDUP_DURATION%60))
            echo "  > Deduplication complete. (Total duration: ${DEDUP_M}min ${DEDUP_S}secs)"
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
    run_demultiplexing "$WORK_INPUT" "$BARCODE_FILE" "$SAMPLE_SIZE" "demux_fastq"
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
        # Dedup already ran on the pooled file, so skip it here
        if preprocess_one "$f" "$sname" "$sdir" "true"; then
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
