#!/bin/bash

# MAPittyMap.sh - Mapping-only module for CLIPittyClip (v3.5)
# Uses the unified lib/modules.sh for consistent behavior

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/modules.sh"
# Source wizard but don't run it immediately
source "${SCRIPT_DIR}/lib/wizard.sh"

# Default Values
THREADS=1
GENOME_INDEX=""
INPUT_FILE=""
EXP_ID=""
ALIGNER="star"
MISMATCH_MAX=2
UMI_LEN=0
GENOME_FASTA=""       # Path to reference FASTA (optional; strongly recommended for CIMS)
WIZARD_MODE="false"
FILTER_REPEAT="false"  # Repeat element filtering is OFF by default; enable with --filter-repeat
FILTER_CHR="true"     # Canonical chromosome filter ON by default; disable with --no-chr-filter
BAM_ONLY="false"      # If true, skip post-mapping (no BED, no bedgraph, no CIMS/CITS, no Clink)
DEDUP_MODE="true"     # PCR dedup (CTK path); always on

# Analysis tracks (opt-in)
RUN_CIMS="false"
RUN_CITS="false"
RUN_CLINK="false"

# CTK CIMS/CITS knobs
CIMS_ITERATIONS=5
CIMS_FDR="0.05"
CITS_PVALUE="0.05"
CITS_GAP=25

# Clink knobs
CLINK_UMI_LEN="-1"    # -1 = auto-detect from read names
CLINK_FDR="0.05"
CLINK_MIN_COV=5
CLINK_MIN_FRAC="0.05"  # internal; not exposed
CLINK_MULTI_MAP="false"

# Motif / flanked BED
RUN_MOTIF="yes"
MOTIF_FLANK=10

function show_usage {
    echo ""
    echo "Usage: $0 -i <input.fastq.gz> -x <index_dir> [options]"
    echo ""
    echo "MAPittyMap v3.5 - Standalone Mapping + Crosslink-Site Analysis Module"
    echo ""
    echo "REQUIRED:"
    echo "  -i <path>      Input FASTQ/FASTA file (gzipped supported)"
    echo "  -x <path>      Path to genome index directory"
    echo ""
    echo "ALIGNMENT OPTIONS:"
    echo "  --aligner <str>       Aligner: 'star' (default) or 'bowtie2'"
    echo "  --genome-fasta <path> Reference FASTA (enables samtools calmd; recommended for CIMS)"
    echo "  -o <str>              Output ID (default: derived from filename)"
    echo "  -t <int>              Number of threads (default: 1)"
    echo "  -m <int>              Max absolute mismatches (default: 2)"
    echo "  -u, --umi-length <int>  UMI length in bp (default: 0)"
    echo "  --filter-repeat       Enable repeat element pre-filtering (off by default)"
    echo "  --no-chr-filter       Skip canonical-chromosome filter (kept ON by default)"
    echo ""
    echo "POST-ALIGNMENT OUTPUTS:"
    echo "  --bam-only            Stop after mapping (no BED / bedgraph / CIMS / CITS / Clink)"
    echo ""
    echo "ANALYSIS TRACKS (opt-in):"
    echo "  --run-cims            CTK CIMS (mutation-site) analysis"
    echo "  --run-cits            CTK CITS (truncation-site) analysis"
    echo "  --run-cims-cits       Both CTK CIMS and CITS"
    echo "  --cims-iter <int>     CIMS permutation iterations (default: 5)"
    echo "  --cims-fdr <float>    CIMS FDR threshold (default: 0.05)"
    echo "  --cits-pval <float>   CITS p-value threshold (default: 0.05)"
    echo "  --cits-gap <int>      CITS clustering gap (default: 25, -1=off)"
    echo "  --run-clink           Clink pipeline (umi_tools dedup + Python pileup → CITS/CIMS)"
    echo "  --clink-umi-len <int>  UMI length for umi_tools (default: auto-detect)"
    echo "  --clink-fdr <float>   Clink FDR threshold (default: 0.05)"
    echo "  --clink-min-cov <int> Clink minimum coverage (default: 5)"
    echo "  --clink-multi-map     EM-rescue of multi-mapped reads before Clink dedup"
    echo "  --no-motif            Skip flanked BED generation for CTK CIMS/CITS"
    echo "  -f, --flank <int>     Flanked BED nucleotides (default: 10)"
    echo ""
    echo "WIZARD:"
    echo "  --wizard, --advanced  Launch interactive configuration wizard"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "OUTPUT TREE ({BASENAME}_mapping/):"
    echo "  1_BAM/                Sorted, indexed BAM"
    echo "  2_BED/                PCR-collapsed BED (feeds PEAKittyPeak)"
    echo "  3_COVERAGE/           Strand-specific normalized bedgraphs"
    echo "  4_CTK_Analysis/       (--run-cims/--run-cits) CIMS / CITS BED files"
    echo "  5_Clink_Analysis/     (--run-clink) Clink truncations + mutations + crosslinks"
    echo "  REPORTS/              Alignment logs"
    echo ""
    echo "EXAMPLES:"
    echo "  # Map only (legacy behavior)"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x star_index/ --bam-only"
    echo ""
    echo "  # Map + collapsed BED + bedgraph (default)"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x star_index/ -t 8"
    echo ""
    echo "  # Add CTK CIMS+CITS"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x star_index/ --run-cims-cits --genome-fasta ref.fa -t 8"
    echo ""
    echo "  # Add Clink track"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x star_index/ --run-clink -u 7 -t 8"
    echo ""
}

if [[ $# -eq 0 ]]; then
    show_usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -i) INPUT_FILE="$2"; shift 2 ;;
        -x) GENOME_INDEX="$2"; shift 2 ;;
        -o) EXP_ID="$2"; shift 2 ;;
        --aligner) ALIGNER=$(echo "$2" | tr '[:upper:]' '[:lower:]'); shift 2 ;;
        --genome-fasta) GENOME_FASTA="$2"; shift 2 ;;
        --wizard|--advanced) WIZARD_MODE="true"; shift 1 ;;
        --filter-repeat) FILTER_REPEAT="true"; shift ;;
        --no-chr-filter) FILTER_CHR="false"; shift ;;
        --bam-only) BAM_ONLY="true"; shift ;;
        -u|--umi-length) UMI_LEN="$2"; shift 2 ;;
        --run-cims) RUN_CIMS="true"; shift ;;
        --run-cits) RUN_CITS="true"; shift ;;
        --run-cims-cits) RUN_CIMS="true"; RUN_CITS="true"; shift ;;
        --cims-iter) CIMS_ITERATIONS="$2"; shift 2 ;;
        --cims-fdr) CIMS_FDR="$2"; shift 2 ;;
        --cits-pval) CITS_PVALUE="$2"; shift 2 ;;
        --cits-gap) CITS_GAP="$2"; shift 2 ;;
        --run-clink) RUN_CLINK="true"; shift ;;
        --clink-umi-len) CLINK_UMI_LEN="$2"; shift 2 ;;
        --clink-fdr) CLINK_FDR="$2"; shift 2 ;;
        --clink-min-cov) CLINK_MIN_COV="$2"; shift 2 ;;
        --clink-multi-map) CLINK_MULTI_MAP="true"; shift ;;
        --no-motif) RUN_MOTIF="no"; shift ;;
        -f|--flank) MOTIF_FLANK="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -m) MISMATCH_MAX="$2"; shift 2 ;;
        -h|--help) show_usage; exit 0 ;;
        *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# CLINK_MULTI_MAP is read from env by run_clink_collapse
export CLINK_MULTI_MAP

# Run Wizard if requested (before validation so it can collect inputs)
if [[ "$WIZARD_MODE" == "true" ]]; then
    run_wizard_mapittymap
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    
    # Apply wizard settings
    INPUT_FILE="$WIZ_INPUT_FILE"
    GENOME_INDEX="$WIZ_GENOME_INDEX"
    [[ -n "$WIZ_GENOME_FASTA" ]] && GENOME_FASTA="$WIZ_GENOME_FASTA"
    [[ "$WIZ_FILTER_REPEAT" == "true" ]] && FILTER_REPEAT="true"
    [[ "$WIZ_NO_CHR_FILTER" == "true" ]] && FILTER_CHR="false"
    ALIGNER="$WIZ_ALIGNER"
    THREADS="$WIZ_THREADS"
    MISMATCH_MAX="${WIZ_ALIGN_MISMATCHES:-2}"
    [[ -n "$WIZ_UMI_LEN" ]] && UMI_LEN="$WIZ_UMI_LEN"
    [[ -n "$WIZ_OUTPUT_NAME" ]] && EXP_ID="$WIZ_OUTPUT_NAME"
    [[ -n "$WIZ_ALIGNER_ARGS" ]] && ADV_ALIGNER_ARGS="$WIZ_ALIGNER_ARGS"

    # Analysis tracks
    [[ "$WIZ_RUN_CIMS" == "true" ]] && RUN_CIMS="true"
    [[ "$WIZ_RUN_CITS" == "true" ]] && RUN_CITS="true"
    [[ -n "$WIZ_CIMS_ITER" ]] && CIMS_ITERATIONS="$WIZ_CIMS_ITER"
    [[ -n "$WIZ_CIMS_FDR" ]]  && CIMS_FDR="$WIZ_CIMS_FDR"
    [[ -n "$WIZ_CITS_PVAL" ]] && CITS_PVALUE="$WIZ_CITS_PVAL"
    [[ -n "$WIZ_CITS_GAP" ]]  && CITS_GAP="$WIZ_CITS_GAP"
    [[ "$WIZ_RUN_CLINK" == "true" ]] && RUN_CLINK="true"
    [[ -n "$WIZ_CLINK_UMI_LEN" ]] && CLINK_UMI_LEN="$WIZ_CLINK_UMI_LEN"
    [[ -n "$WIZ_CLINK_FDR" ]]     && CLINK_FDR="$WIZ_CLINK_FDR"
    [[ -n "$WIZ_CLINK_MIN_COV" ]] && CLINK_MIN_COV="$WIZ_CLINK_MIN_COV"
    [[ "$WIZ_CLINK_MULTI_MAP" == "true" ]] && CLINK_MULTI_MAP="true"
    [[ "$WIZ_NO_MOTIF" == "true" ]] && RUN_MOTIF="no"
    [[ -n "$WIZ_FLANK" ]] && MOTIF_FLANK="$WIZ_FLANK"
fi

if [[ -z "$INPUT_FILE" ]] || [[ -z "$GENOME_INDEX" ]]; then
    log_error "Missing arguments."
    show_usage
    exit 1
fi

check_file "$INPUT_FILE" || exit 1

# Resolve absolute paths
INPUT_FILE="$(cd "$(dirname "$INPUT_FILE")" && pwd)/$(basename "$INPUT_FILE")"
GENOME_INDEX="$(cd "$GENOME_INDEX" && pwd)"

# Resolve and export GENOME_FASTA if provided (used by run_parse_alignment in modules.sh)
if [[ -n "$GENOME_FASTA" ]]; then
    if [[ ! -f "$GENOME_FASTA" ]]; then
        log_error "Genome FASTA not found: $GENOME_FASTA"
        exit 1
    fi
    GENOME_FASTA="$(cd "$(dirname "$GENOME_FASTA")" && pwd)/$(basename "$GENOME_FASTA")"
    export GENOME_FASTA
    log_info "Genome FASTA: $GENOME_FASTA"
fi

# Check Dependencies based on Aligner
if [[ "$ALIGNER" == "bowtie2" ]]; then
    check_dependency bowtie2
    check_bowtie_index "$GENOME_INDEX" || exit 1
else
    check_dependency STAR
    check_star_index "$GENOME_INDEX" || exit 1
fi
check_dependency samtools

show_header
log_info "MAPittyMap: Standalone Mapping Module"
log_info "Mode:      $(if [[ "$ADVANCED_MODE" == "true" ]]; then echo "ADVANCED"; else echo "STANDARD"; fi)"
log_info "Aligner:   $ALIGNER"
log_info "Input:     $INPUT_FILE"
log_info "Genome:    $GENOME_INDEX"

if [[ "$ADVANCED_MODE" == "true" ]]; then
    log_info "Added Args: ${ADV_ALIGNER_ARGS:-(None)}"
fi

# Name Handling
BASENAME=$(basename "$INPUT_FILE" .fastq.gz)
if [[ "$BASENAME" == "$INPUT_FILE" ]]; then BASENAME=$(basename "$INPUT_FILE" .fa); fi
if [[ "$BASENAME" == "$INPUT_FILE" ]]; then BASENAME=$(basename "$INPUT_FILE" .fasta); fi
if [[ -n "$EXP_ID" ]]; then BASENAME="$EXP_ID"; fi

# Structure: Create 1_BAM inside a named folder
if [[ "$BASENAME" == /* ]]; then
    OUT_DIR="${BASENAME}_mapping"
    BASENAME=$(basename "$BASENAME")
else
    OUT_DIR="${BASENAME}_mapping"
fi
mkdir -p "${OUT_DIR}/1_BAM"
mkdir -p "${OUT_DIR}/REPORTS"

# Log Logic: Redirect all subsequent output to a log file in REPORTS
LOG_FILE="${OUT_DIR}/REPORTS/${BASENAME}_mapping.log"
touch "$LOG_FILE"
log_info "Logging to: $LOG_FILE"

# Run Mapping
cd "${OUT_DIR}" || exit 1

# Run repeat element pre-filtering if enabled and index exists
MAPPING_INPUT="$INPUT_FILE"
if [[ "$FILTER_REPEAT" == "true" ]]; then
    REPEAT_INDEX_DIR=$(check_repeat_index "$GENOME_INDEX")
    if [[ -n "$REPEAT_INDEX_DIR" ]]; then
        mkdir -p "OTHERS/Repeat_Mapping"
        REPEAT_UNMAPPED="OTHERS/Repeat_Mapping/${BASENAME}_repeat_filtered.fastq.gz"
        run_repeat_filter "$INPUT_FILE" "$REPEAT_UNMAPPED" "OTHERS/Repeat_Mapping" "$REPEAT_INDEX_DIR" "$THREADS" "$BASENAME"
        run_repeat_quantify \
            "OTHERS/Repeat_Mapping/${BASENAME}_repeat.bam" \
            "OTHERS/Repeat_Mapping/${BASENAME}_repeat_stats.txt" \
            "OTHERS/Repeat_Mapping" "$BASENAME"
        MAPPING_INPUT="$REPEAT_UNMAPPED"
    else
        log_warning "Repeat index not found in $GENOME_INDEX or $GENOME_INDEX/Repeat. Skipping repeat pre-filtering."
    fi
fi

if [[ "$ALIGNER" == "bowtie2" ]]; then
    run_mapping_bowtie2 "$MAPPING_INPUT" "1_BAM/${BASENAME}" "$GENOME_INDEX" "$THREADS"
else
    # Only Star uses Default Mismatches in the arguments, Bowtie2 uses standard sensitivity
    run_mapping_star "$MAPPING_INPUT" "1_BAM/${BASENAME}" "$GENOME_INDEX" "$THREADS" "$MISMATCH_MAX"
fi

# Validate Output
if [[ "$ALIGNER" == "bowtie2" ]]; then
    BAM_FILE="1_BAM/${BASENAME}.bam"
else
    BAM_FILE="1_BAM/${BASENAME}.Aligned.sortedByCoord.out.bam"
fi
if [[ -f "$BAM_FILE" ]]; then
    log_info "Mapping finished successfully."
    log_info "BAM File: ${OUT_DIR}/${BAM_FILE}"

    # Optional: If Advanced Mode, print reminder
    if [[ -n "$ADV_ALIGNER_ARGS" ]]; then
        log_info "Arguments used: $ADV_ALIGNER_ARGS"
    fi
else
    log_error "Mapping failed (BAM not found)."
    exit 1
fi

# ─── POST-MAPPING ANALYSIS ──────────────────────────────────────────────────
# Skip everything if --bam-only.  Otherwise: collapsed BED (2_BED) + bedgraph
# (3_COVERAGE) are always produced; CTK CIMS/CITS (4_CTK_Analysis) and Clink
# (5_Clink_Analysis) run only when their --run-* flag is set.
if [[ "$BAM_ONLY" == "true" ]]; then
    log_info "--bam-only set: skipping post-alignment analysis."
    exit 0
fi

# Canonical chromosome filter (ON by default; --no-chr-filter skips)
if [[ "$FILTER_CHR" == "true" ]]; then
    FILTERED_BAM="1_BAM/${BASENAME}.filtered.bam"
    filter_canonical_chromosomes "$BAM_FILE" "$FILTERED_BAM"
    BAM_FILE="$FILTERED_BAM"
fi

check_dependency bedtools

# PCR dedup via Clink (umi_tools-based) → dedup BAM + collapsed BED (5' positions)
mkdir -p 2_BED 3_COVERAGE 5_Clink_Analysis
CLINK_DEDUP_BAM="5_Clink_Analysis/${BASENAME}_dedup.bam"
log_info "PCR dedup: collapse.py (Clink) → bedtools bamtobed → 2_BED"
run_clink_collapse "$BAM_FILE" "$CLINK_DEDUP_BAM" "$CLINK_UMI_LEN" "$THREADS" || true

COLLAPSED_BED="2_BED/${BASENAME}_collapsed.bed"
bedtools bamtobed -i "$CLINK_DEDUP_BAM" -split 2>/dev/null \
    | sort -k1,1 -k2,2n > "$COLLAPSED_BED"
log_info "Collapsed BED: $COLLAPSED_BED ($(wc -l < "$COLLAPSED_BED") reads)"

# Strand-specific normalized bedgraphs
if [[ -f "$GENOME_INDEX/chrom.sizes" ]]; then
    run_coverage "$COLLAPSED_BED" "3_COVERAGE/${BASENAME}" \
        "$GENOME_INDEX/chrom.sizes" "$CLINK_DEDUP_BAM"
else
    log_warning "chrom.sizes not found in $GENOME_INDEX. Skipping bedgraph."
fi

# ── CTK CIMS/CITS track (opt-in) ──
if [[ "$RUN_CIMS" == "true" || "$RUN_CITS" == "true" ]]; then
    log_info "CTK preprocessing: parseAlignment.pl + tag2collapse.pl"
    check_dependency parseAlignment.pl
    check_dependency tag2collapse.pl
    mkdir -p 4_CTK_Analysis
    CTK_COLLAPSED_BED="4_CTK_Analysis/${BASENAME}_ctk_collapsed.bed"
    MUTATION_FILE="4_CTK_Analysis/${BASENAME}_mutations.txt"
    run_parse_alignment "$BAM_FILE" "4_CTK_Analysis/${BASENAME}_parsed.bed" \
        "$MUTATION_FILE" "$GENOME_INDEX"
    run_collapse_pcr "4_CTK_Analysis/${BASENAME}_parsed.bed" "$CTK_COLLAPSED_BED" \
        "$UMI_LEN" "$DEDUP_MODE"

    # Resolve reference FASTA for motif analysis (same priority as CLIPittyClip.sh)
    ref_fasta="$GENOME_FASTA"
    if [[ -z "$ref_fasta" ]]; then
        ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \
            \( -name "*genome*.fa" -o -name "*genome*.fasta" \) 2>/dev/null | head -n 1)
    fi
    if [[ -z "$ref_fasta" ]]; then
        ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \
            \( -name "*primary*.fa" -o -name "*primary*.fasta" \) 2>/dev/null | head -n 1)
    fi
    if [[ -z "$ref_fasta" ]]; then
        ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \
            \( -name "*.fa" -o -name "*.fasta" \) ! -name "*repeat*" 2>/dev/null | head -n 1)
    fi
    [[ -z "$ref_fasta" ]] && log_warning "Reference FASTA not found. Motif analysis may be skipped."

    run_ctk_analysis \
        "$CTK_COLLAPSED_BED" "$MUTATION_FILE" "4_CTK_Analysis" "$ref_fasta" "$BASENAME" \
        "$CIMS_ITERATIONS" "$CIMS_FDR" "$CITS_PVALUE" "$CITS_GAP" \
        "$MOTIF_FLANK" "$RUN_MOTIF" "$RUN_CIMS" "$RUN_CITS"
    log_info "CTK Analysis complete. Output: 4_CTK_Analysis/"
fi

# ── Clink track (opt-in) ──
if [[ "$RUN_CLINK" == "true" ]]; then
    log_info "Running Clink pipeline..."
    if ! check_clink_deps; then
        log_error "Clink dependencies missing. Skipping Clink pipeline."
    else
        # run_clink_full receives the pre-built dedup BAM so it skips re-deduplication
        run_clink_full "$BAM_FILE" "5_Clink_Analysis" "$BASENAME" \
            "$CLINK_UMI_LEN" "$THREADS" "true" "true" \
            "$CLINK_MIN_COV" "$CLINK_MIN_FRAC" "$CLINK_FDR" "$CLINK_DEDUP_BAM"
        log_info "Clink Analysis complete. Output: 5_Clink_Analysis/"
    fi
fi

log_info "MAPittyMap done."
