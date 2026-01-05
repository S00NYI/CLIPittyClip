#!/bin/bash

# MAPittyMap.sh - Mapping-only module for CLIPittyClip (v3.0)
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
WIZARD_MODE="false"

function show_usage {
    echo ""
    echo "Usage: $0 -i <input.fastq.gz> -x <index_dir> [options]"
    echo ""
    echo "MAPittyMap v3.0 - Standalone Mapping Module for CLIPittyClip"
    echo ""
    echo "REQUIRED:"
    echo "  -i <path>      Input FASTQ/FASTA file (gzipped supported)"
    echo "  -x <path>      Path to genome index directory"
    echo ""
    echo "OPTIONS:"
    echo "  --aligner <str>  Aligner: 'star' (default) or 'bowtie2'"
    echo "  -o <str>         Output ID (default: derived from filename)"
    echo "  -t <int>         Number of threads (default: 1)"
    echo "  -m <int>         Max mismatches (default: 2)"
    echo "  --wizard         Launch interactive configuration wizard
  --advanced       Alias for --wizard (backward compatibility)
  -h, --help       Show this help message"
    echo ""
    echo "OUTPUT:"
    echo "  Creates 1_BAM/ directory with sorted, indexed BAM file"
    echo "  Creates REPORTS/ directory with alignment logs"
    echo ""
    echo "EXAMPLES:"
    echo "  # Map with STAR"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x /path/to/star_index -t 8"
    echo ""
    echo "  # Map with Bowtie2"
    echo "  MAPittyMap.sh -i reads.fastq.gz -x /path/to/bt2_index -t 8 --aligner bowtie2"
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
        --wizard|--advanced) WIZARD_MODE="true"; shift 1 ;;
        -t) THREADS="$2"; shift 2 ;;
        -m) MISMATCH_MAX="$2"; shift 2 ;;
        -h|--help) show_usage; exit 0 ;;
        *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# Run Wizard if requested (before validation so it can collect inputs)
if [[ "$WIZARD_MODE" == "true" ]]; then
    run_wizard_mapittymap
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    
    # Apply wizard settings
    INPUT_FILE="$WIZ_INPUT_FILE"
    GENOME_INDEX="$WIZ_GENOME_INDEX"
    ALIGNER="$WIZ_ALIGNER"
    THREADS="$WIZ_THREADS"
    MISMATCH_MAX="$WIZ_MISMATCHES"
    [[ -n "$WIZ_OUTPUT_NAME" ]] && EXP_ID="$WIZ_OUTPUT_NAME"
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
OUT_DIR="${BASENAME}_mapping"
mkdir -p "${OUT_DIR}/1_BAM"
mkdir -p "${OUT_DIR}/REPORTS"

# Log Logic: Redirect all subsequent output to a log file in REPORTS
LOG_FILE="${OUT_DIR}/REPORTS/${BASENAME}_mapping.log"
touch "$LOG_FILE"
log_info "Logging to: $LOG_FILE"

# Run Mapping
cd "${OUT_DIR}" || exit 1

if [[ "$ALIGNER" == "bowtie2" ]]; then
    run_mapping_bowtie2 "$INPUT_FILE" "1_BAM/${BASENAME}" "$GENOME_INDEX" "$THREADS"
else
    # Only Star uses Default Mismatches in the arguments, Bowtie2 uses standard sensitivity
    run_mapping_star "$INPUT_FILE" "1_BAM/${BASENAME}" "$GENOME_INDEX" "$THREADS" "$MISMATCH_MAX"
fi

# Validate Output
BAM_FILE="1_BAM/${BASENAME}.Aligned.sortedByCoord.out.bam"
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
