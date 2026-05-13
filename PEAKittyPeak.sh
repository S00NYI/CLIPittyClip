#!/bin/bash

# PEAKittyPeak.sh - Peak Calling Module for CLIPittyClip (v3.0)
# Uses the unified lib/utils.sh for logging

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Source utils, modules, and wizard
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/modules.sh"
source "${SCRIPT_DIR}/lib/wizard.sh"

# Fallback: if LOG_FILE not set by parent pipeline, discard log output
LOG_FILE="${LOG_FILE:-/dev/null}"

# Default Values
PEAK_DIST=50
PEAK_SIZE=20
FRAG_LEN=25
BASE_NAME="Combined"
ADV_PEAK_CALLER_ARGS="" # Additional arguments for the peak caller
PEAK_CALLER="homer" # Peak caller: homer (default) or ctk
LOG_FILE="${LOG_FILE:-/dev/null}" # Use parent's log if set, otherwise discard
WIZARD_MODE="false"
CTK_DIR=""         # Optional: CTK analysis output directory
CTK_GROUPS_FILE="" # Optional: groups file for CTK aggregation
CIMS_FDR="0.05"    # Default FDR threshold for CIMS
CITS_PVALUE="0.05" # Default p-value threshold for CITS

function show_usage {
    echo ""
    echo "Usage: PEAKittyPeak.sh [options]"
    echo ""
    echo "PEAKittyPeak v3.0 - Peak Calling Module for CLIPittyClip"
    echo ""
    echo "CONTEXT:"
    echo "  Run this in a directory containing a 'BED' folder with collapsed .bed files."
    echo "  Typically called automatically by CLIPittyClip.sh after sample processing."
    echo ""
    echo "OPTIONS:"
    echo "  -i, --input <dir> Input directory containing BED files (default: BED)"
    echo "  --aggregate    Enable Aggregation Mode (combine all inputs)"
    echo "  --no-aggregate Disable Aggregation Mode (process individually)"
    echo "  -p <int>       Min distance between peaks (default: 50)"
    echo "  -z <int>       Peak size (default: 20)"
    echo "  -f <int>       Fragment length (default: 25)"
    echo "  -n <str>       Base name for output (default: 'Combined')"
    echo "  --peak-caller-args <str> Additional peak caller arguments (quoted string)"
    echo "  --peak-caller <str> Peak caller: homer (default) or ctk"
    echo "  --ctk-dir <path>   Add CIMS/CITS site counts from CTK analysis"
    echo "  -g, --groups <file> Groups file for aggregation metrics"
    echo "  --cims-fdr <float> CIMS FDR threshold (default: 0.05)"
    echo "  --cits-pval <float> CITS p-value threshold (default: 0.05)"
    echo "  --wizard       Launch interactive configuration wizard"
    echo "  --advanced     Alias for --wizard (backward compatibility)"
    echo "  -h, --help     Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  # Basic peak calling (Individual Mode)"
    echo "  PEAKittyPeak.sh -i ./bed_files/"
    echo ""
    echo "  # Aggregated peak calling"
    echo "  PEAKittyPeak.sh -i ./bed_files/ --aggregate -n AllSamples"
    echo ""
}

if [[ $# -eq 0 ]]; then
    show_usage
    exit 1
fi

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Parse Options using while loop for better long option handling
INPUT_BED_DIR="BED"  # Default used if -i not provided
AGGREGATE="false"    # Default: Individual mode
GROUPS_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input) INPUT_BED_DIR="$2"; shift 2 ;;
        --aggregate) AGGREGATE="true"; shift ;;
        --no-aggregate) AGGREGATE="false"; shift ;;
        -p) PEAK_DIST="$2"; shift 2 ;;
        -z) PEAK_SIZE="$2"; shift 2 ;;
        -f) FRAG_LEN="$2"; shift 2 ;;
        -n) BASE_NAME="$2"; shift 2 ;;
        --peak-caller-args) ADV_PEAK_CALLER_ARGS="$2"; shift 2 ;;
        --peak-caller) PEAK_CALLER="$2"; shift 2 ;;
        --ctk-dir) CTK_DIR="$2"; shift 2 ;;
        -g|--groups|--group-xlsite) GROUPS_FILE="$2"; shift 2 ;;
        --cims-fdr) CIMS_FDR="$2"; shift 2 ;;
        --cits-pval) CITS_PVALUE="$2"; shift 2 ;;
        --wizard|--advanced) WIZARD_MODE="true"; shift ;;
        -h|--help) show_usage; exit 0 ;;
        *) log_error "Invalid option: $1"; show_usage; exit 1 ;;
    esac
done

# Run Wizard if requested
if [[ "$WIZARD_MODE" == "true" ]]; then
    run_wizard_peakittypeak
    if [[ $? -ne 0 ]]; then exit 1; fi
    [[ -n "$WIZ_WORK_DIR" ]] && cd "$WIZ_WORK_DIR"
    PEAK_DIST="$WIZ_PEAK_DIST"
    PEAK_SIZE="$WIZ_PEAK_SIZE"
    FRAG_LEN="$WIZ_FRAG_LEN"
    BASE_NAME="$WIZ_OUTPUT_NAME"
    ADV_PEAK_CALLER_ARGS="$WIZ_HOMER_ARGS"
    
    if [[ -n "$WIZ_PEAK_CALLER" ]]; then
        PEAK_CALLER="$WIZ_PEAK_CALLER"
    fi
    
    if [[ -n "$WIZ_CTK_PEAK_ARGS" ]]; then
        ADV_PEAK_CALLER_ARGS="$WIZ_CTK_PEAK_ARGS"
    fi
    # Wizard currently assumes CWD/BED, can be updated later if needed
fi

# Check Requirements
if [[ ! -d "$INPUT_BED_DIR" ]]; then
    log_error "Input directory '$INPUT_BED_DIR' not found."
    exit 1
fi

BED_COUNT=$(ls "$INPUT_BED_DIR"/*.bed 2>/dev/null | wc -l)
if [[ $BED_COUNT -eq 0 ]]; then
    log_error "No .bed files found in '$INPUT_BED_DIR'."
    exit 1
fi

log_info "PEAKittyPeak: Peak Calling Module"
log_info "Input Dir:   $INPUT_BED_DIR"
log_info "Input Files: $BED_COUNT"
log_info "Mode:        $(if [[ "$AGGREGATE" == "true" ]]; then echo "AGGREGATE"; else echo "INDIVIDUAL"; fi)"
log_info "Caller:      $PEAK_CALLER"
log_info "Parameters:  Dist=$PEAK_DIST, Size=$PEAK_SIZE, Frag=$FRAG_LEN"

# --- Function: Call Peaks ---
call_peaks() {
    local input_file="$1"
    local output_name="$2"
    local out_dir="${output_name}_peaks"

    log_info "Processing: $output_name"

    # ----------------------------------------------------------------
    # Peak calling: HOMER or CTK
    # ----------------------------------------------------------------
    if [[ "${PEAK_CALLER:-homer}" == "ctk" ]]; then
        # --- CTK tag2peak.pl ---
        local raw_peaks="peaks_raw.bed"
        local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/tag2peak_cache.XXXXXX")

        log_info "Running: tag2peak.pl ..."
        $CONDA_PREFIX/bin/perl $(which tag2peak.pl) -big -ss --valley-seeking -minPH 2 -gap ${PEAK_DIST} \
            ${ADV_PEAK_CALLER_ARGS} -c "${cache_dir}" "${input_file}" "${raw_peaks}" 2>&1 | grep -v "^CMD="
        local exit_code=${PIPESTATUS[0]}
        rm -rf "$cache_dir"

        if [[ $exit_code -ne 0 || ! -s "$raw_peaks" ]]; then
            log_error "No peaks generated for $output_name"
            rm -f "$raw_peaks"
            return
        fi

        sort -k1,1 -k2,2n "$raw_peaks" > peaks_Sorted.bed
        rm -f "$raw_peaks"

        mkdir -p "$out_dir"
        mkdir -p "${out_dir}/peakCoverage"
        cp "$input_file" "$out_dir/${output_name}.bed"
        mv peaks_Sorted.bed "$out_dir/"

    else
        # --- HOMER ---
        local tag_dir="${output_name}_TagDir"

        makeTagDirectory "$tag_dir/" "$input_file" -single -format bed > /dev/null 2>&1

        local cmd="findPeaks $tag_dir/ -o auto -style factor -L 2 -localSize 1000 -strand separate \
            -minDist ${PEAK_DIST} -size ${PEAK_SIZE} -fragLength ${FRAG_LEN} ${ADV_PEAK_CALLER_ARGS}"

        log_info "Running: findPeaks ..."
        eval "$cmd" 2>&1 | grep -v "Job finished"

        if [[ ! -f "$tag_dir/peaks.txt" ]]; then
            log_error "No peaks generated for $output_name"
            return
        fi

        sed '/^[[:blank:]]*#/d;s/#.*//' "$tag_dir/peaks.txt" > peaksTemp.bed
        awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
        rm peaksTemp.bed
        sort -k1,1 -k2,2n peaks.bed > peaks_Sorted.bed

        mkdir -p "$out_dir"
        mkdir -p "${out_dir}/peakCoverage"
        mv "$tag_dir" "$out_dir/"
        cp "$input_file" "$out_dir/${output_name}.bed"
        mv peaks.bed "$out_dir/"
        mv peaks_Sorted.bed "$out_dir/"
    fi

    # ----------------------------------------------------------------
    # Shared: Coverage Analysis (same for both callers)
    # ----------------------------------------------------------------
    local coverage_file="${out_dir}/${output_name}_PEAK_MATRIX.txt"
    cp "${out_dir}/peaks_Sorted.bed" "$coverage_file"

    HEADER_STR="chr\tstart\tend\tname\tscore\tstrand"

    log_info "Calculating per-sample coverage (bedtools)..."
    for bed_file in "$INPUT_BED_DIR"/*.bed; do
        if [[ -f "$bed_file" ]]; then
            local s_name=$(basename "$bed_file" .bed)
            s_name=${s_name%_collapsed}
            if [[ "$s_name" == "$output_name" ]]; then continue; fi

            bedtools coverage -s -a "${out_dir}/peaks_Sorted.bed" -b "$bed_file" > "temp_cov.txt"
            awk '{print $7}' "temp_cov.txt" > "col_count.txt"
            paste "$coverage_file" "col_count.txt" > "temp_paste.txt"
            mv "temp_paste.txt" "$coverage_file"
            HEADER_STR="${HEADER_STR}\tTC_${s_name}"
            rm "temp_cov.txt" "col_count.txt" 2>/dev/null
        fi
    done

    echo -e "$HEADER_STR" > colnames.txt
    cat colnames.txt "$coverage_file" > temp_final.txt
    mv temp_final.txt "$coverage_file"
    rm colnames.txt

    if [[ -n "$CTK_DIR" ]] && [[ "${MATRIX_CTK_COLS:-false}" == "true" ]]; then
        log_info "Adding CTK columns..."
        add_ctk_columns_to_peak_matrix "$coverage_file" "${out_dir}/peaks_Sorted.bed" "$CTK_DIR" "$CIMS_FDR" "$CITS_PVALUE" "$GROUPS_FILE"
    elif [[ -n "$CTK_DIR" ]]; then
        log_info "Skipping CTK columns (MATRIX_CTK_COLS=false)"
    fi

    log_info "Peak calling for $output_name complete: $out_dir/"
}


# --- Execution Flow ---

if [[ "$AGGREGATE" == "true" ]]; then
    # Aggregate Mode: Combine all BEDs and call once
    log_info "Mode: AGGREGATE - Combining all .bed files..."
    cat "$INPUT_BED_DIR"/*.bed > "${BASE_NAME}.bed"
    call_peaks "${BASE_NAME}.bed" "${BASE_NAME}"
    rm "${BASE_NAME}.bed" 2>/dev/null
    
else
    # Individual Mode: Loop through files
    log_info "Mode: INDIVIDUAL - Processing files separately..."
    # Ensure shell expansion works
    shopt -s nullglob
    for bed_file in "$INPUT_BED_DIR"/*.bed; do
        if [[ -f "$bed_file" ]]; then
            bname=$(basename "$bed_file" .bed)
            call_peaks "$bed_file" "$bname"
        fi
    done
    shopt -u nullglob
fi

log_info "PEAKittyPeak analysis complete."
exit 0