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
ADV_HOMER_ARGS="" # Additional arguments for findPeaks
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
    echo "  -a <str>       Additional HOMER findPeaks arguments (quoted string)"
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
        -a) ADV_HOMER_ARGS="$2"; shift 2 ;;
        --ctk-dir) CTK_DIR="$2"; shift 2 ;;
        -g|--groups|--ctk-group) GROUPS_FILE="$2"; shift 2 ;;
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
    ADV_HOMER_ARGS="$WIZ_HOMER_ARGS"
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
log_info "Parameters:  Dist=$PEAK_DIST, Size=$PEAK_SIZE, Frag=$FRAG_LEN"

# --- Function: Call Peaks ---
call_peaks() {
    local input_file="$1"
    local output_name="$2"
    
    log_info "Processing: $output_name"
    
    # 1. Provide a temporary directory for tag outputs to avoid conflicts
    local tag_dir="${output_name}_TagDir"
    
    # 2. Call Peaks (HOMER)
    makeTagDirectory "$tag_dir/" "$input_file" -single -format bed > /dev/null 2>&1
    
    local cmd="findPeaks $tag_dir/ -o auto -style factor -L 2 -localSize 10000 -strand separate \
        -minDist ${PEAK_DIST} -size ${PEAK_SIZE} -fragLength ${FRAG_LEN} ${ADV_HOMER_ARGS}"
        
    log_info "Running: findPeaks ..."
    # Capture output but reduce noise
    eval "$cmd" 2>&1 | grep -v "Job finished"
    
    # 3. Process Peaks
    # findPeaks outputs to peaks.txt in the tag directory usually, unless -o specified
    # Using -o auto, HOMER usually outputs to tagDir/peaks.txt
    
    if [[ -f "$tag_dir/peaks.txt" ]]; then
        # Convert peaks.txt to BED format
        sed '/^[[:blank:]]*#/d;s/#.*//' "$tag_dir/peaks.txt" > peaksTemp.bed
        awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
        rm peaksTemp.bed
        
        # Sort
        sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed
        
        # Organize Output
        local out_dir="${output_name}_peaks"
        mkdir -p "$out_dir"
        mkdir -p "${out_dir}/peakCoverage"
        
        mv "$tag_dir" "$out_dir/"
        # Copy input bed to output to preserve input
        cp "$input_file" "$out_dir/${output_name}.bed"
        mv peaks.bed "$out_dir/"
        mv peaks_Sorted.bed "$out_dir/"
        
        # 4. Coverage Analysis
        local coverage_file="${out_dir}/${output_name}_peakMatrix.txt"
        cp "${out_dir}/peaks_Sorted.bed" "$coverage_file"
        
        # Header String (TAB separated)
        # Note: We construct this now but apply it at the very end
        # Header String (TAB separated) - Initialize with BED fields
        HEADER_STR="chr\tstart\tend\tname\tscore\tstrand"
        
        # 4a. Per-Sample Coverage (Restored v2 Logic)
        log_info "Calculating per-sample coverage (bedtools)..."
        
        # Iterate over all BED files in the input directory
        # This restores the "Read Counts per Sample" columns
        for bed_file in "$INPUT_BED_DIR"/*.bed; do
            if [[ -f "$bed_file" ]]; then
                # Extract Sample Name
                local s_name=$(basename "$bed_file" .bed)
                s_name=${s_name%_collapsed} # Strip _collapsed if present
                
                # Skip if it matches the Combined output itself (to avoid self-inclusion)
                if [[ "$s_name" == "$output_name" ]]; then continue; fi
                
                # bedtools coverage -s (force strandedness per user v2)
                # -a: Peaks, -b: Sample Reads
                bedtools coverage -s -a "${out_dir}/peaks_Sorted.bed" -b "$bed_file" > "temp_cov.txt"
                
                # Extract count column (Column 7 in bedtools coverage default output)
                awk '{print $7}' "temp_cov.txt" > "col_count.txt"
                
                # Paste to Matrix
                paste "$coverage_file" "col_count.txt" > "temp_paste.txt"
                mv "temp_paste.txt" "$coverage_file"
                
                # Append Header Name with TC_ prefix (Tag Count)
                HEADER_STR="${HEADER_STR}\tTC_${s_name}"
                
                # Cleanup
                rm "temp_cov.txt" "col_count.txt" 2>/dev/null
            fi
        done
        
        # [BUG FIX] Prepend header NOW so downstream functions receive valid matrix
        # Use simple echo instead of file
        echo -e "$HEADER_STR" > colnames.txt
        cat colnames.txt "$coverage_file" > temp_final.txt
        mv temp_final.txt "$coverage_file"
        rm colnames.txt
        
        # 5. [DISABLED] Group Coverage Columns (Sum/Avg) - Removed per user request
        # if [[ -n "$GROUPS_FILE" ]]; then ... fi
        
        # 6. Add CTK CIMS/CITS columns if requested
        if [[ -n "$CTK_DIR" ]]; then
            log_info "Adding CTK columns..."
            # Pass correct arguments: Matrix, PeaksBED, CTK_Dir, Thresholds, AND GroupsFile
            add_ctk_columns_to_peak_matrix "$coverage_file" "${out_dir}/peaks_Sorted.bed" "$CTK_DIR" "$CIMS_FDR" "$CITS_PVALUE" "$GROUPS_FILE"
        fi
        
        log_info "Peak calling for $output_name complete: $out_dir/"
    else
        log_error "No peaks generated for $output_name"
    fi
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