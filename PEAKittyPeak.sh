#!/bin/bash

# PEAKittyPeak.sh - Peak Calling Module for CLIPittyClip (v3.0)
# Uses the unified lib/utils.sh for logging

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Source utils and WIZARD (just in case)
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/wizard.sh"

# Default Values
PEAK_DIST=50
PEAK_SIZE=20
FRAG_LEN=25
BASE_NAME="Combined"
ADV_HOMER_ARGS="" # Additional arguments for findPeaks
WIZARD_MODE="false"

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
    echo "  -p <int>       Min distance between peaks (default: 50)"
    echo "  -z <int>       Peak size (default: 20)"
    echo "  -f <int>       Fragment length (default: 25)"
    echo "  -n <str>       Base name for output (default: 'Combined')"
    echo "  -a <str>       Additional HOMER findPeaks arguments (quoted string)"
    echo "  --wizard       Launch interactive configuration wizard"
    echo "  --advanced     Alias for --wizard (backward compatibility)"
    echo "  -h, --help     Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  # Basic peak calling"
    echo "  PEAKittyPeak.sh -p 50 -z 20 -n MyExperiment"
    echo ""
    echo "  # With custom HOMER args"
    echo "  PEAKittyPeak.sh -n Combined -a '-style factor -L 2'"
    echo ""
}

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Parse Options
while getopts "p:z:f:n:a:h-:" opt; do
  # Handle long options manually
  if [[ "$opt" == "-" ]]; then
      case "${OPTARG}" in
          wizard|advanced) WIZARD_MODE="true" ;;
          *) log_error "Invalid option --${OPTARG}"; show_usage; exit 1 ;;
      esac
      continue
  fi

  case $opt in
    p) PEAK_DIST="$OPTARG" ;;
    z) PEAK_SIZE="$OPTARG" ;;
    f) FRAG_LEN="$OPTARG" ;;
    n) BASE_NAME="$OPTARG" ;;
    a) ADV_HOMER_ARGS="$OPTARG" ;;
    h) show_usage; exit 0 ;;
    \?) log_error "Invalid option -$OPTARG"; show_usage; exit 1 ;;
  esac
done

# Run Wizard if requested (before validation so it can collect inputs)
if [[ "$WIZARD_MODE" == "true" ]]; then
    run_wizard_peakittypeak
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    
    # Apply wizard settings
    [[ -n "$WIZ_WORK_DIR" ]] && cd "$WIZ_WORK_DIR"
    PEAK_DIST="$WIZ_PEAK_DIST"
    PEAK_SIZE="$WIZ_PEAK_SIZE"
    FRAG_LEN="$WIZ_FRAG_LEN"
    BASE_NAME="$WIZ_OUTPUT_NAME"
    ADV_HOMER_ARGS="$WIZ_HOMER_ARGS"
fi

# Check Requirements
if [[ ! -d "BED" ]]; then
    log_error "Directory 'BED' not found. Please run this script in the parent folder of your BED files."
    exit 1
fi

BED_COUNT=$(ls BED/*.bed 2>/dev/null | wc -l)
if [[ $BED_COUNT -eq 0 ]]; then
    log_error "No .bed files found in BED/ folder."
    exit 1
fi

log_info "PEAKittyPeak: Peak Calling Module"
log_info "Mode:        $(if [[ "$ADVANCED_MODE" == "true" ]]; then echo "ADVANCED"; else echo "STANDARD"; fi)"
log_info "Input Files: $BED_COUNT"
log_info "Parameters:  Dist=$PEAK_DIST, Size=$PEAK_SIZE, Frag=$FRAG_LEN"
if [[ -n "$ADV_HOMER_ARGS" ]]; then
    log_info "Advanced Args: $ADV_HOMER_ARGS"
fi

# 1. Combine BED Files
log_info "Combining BED files..."
cat BED/*.bed > "${BASE_NAME}.bed"

# 2. Call Peaks (HOMER)
log_info "Calling peaks with HOMER..."
makeTagDirectory "${BASE_NAME}_TagDir/" "${BASE_NAME}.bed" -single -format bed > /dev/null 2>&1

CMD="findPeaks ${BASE_NAME}_TagDir/ -o auto -style factor -L 2 -localSize 10000 -strand separate \
    -minDist ${PEAK_DIST} -size ${PEAK_SIZE} -fragLength ${FRAG_LEN} ${ADV_HOMER_ARGS}"

log_info "Running: $CMD"
eval "$CMD" 2>&1 | grep -v "Job finished" # Reduce noise

# 3. Process Peaks
log_info "Processing peaks output..."
if [[ -f "${BASE_NAME}_TagDir/peaks.txt" ]]; then
    # Convert peaks.txt to BED format
    sed '/^[[:blank:]]*#/d;s/#.*//' "${BASE_NAME}_TagDir/peaks.txt" > peaksTemp.bed
    awk 'OFS="\t" {print $2, $3, $4, $1, $6, $5}' peaksTemp.bed > peaks.bed
    rm peaksTemp.bed
    
    # Sort
    sort -k 1,1 -k2,2n peaks.bed > peaks_Sorted.bed
    
    # Organize Output
    OUT_DIR="${BASE_NAME}_peaks"
    mkdir -p "$OUT_DIR"
    mkdir -p "${OUT_DIR}/peakCoverage"
    
    mv "${BASE_NAME}_TagDir" "$OUT_DIR"
    mv "${BASE_NAME}.bed" "$OUT_DIR"
    mv peaks.bed "$OUT_DIR"
    mv peaks_Sorted.bed "$OUT_DIR" # Keep sorted version in root
    
    # 4. Coverage Analysis
    log_info "Calculating coverage..."
    # Copy sorted bed to serve as base table
    cp "${OUT_DIR}/peaks_Sorted.bed" "${OUT_DIR}/${BASE_NAME}_peakCoverage.txt"
    
    # Header
    echo -e "chr\tstart\tend\tname\tscore\tstrand" > colnames.txt
    cat "${OUT_DIR}/${BASE_NAME}_peakCoverage.txt" >> colnames.txt
    mv colnames.txt "${OUT_DIR}/${BASE_NAME}_peakCoverage.txt"
    
    for bed_file in BED/*.bed; do
        s_name=$(basename "$bed_file" .bed)
        s_name=$(basename "$s_name" .collapsed) # Strip .collapsed if present
        
        bedtools coverage -s -a "${OUT_DIR}/peaks_Sorted.bed" -b "$bed_file" > "coverage_${s_name}.txt"
        
        # Extract coverage column (7th column in bedtools coverage default output)
        awk 'FNR>0 {print $7}' "coverage_${s_name}.txt" > temp_cov.txt
        
        # Append name to temp header
        echo "$s_name" > temp_col.txt
        cat temp_cov.txt >> temp_col.txt
        
        # Paste to main table
        paste "${OUT_DIR}/${BASE_NAME}_peakCoverage.txt" temp_col.txt > temp_table.txt
        mv temp_table.txt "${OUT_DIR}/${BASE_NAME}_peakCoverage.txt"
        
        # Cleanup temp
        rm temp_cov.txt temp_col.txt
        
        # Move raw coverage file
        mv "coverage_${s_name}.txt" "${OUT_DIR}/peakCoverage/"
    done
    
    log_info "Analysis Finished."
    log_info "Outputs stored in: $OUT_DIR"
else
    log_error "Peak calling failed (peaks.txt not found)."
    exit 1
fi