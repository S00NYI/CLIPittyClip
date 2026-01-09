#!/bin/bash

# CLIPittyClip.sh - Modernized CLIP-seq Analysis Pipeline
# Version 3.0.0

# ------------------------------------------------------------------
# Initialization & Setup
# ------------------------------------------------------------------

# Resolve script directory to source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/modules.sh"
source "${SCRIPT_DIR}/lib/wizard.sh"

# Default Values
THREADS=1
KEEP_INTERMEDIATE="no"
GENOME_INDEX=""
INPUT_FILE=""
EXP_ID=""     # Optional, derived from filename if not provided
EXP_TYPE=""   # Optional, derived from filename if not provided
DEMUX="no"    # handled by fastp automatically for UMI, but for sample demux we stick to single sample default per run for now or need bc file
BARCODE_FILE=""
INPUT_DIR=""  # Directory mode for pre-demultiplexed FASTQs
UMI_LEN=0
ADAPTER_3="GTGTCAGTCACTTCCAGCGG" # L32 default
PEAK_DIST=50
PEAK_SIZE=20
FRAG_LEN=25
MISMATCH_MAX=2 # STAR default
# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------

function show_usage {
    echo ""
    echo "Usage: $0 [-i <input.fastq.gz> | -d <input_dir>] -x <index_dir> [options]"
    echo ""
    echo "CLIPittyClip v3.0 - Modern CLIP-seq Analysis Pipeline"
    echo ""
    echo "REQUIRED INPUT (Choose one):"
    echo "  -i <path>          Input FASTQ file (gzipped supported)"
    echo "  -d <dir>           Input directory containing .fastq.gz files (Batch Mode)"
    echo ""
    echo "REQUIRED REFERENCE:"
    echo "  -x <path>          Path to genome index directory (STAR or Bowtie2)"
    echo ""
    echo "GENERAL OPTIONS:"
    echo "  -o <str>           Output ID / prefix (default: derived from filename)"
    echo "  -t <int>           Number of threads (default: 1)"
    echo "  --aligner <str>    Aligner: 'star' (default) or 'bowtie2'"
    echo "  -h, --help         Show this help message"
    echo "  --verbose          Enable verbose logging"
    echo ""
    echo "PREPROCESSING OPTIONS:"
    echo "  -u <int>           UMI length (e.g., 7 for CoCLIP)"
    echo "  -a <str>           3' adapter sequence (default: L32)"
    echo "  --dedup            Enable FASTQ deduplication (default: ON)"
    echo "  --no-dedup         Disable FASTQ deduplication"
    echo ""
    echo "DEMULTIPLEXING OPTIONS (for -i mode):"
    echo "  -b, --barcode <path>   Barcode file for demultiplexing"
    echo "  --mismatches <int>     Max barcode mismatches (default: 1)"
    echo ""
    echo "ANALYSIS OPTIONS:"
    echo "  --run-ctk          Enable full CTK CIMS+CITS analysis"
    echo "  --cims             Enable CIMS analysis (mutation sites only)"
    echo "  --cits             Enable CITS analysis (truncation sites only)"
    echo "  --cims-iter <int>  CIMS permutation iterations (default: 10)"
    echo "  --cims-fdr <float> CIMS FDR threshold (default: 0.05)"
    echo "  --cits-pval <float> CITS p-value threshold (default: 0.05)"
    echo "  --cits-gap <int>   CITS clustering gap (default: 25, -1=no cluster)"
    echo "  --motif-flank <int> Motif flanking nucleotides (default: 10)"
    echo "  --no-motif         Skip motif enrichment analysis"
    echo "  -g, --ctk-group <file> Aggregate samples by group for CIMS/CITS"
    echo "                      Format: sample_name<TAB>group_name"
    echo "  --sample <int>     Test mode: process only first N reads"
    echo "  --skip-ncrna       Disable ncRNA pre-filtering (on by default)"
    echo ""
    echo "OUTPUT OPTIONS:"
    echo "  -k                 Keep intermediate files"
    echo "  --notification     Enable system notifications on completion"
    echo "  --wizard           Launch interactive configuration wizard
  --advanced         Alias for --wizard (backward compatibility)"
    echo ""
    echo "EXAMPLES:"
    echo "  # Standard run (Single File)"
    echo "  $0 -i reads.fastq.gz -x /path/to/star_index -t 8 -u 7"
    echo ""
    echo "  # Demultiplexing run"
    echo "  $0 -i pool.fastq.gz -b barcodes.txt -x /path/to/star_index -t 8"
    echo ""
    echo "  # Directory Batch Mode"
    echo "  $0 -d /path/to/samples/ -x /path/to/star_index -t 8 --cims"
    echo ""
    echo "  # Bowtie2 Alignment"
    echo "  $0 -i reads.fastq.gz -x /path/to/bt2_index -t 8 --aligner bowtie2"
    echo ""
}

# ------------------------------------------------------------------
# Argument Parsing
# ------------------------------------------------------------------

RUN_CIMS=false
RUN_CITS=false
VERBOSE=false
SAMPLE_SIZE=0 
CHILD_MODE="false"
DEDUP_MODE="true" # Always on by default per user request
NOTIFY_MODE="false"
WIZARD_MODE="false"
SKIP_NCRNA="false"  # ncRNA filtering is ON by default
ALIGNER="star" # Default match

# CTK CIMS/CITS Parameters (with defaults)
CIMS_ITERATIONS="10"
CIMS_FDR="0.05"
CITS_PVALUE="0.05"
CITS_GAP="25"
RUN_MOTIF="yes"
MOTIF_FLANK="10"
CTK_GROUPS_FILE=""  # Optional: group samples for CIMS/CITS aggregation

# Capture Start Time (Seconds) for duration calculation
PIPELINE_START=$(date +%s)

if [[ $# -eq 0 ]]; then
    # Only show header if we are going to exit anyway
    echo "$separator_line"
    echo -e "${BLUE}CLIPittyClip: Modern CLIP-seq Analysis Pipeline${NC}"
    echo "Version 3.0.0"
    echo "Author: Soon Yi (Updated by Antigravity)"
    echo "Last updated: $(date +'%Y-%m-%d')"
    echo "$separator_line"
    
    echo "Welcome to CLIPittyClip! No arguments provided."
    echo "Please provide arguments via command line."
    show_usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    # Skip empty arguments
    if [[ -z "$1" ]]; then shift; continue; fi
    case $1 in
        -i) INPUT_FILE="$2"; shift 2 ;;
        -x) GENOME_INDEX="$2"; shift 2 ;;
        --aligner) ALIGNER=$(echo "$2" | tr '[:upper:]' '[:lower:]'); shift 2 ;;
        -o) EXP_ID="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -u) UMI_LEN="$2"; shift 2 ;;
        -a) ADAPTER_3="$2"; shift 2 ;;
        -k) KEEP_INTERMEDIATE="yes"; shift ;;
        --cims) RUN_CIMS=true; RUN_CTK="yes"; shift ;;
        --cits) RUN_CITS=true; RUN_CTK="yes"; shift ;;
        --run-ctk) RUN_CTK="yes"; RUN_CIMS=true; RUN_CITS=true; shift ;;
        --cims-iter) CIMS_ITERATIONS="$2"; shift 2 ;;
        --cims-fdr) CIMS_FDR="$2"; shift 2 ;;
        --cits-pval) CITS_PVALUE="$2"; shift 2 ;;
        --cits-gap) CITS_GAP="$2"; shift 2 ;;
        --motif-flank) MOTIF_FLANK="$2"; shift 2 ;;
        --no-motif) RUN_MOTIF="no"; shift ;;
        -g|--ctk-group) CTK_GROUPS_FILE="$2"; shift 2 ;;
        --sample) SAMPLE_SIZE="$2"; shift 2 ;;
        -b|--barcode) BARCODE_FILE="$2"; DEMUX="yes"; shift 2 ;;
        -d|--input-dir) INPUT_DIR="$2"; shift 2 ;;
        --mismatches) MISMATCHES="$2"; shift 2 ;;
        --no-dedup) DEDUP_MODE="false"; shift ;;
        --dedup) DEDUP_MODE="true"; shift ;; # Keep for compat
        --skip-ncrna) SKIP_NCRNA="true"; shift ;;
        --notification) NOTIFY_MODE="true"; shift ;;
        --child) CHILD_MODE="true"; shift ;;
        --wizard|--advanced)
            WIZARD_MODE="true"
            shift
            ;;
        --verbose) VERBOSE="true"; shift ;;
        -h|--help) show_usage; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# ------------------------------------------------------------------
# Log File Initialization (Post-Parsing)
# ------------------------------------------------------------------
# If CHILD_MODE, we suppress main log creation to avoid spamming root dir.
if [[ "$CHILD_MODE" == "true" ]]; then
    LOG_FILE="/dev/null"
else
    # Use absolute path so redirects work correctly even if we cd later
# Check for existing config in current directory (e.g. from parent process or previous run)
if [[ -f "analysis_config.env" ]]; then
    # We only source if we are NOT running the wizard right now (child process)
    # OR if we want to load defaults. 
    # Logic: If --wizard, run wizard (which overwrites config). If not --wizard, try to load config.
    if [[ "$WIZARD_MODE" == "false" ]]; then
        source "analysis_config.env"
    fi
fi

if [[ "$WIZARD_MODE" == "true" ]]; then
    # Run the comprehensive wizard
    run_wizard_clipittyclip
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
    
    # Apply wizard settings to main script variables
    if [[ "$WIZ_MODE" == "pooled" ]]; then
        INPUT_FILE="$WIZ_INPUT_FILE"
        BARCODE_FILE="$WIZ_BARCODE_FILE"
    elif [[ "$WIZ_MODE" == "single" ]]; then
        INPUT_FILE="$WIZ_INPUT_FILE"
    elif [[ "$WIZ_MODE" == "directory" ]]; then
        INPUT_DIR="$WIZ_INPUT_DIR"
    fi
    
    GENOME_INDEX="$WIZ_GENOME_INDEX"
    ALIGNER="$WIZ_ALIGNER"
    THREADS="$WIZ_THREADS"
    UMI_LEN="$WIZ_UMI_LEN"
    ADAPTER="$WIZ_ADAPTER"
    [[ "$WIZ_CIMS" == "y" ]] && RUN_CIMS="true"
    [[ "$WIZ_CITS" == "y" ]] && RUN_CITS="true"
    PEAK_DIST="$WIZ_PEAK_DIST"
    PEAK_SIZE="$WIZ_PEAK_SIZE"
    FRAG_LEN="$WIZ_FRAG_LEN"
    ADV_HOMER_ARGS="$WIZ_HOMER_ARGS"
fi

# Log file setup
LOG_FILE="$(pwd)/CLIPittyClip_$(date +%Y%m%d_%H%M%S).log"
    touch "${LOG_FILE}"
fi

# Validation: Need either -i or -d (but not both) - skip if wizard already set them
if [[ -n "$INPUT_DIR" ]] && [[ -n "$INPUT_FILE" ]]; then
    log_error "Cannot use both -i and -d. Choose one input mode."
    show_usage
    exit 1
fi

if [[ -z "$INPUT_FILE" ]] && [[ -z "$INPUT_DIR" ]] && [[ "$WIZARD_MODE" == "false" ]]; then
    log_error "Missing required input (-i or -d). Provide a single file or input directory."
    show_usage
    exit 1
fi

if [[ -z "$GENOME_INDEX" ]]; then
    log_error "Missing required argument (-x genome index)."
    show_usage
    exit 1
fi

# Resolve absolute paths
if [[ -n "$INPUT_FILE" ]]; then
    check_file "$INPUT_FILE" || exit 1
    INPUT_FILE="$(cd "$(dirname "$INPUT_FILE")" && pwd)/$(basename "$INPUT_FILE")"
fi

if [[ -n "$INPUT_DIR" ]]; then
    if [[ ! -d "$INPUT_DIR" ]]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
    INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
    # Check for FASTQ files
    FASTQ_COUNT=$(ls "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz 2>/dev/null | wc -l)
    if [[ "$FASTQ_COUNT" -eq 0 ]]; then
        log_error "No .fastq.gz or .fq.gz files found in: $INPUT_DIR"
        exit 1
    fi
fi

GENOME_INDEX="$(cd "$GENOME_INDEX" && pwd)"

# Validate groups file (requires --cims or --cits)
if [[ -n "$CTK_GROUPS_FILE" ]]; then
    if [[ "$RUN_CIMS" != "true" ]] && [[ "$RUN_CITS" != "true" ]]; then
        log_error "-g/--groups requires --cims and/or --cits"
        show_usage
        exit 1
    fi
    if [[ ! -f "$CTK_GROUPS_FILE" ]]; then
        log_error "Groups file not found: $CTK_GROUPS_FILE"
        exit 1
    fi
    CTK_GROUPS_FILE="$(cd "$(dirname "$CTK_GROUPS_FILE")" && pwd)/$(basename "$CTK_GROUPS_FILE")"
    log_info "Groups file: $CTK_GROUPS_FILE"
fi

# ------------------------------------------------------------------
# Log Configuration (Standard or Advanced Results)
# ------------------------------------------------------------------
log_info "------------------------------------------------------------------"
log_info "Configuration Summary"
log_info "------------------------------------------------------------------"

# Define defaults for display (matching lib/modules.sh architecture)
DEF_FASTP="--length_required 16 --average_qual 30"
DEF_STAR="--outFilterMultimapNmax 10 --outFilterMismatchNmax ${MISMATCHES:-2} --alignEndsType EndToEnd --outSAMattributes ... MD"
DEF_BT2="--md --end-to-end (Standard Sensitivity)"
DEF_HOMER="-style factor -L 2 -localSize 10000 -minDist ${PEAK_DIST:-50}"

if [[ "$ADVANCED_MODE" == "true" ]]; then
    log_info "Mode:           ADVANCED (Using Overrides)"
    log_info "Fastp Defaults: $DEF_FASTP"
    log_info "Fastp Added:    ${ADV_FASTP_ARGS:-(None)}"
    
    log_info "Aligner:        $ALIGNER"
    if [[ "$ALIGNER" == "star" ]]; then log_info "Aligner Def:    $DEF_STAR"; else log_info "Aligner Def:    $DEF_BT2"; fi
    log_info "Aligner Added:  ${ADV_ALIGNER_ARGS:-(None)}"
    
    log_info "Homer Defaults: $DEF_HOMER"
    log_info "Homer Added:    ${ADV_HOMER_ARGS:-(None)}"
else
    log_info "Mode:           STANDARD"
    log_info "Fastp Config:   $DEF_FASTP"
    log_info "Aligner:        $ALIGNER"
    if [[ "$ALIGNER" == "star" ]]; then log_info "Aligner Config: $DEF_STAR"; else log_info "Aligner Config: $DEF_BT2"; fi
    log_info "Homer Config:   $DEF_HOMER"
fi
log_info "------------------------------------------------------------------"

# Dependencies Check
check_dependency fastp
check_dependency fastp
if [[ "$ALIGNER" == "bowtie2" ]]; then
    check_dependency bowtie2
else
    check_dependency STAR
fi
check_dependency bedtools
check_dependency samtools
check_dependency tag2collapse.pl
check_dependency seqkit

# Validate STAR Index
# Validate Index
if [[ "$ALIGNER" == "bowtie2" ]]; then
    check_bowtie_index "$GENOME_INDEX" || exit 1
else
    check_star_index "$GENOME_INDEX" || exit 1
fi

# Increase file limit for STAR sorting
ulimit -n 2048 2>/dev/null || true

# ------------------------------------------------------------------
# Display Clean Banner (If not Child)
# ------------------------------------------------------------------
if [[ "$CHILD_MODE" != "true" ]]; then
    # Print Clean Banner to Console
    console_msg "********************************************************************************"
    console_msg "CLIPittyClip Standard Pipeline v3.0"
    console_msg "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    if [[ -n "$INPUT_DIR" ]]; then
        console_msg "Input Directory: $INPUT_DIR"
    else
        console_msg "Input File: $(basename "$INPUT_FILE")"
    fi
    console_msg "Index: $GENOME_INDEX"
    console_msg "Threads: $THREADS | UMI: ${UMI_LEN}bp | Adapter: L32"
    console_msg "********************************************************************************"
    
    # Log Start (File Only)
    {
        echo "================================================================================"
        echo "CLIPittyClip Analysis Run Started: $(date)"
        if [[ -n "$INPUT_DIR" ]]; then
            echo "Input Directory: $INPUT_DIR"
        else
            echo "Input File: $INPUT_FILE"
        fi
        echo "Command Line: $0 $@"
        echo "--------------------------------------------------------------------------------"
    } >> "$LOG_FILE"
fi


# ------------------------------------------------------------------
if [[ "$SAMPLE_SIZE" -gt 0 ]] && [[ -n "$INPUT_FILE" ]]; then
    if [[ "$CHILD_MODE" != "true" ]]; then
        console_msg "\n[TEST DRIVE] Sampling first $SAMPLE_SIZE reads..."
    fi
    log_info "Test Drive Mode: Sampling first $SAMPLE_SIZE reads from input..."
    
    # Define sampled filename
    SAMPLED_INPUT="${INPUT_FILE%.fastq.gz}_sampled_${SAMPLE_SIZE}.fastq.gz"
    if [[ "$INPUT_FILE" == *.fq.gz ]]; then
        SAMPLED_INPUT="${INPUT_FILE%.fq.gz}_sampled_${SAMPLE_SIZE}.fastq.gz"
    fi
    SAMPLED_INPUT="$(basename "$SAMPLED_INPUT")" # Keep in CWD
    
    # Calculate lines: 4 lines per read
    LINES=$((SAMPLE_SIZE * 4))
    
    # Stream process: gunzip -> head -> gzip
    # This is efficient and avoids unzipping the whole file
    gzip -cd "$INPUT_FILE" | head -n "$LINES" | gzip > "$SAMPLED_INPUT"
    
    if [[ $? -eq 0 && -s "$SAMPLED_INPUT" ]]; then
        if [[ "$CHILD_MODE" != "true" ]]; then
            console_msg "  > Done. Created: $SAMPLED_INPUT"
        fi
        log_info "Sampling complete. Created: $SAMPLED_INPUT"
        # SWITCH INPUT to the small file for the rest of the pipeline
        INPUT_FILE="$(pwd)/$SAMPLED_INPUT"
        # Reset SAMPLE_SIZE to 0 so we don't try to resample inside child processes (redundant)
        SAMPLE_SIZE=0 
    else
        log_error "Sampling failed."
        exit 1
    fi
fi

# ------------------------------------------------------------------
# Main Pipeline
# ------------------------------------------------------------------

# 0a. Directory Input Mode (Pre-demultiplexed FASTQs)
if [[ -n "$INPUT_DIR" ]]; then
    console_msg "\n[DIRECTORY INPUT MODE]"
    console_msg "  > Input Directory: $INPUT_DIR"
    
    # Count files
    SAMPLE_FILES=("$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz)
    # Filter to only existing files
    SAMPLE_FILES=($(ls "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz 2>/dev/null))
    total_samples=${#SAMPLE_FILES[@]}
    
    console_msg "  > Found $total_samples sample files"
    
    # Build extra flags for child processes
    EXTRA_FLAGS=""
    # Only pass CIMS/CITS to children if NOT using groups file (group CTK runs after batch)
    if [[ -z "$CTK_GROUPS_FILE" ]]; then
        if [[ "$RUN_CIMS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --cims"; fi
        if [[ "$RUN_CITS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --cits"; fi
    fi
    if [[ "$VERBOSE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --verbose"; fi
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS -k"; fi
    if [[ "$SAMPLE_SIZE" -gt 0 ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --sample $SAMPLE_SIZE"; fi
    EXTRA_FLAGS="$EXTRA_FLAGS --aligner $ALIGNER"
    EXTRA_FLAGS="$EXTRA_FLAGS --no-dedup"  # Assume pre-processed, skip dedup
    EXTRA_FLAGS="$EXTRA_FLAGS --child"
    
    console_msg "\n[BATCH ANALYSIS]"
    
    current_sample=0
    for sample in "${SAMPLE_FILES[@]}"; do
        if [ -f "$sample" ]; then
            ((current_sample++))
            sample_name=$(basename "$sample")
            sample_name="${sample_name%.fastq.gz}"
            sample_name="${sample_name%.fq.gz}"
            
            printf "  %2d/%-2d %-20s : " "$current_sample" "$total_samples" "$sample_name"
            
            echo "[BATCH] Launching analysis for sample: $sample_name" >> "$LOG_FILE"
            
            self_script="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
            
            cmd="bash $self_script \
                -i $sample \
                -x $GENOME_INDEX \
                -t $THREADS \
                -u $UMI_LEN \
                -a $ADAPTER_3 \
                -o ${sample_name} \
                $EXTRA_FLAGS"
            
            $cmd
            
            if [ $? -eq 0 ]; then
                update_status_done
                echo "[BATCH] Sample $sample_name analysis complete." >> "$LOG_FILE"
            else
                echo -e "${RED}FAILED${NC}"
                log_error ">>> Sample $sample_name analysis FAILED."
            fi
        fi
    done
    
    # Aggregation - same as demux path
    INPUT_BASENAME=$(basename "$INPUT_DIR")
    OUTPUT_ROOT="${INPUT_BASENAME}_output"
    
    DIR_BAM="1_BAM"
    DIR_BED="2_COLLAPSED_BED"
    DIR_BG="3_BEDGRAPH"
    DIR_PEAKS="4_PEAKS"
    # OTHERS folder number depends on whether CTK analysis is enabled
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        DIR_OTHERS="6_OTHERS"
    else
        DIR_OTHERS="5_OTHERS"
    fi
    DIR_REPORTS="REPORTS"
    DIR_PEAK_LOGS="REPORTS/PEAK"
    DIR_IND_PEAK_LOGS="REPORTS/PEAK/INDIVIDUAL_SAMPLES"
    
    # Create base directories (CTK folders created conditionally during aggregation)
    mkdir -p "$OUTPUT_ROOT/$DIR_BAM" \
             "$OUTPUT_ROOT/$DIR_BED" \
             "$OUTPUT_ROOT/$DIR_BG" \
             "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS" \
             "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS" \
             "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/STAR_LOGS/DETAILED_LOGS_CAN_BE_DELETED" \
             "$OUTPUT_ROOT/$DIR_PEAK_LOGS" \
             "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS"
    
    # Run Group-based CTK Analysis (if groups file provided)
    if [[ -n "$CTK_GROUPS_FILE" ]]; then
        run_group_ctk_analysis "$CTK_GROUPS_FILE" "$OUTPUT_ROOT" "$GENOME_INDEX" \
            "$CIMS_ITERATIONS" "$CIMS_FDR" "$CITS_PVALUE" "$CITS_GAP" \
            "$MOTIF_FLANK" "$RUN_MOTIF" "$RUN_CIMS" "$RUN_CITS"
    fi
    
    console_msg "\n[AGGREGATION]"
    
    # Collect outputs from each sample's output directory
    for sample in "${SAMPLE_FILES[@]}"; do
        sample_name=$(basename "$sample")
        sample_name="${sample_name%.fastq.gz}"
        sample_name="${sample_name%.fq.gz}"
        sample_out="${sample_name}_analysis"
        
        if [[ -d "$sample_out" ]]; then
            console_msg "  > Collecting: $sample_name"
            
            # Move BAMs
            mv "$sample_out"/1_BAM/*.bam* "$OUTPUT_ROOT/$DIR_BAM/" 2>/dev/null
            
            # Move CIMS/CITS (if they exist)
            if [[ -d "$sample_out/5_CIMS" ]]; then
                 mv "$sample_out"/5_CIMS/* "$OUTPUT_ROOT/$DIR_CIMS/" 2>/dev/null
            else
                 # If child process produced flat files (likely)
                 mv "$sample_out"/*_CIMS.txt "$OUTPUT_ROOT/$DIR_CIMS/" 2>/dev/null
            fi
            
            if [[ -d "$sample_out/6_CITS" ]]; then
                 mv "$sample_out"/6_CITS/* "$OUTPUT_ROOT/$DIR_CITS/" 2>/dev/null
            else
                 mv "$sample_out"/*_CITS.bed "$OUTPUT_ROOT/$DIR_CITS/" 2>/dev/null
            fi
            # Move BEDs
            mv "$sample_out"/2_COLLAPSED_BED/*.bed "$OUTPUT_ROOT/$DIR_BED/" 2>/dev/null
            # Move Bedgraphs
            mv "$sample_out"/3_BEDGRAPH/*.bedgraph "$OUTPUT_ROOT/$DIR_BG/" 2>/dev/null
            cp "$sample_out"/3_BEDGRAPH/chrom.sizes "$OUTPUT_ROOT/$DIR_BG/" 2>/dev/null
            # Move Peak folders
            mv "$sample_out"/4_PEAKS/* "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/${sample_name}_peaks/" 2>/dev/null || true
            # Move Reports
            mv "$sample_out"/REPORTS/FASTP_REPORT/* "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
            mv "$sample_out"/REPORTS/*_pipeline.log "$OUTPUT_ROOT/$DIR_REPORTS/" 2>/dev/null
            
            # Move ncRNA Mapping outputs
            if [[ -d "$sample_out/OTHERS/ncRNA_Mapping" ]]; then
                mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping"
                mv "$sample_out"/OTHERS/ncRNA_Mapping/*_ncrna_stats.txt "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                mv "$sample_out"/OTHERS/ncRNA_Mapping/*_ncrna.bam* "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
            fi
            
            # Cleanup sample dir (keep if -k)
            if [[ "$KEEP_INTERMEDIATE" != "yes" ]]; then
                rm -rf "$sample_out"
            fi
        fi
    done
    
    # ncRNA Filtering Summary
    console_msg "\n[ncRNA FILTERING SUMMARY]"
    printf "  %-25s %-15s %-15s %s\n" "Sample" "ncRNA Reads" "Total Reads" "% Filtered"
    console_msg "  ----------------------------------------------------------------"
    
    for sample in "${SAMPLE_FILES[@]}"; do
        sample_name=$(basename "$sample")
        sample_name="${sample_name%.fastq.gz}"
        sample_name="${sample_name%.fq.gz}"
        sample_out="${sample_name}_analysis"
        ncrna_stats="${sample_out}/OTHERS/ncRNA_Mapping/${sample_name}_ncrna_stats.txt"
        
        if [[ -f "$ncrna_stats" ]]; then
            align_rate=$(grep "overall alignment rate" "$ncrna_stats" | grep -oE "[0-9]+\.[0-9]+%" || echo "N/A")
            total=$(grep "reads; of these:" "$ncrna_stats" | grep -oE "^[0-9]+" || echo "N/A")
            aligned=$(grep "aligned exactly 1 time" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
            multi=$(grep "aligned >1 times" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
            ncrna=$((aligned + multi))
            printf "  %-25s %-15s %-15s %s\n" "$sample_name" "$ncrna" "$total" "$align_rate"
        else
            printf "  %-25s %-15s %-15s %s\n" "$sample_name" "-" "-" "SKIPPED"
        fi
    done
    console_msg "  ----------------------------------------------------------------"
    
    # Combined Peak Calling
    console_msg "\n[COMBINED PEAK CALLING]"
    BED_DIR="$OUTPUT_ROOT/$DIR_BED"
    BED_SYMLINK="BED"
    
    if [[ -d "$BED_SYMLINK" ]]; then rm -rf "$BED_SYMLINK"; fi
    ln -s "$BED_DIR" "$BED_SYMLINK"
    
    PEAK_LOG="$OUTPUT_ROOT/$DIR_PEAK_LOGS/${INPUT_BASENAME}_PeakCalling.log"
    "$SCRIPT_DIR/PEAKittyPeak.sh" -n "Combined" -p "$PEAK_DIST" > "$PEAK_LOG" 2>&1
    
    rm -f "$BED_SYMLINK"
    
    # Move combined results
    if [[ -d "PEAKS" ]]; then
        mv PEAKS/* "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/" 2>/dev/null
        rm -rf PEAKS
    fi
    
    console_msg "  > Combined peaks: $OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/"
    
    # Final Summary
    PIPELINE_END=$(date +%s)
    DURATION=$((PIPELINE_END - PIPELINE_START))
    H=$((DURATION/3600))
    M=$(( (DURATION%3600)/60 ))
    S=$((DURATION%60))
    
    console_msg "\n[COMPLETE]"
    console_msg "  > Duration: ${H}h ${M}m ${S}s"
    console_msg "  > Output: $OUTPUT_ROOT/"
    
    # Save console summary log at output folder level
    CONSOLE_LOG="${OUTPUT_ROOT}_summary.log"
    {
        echo "========================================"
        echo "CLIPittyClip Analysis Summary"
        echo "========================================"
        echo "Date: $(date)"
        echo "Input Directory: ${INPUT_DIR}"
        echo "Samples: $total_samples"
        echo "Duration: ${H}h ${M}m ${S}s"
        echo ""
        echo "Output: $OUTPUT_ROOT/"
        echo "  - BAM files:      $OUTPUT_ROOT/$DIR_BAM/"
        echo "  - BED files:      $OUTPUT_ROOT/$DIR_BED/"
        echo "  - Peaks:          $OUTPUT_ROOT/$DIR_PEAKS/"
        echo "  - Detailed logs:  $OUTPUT_ROOT/$DIR_REPORTS/"
        echo "========================================"
    } > "$CONSOLE_LOG"
    console_msg "  > Summary log: $CONSOLE_LOG"
    
    send_notification "CLIPittyClip" "Directory batch analysis complete: $total_samples samples"
    
    exit 0
fi

# 0b. Demultiplexing (Recursive Branch)
if [[ "$DEMUX" == "yes" ]]; then
    
    # 1a. Deduplication (if enabled)
    if [[ "$DEDUP_MODE" == "true" ]]; then
        console_msg "\n[DEDUPLICATING]"
        print_section_item "Deduplicating Pooled Reads (SeqKit)"
        
        dedup_temp="pooled_dedup.fastq.gz"
        seqkit rmdup -s -o "$dedup_temp" "$INPUT_FILE" 2>> "${LOG_FILE}"
        
        if [[ $? -eq 0 && -s "$dedup_temp" ]]; then
            WORK_INPUT="$dedup_temp"
            print_section_item "Deduplicating Complete"
        else
            log_warning "Deduplication failed. Using original file."
            rm -f "$dedup_temp"
            WORK_INPUT="$INPUT_FILE"
        fi
    else
        WORK_INPUT="$INPUT_FILE"
    fi
    
    # 1b. Demultiplexing
    console_msg "\n[DEMULTIPLEXING]"
    print_section_item "Barcodes: $(basename "$BARCODE_FILE")"
    print_section_item "Mismatches Allowed: $MISMATCHES"
    
    # Run demultiplexing (with dedup disabled since we already did it)
    run_demultiplexing "$WORK_INPUT" "$BARCODE_FILE" "$SAMPLE_SIZE" "false"
    
    print_section_item "Demultiplexing Complete"
    
    # Cleanup dedup temp file
    if [[ -f "pooled_dedup.fastq.gz" ]]; then
        rm -f "pooled_dedup.fastq.gz"
        log_info "Cleaned up pooled dedup temp file."
    fi
    
    send_notification "CLIPittyClip" "Demultiplexing complete for $(basename "$INPUT_FILE")"
    
    echo "" # Newline after progress

    # 1.5. Demux Summary Table
    console_msg "\n[DEMUX SUMMARY]"
    printf "  %-20s %-12s %s\n" "Sample" "Reads" "% of Total"
    console_msg "  ------------------------------------------------"
    
    total_reads=0
    # Calculate totals first
    for sample in demux_fastq/*.fastq.gz; do
        if [ -f "$sample" ]; then
            # Fast counting using zcat/wc is slow for huge files.
            # But fastp json output from demux might be better?
            # Cutadapt demux doesn't give easy per-sample JSON unless parsed.
            # We will use wc -l / 4. 
            # Ideally this should be optimized but for now is robust.
            # Actually, `fastq` is 4 lines per read.
            lines=$(gzip -dc "$sample" | wc -l)
            count=$((lines / 4))
            total_reads=$((total_reads + count))
        fi
    done

    # Print table
    for sample in demux_fastq/*.fastq.gz; do
        if [ -f "$sample" ]; then
            sample_name=$(basename "$sample" .fastq.gz)
            lines=$(gzip -dc "$sample" | wc -l)
            count=$((lines / 4))
            
            if [ "$total_reads" -gt 0 ]; then
                pct=$(awk "BEGIN {printf \"%.1f\", ($count / $total_reads) * 100}")
            else
                pct="0.0"
            fi
            
            printf "  %-20s %-12s %s%%\n" "$sample_name" "$count" "$pct"
            # Log the stats too
            echo "[STATS] $sample_name: $count reads ($pct%)" >> "$LOG_FILE"
        fi
    done
    console_msg "  ------------------------------------------------"
    console_msg "  Total Processed: $total_reads reads"
    
    # Check for unknown samples and notify
    if [[ -f "demux_fastq/unknown.fastq.gz" ]]; then
        console_msg "  ${YELLOW}Note: 'unknown' samples will not be included in batch analysis.${NC}"
    fi

    # 2. Iterate and Recurse
    console_msg "\n[BATCH ANALYSIS]"
    
    # Check if we need to pass CIMS/CITS flags
    EXTRA_FLAGS=""
    # Only pass CIMS/CITS to children if NOT using groups file (group CTK runs after batch)
    if [[ -z "$CTK_GROUPS_FILE" ]]; then
        if [[ "$RUN_CIMS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --cims"; fi
        if [[ "$RUN_CITS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --cits"; fi
    fi
    if [[ "$VERBOSE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --verbose"; fi
    if [[ "$SAMPLE_SIZE" -gt 0 ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --sample $SAMPLE_SIZE"; fi
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS -k"; fi
    
    # Pass Aligner choice
    EXTRA_FLAGS="$EXTRA_FLAGS --aligner $ALIGNER"
    
    # We already deduped the pooled file before demux (if DEDUP_MODE was true).
    # So we explicitly tell child processes NOT to dedup again to save time.
    # Note: If user passed --no-dedup, DEDUP_MODE is false, so we don't care (default child is true though).
    # Wait, default child is true. So we MUST pass --no-dedup regardless of parent mode?
    # Yes, because the split files are coming from a potentially deduped source.
    # Even if they weren't deduped (parent --no-dedup), deduping individual split files is less efficient 
    # than pool dedup, but maybe user wants it? 
    # Logic: Parent DEDUP_MODE controls the whole workflow.
    # If Parent DEDUP=true -> Pool Dedup -> Split -> Child (should not dedup).
    # If Parent DEDUP=false -> No Pool Dedup -> Split -> Child (should not dedup? or follow parent?)
    # Generally, pass --no-dedup to child to avoid double-dedup.
    
    EXTRA_FLAGS="$EXTRA_FLAGS --no-dedup"
    
    # Pass --child to suppress header in sub-calls
    EXTRA_FLAGS="$EXTRA_FLAGS --child"

    # Get sample count for progress (excluding unknown)
    total_samples=$(ls demux_fastq/*.fastq.gz 2>/dev/null | grep -v "/unknown.fastq.gz" | wc -l | tr -d ' ')
    current_sample=0

    # Loop through all generated files (skip unknown)
    for sample in demux_fastq/*.fastq.gz; do
        if [ -f "$sample" ]; then
            sample_name=$(basename "$sample" .fastq.gz)
            
            # Skip unknown samples
            if [[ "$sample_name" == "unknown" ]]; then
                continue
            fi
            
            ((current_sample++))
            
            # Print the leader for this sample. 
            # We use printf without newline, or just let update_status handle the rest?
            # The Child process will overwrite the same line.
            # So we print the static prefix here:
            printf "  %2d/%-2d %-20s : " "$current_sample" "$total_samples" "$sample_name"
            
            # Log launch
            echo "[BATCH] Launching analysis for sample: $sample_name" >> "$LOG_FILE"
            
            # Construct command for recursive call
            self_script="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
            
            # We must NOT pass -b again
            cmd="bash $self_script \
                -i $sample \
                -x $GENOME_INDEX \
                -t $THREADS \
                -u $UMI_LEN \
                -a $ADAPTER_3 \
                -o ${sample_name} \
                $EXTRA_FLAGS"
            
            # Run the command. 
            # In CHILD_MODE, the script will produce `update_status` outputs (echo -ne ... \r).
            # This should appear on the console right after our printf above.
            $cmd
            
            # When child finishes, we need to print "Done." and a newline to finalize the line.
            if [ $? -eq 0 ]; then
                update_status_done
                echo "[BATCH] Sample $sample_name analysis complete." >> "$LOG_FILE"
            else
                echo -e "${RED}FAILED${NC}"
                log_error ">>> Sample $sample_name analysis FAILED."
            fi
        fi
    done
    
    # 3. Aggregation and Filing
    # Define Final Output Structure
    INPUT_BASENAME=$(basename "$INPUT_FILE" .fastq.gz) 
    if [[ "$INPUT_BASENAME" == *"_sampled_"* ]]; then
        # If sampled was used, try to get original name or keep it
        # Actually safer to use the original argument input name if possible, 
        # but here we rely on what we have.
        # Let's just use a clean name.
        INPUT_BASENAME="${INPUT_BASENAME%_sampled_*}"
    fi
    OUTPUT_ROOT="${INPUT_BASENAME}_output"
    
    DIR_DEMUX="0_DEMUX_FASTQ"
    DIR_BAM="1_BAM"
    DIR_BED="2_COLLAPSED_BED"
    DIR_BG="3_BEDGRAPH"
    DIR_PEAKS="4_PEAKS"
    # OTHERS folder number depends on whether CTK analysis is enabled
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        DIR_OTHERS="6_OTHERS"
    else
        DIR_OTHERS="5_OTHERS"
    fi
    DIR_REPORTS="REPORTS"
    DIR_PEAK_LOGS="REPORTS/PEAK"
    DIR_IND_PEAK_LOGS="REPORTS/PEAK/INDIVIDUAL_SAMPLES"
    
    # Create Central Output Directories (CTK folders created conditionally during aggregation)
    mkdir -p "$OUTPUT_ROOT/$DIR_DEMUX" \
             "$OUTPUT_ROOT/$DIR_BAM" \
             "$OUTPUT_ROOT/$DIR_BED" \
             "$OUTPUT_ROOT/$DIR_BG" \
             "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS" \
             "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS" \
             "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/STAR_LOGS/DETAILED_LOGS_CAN_BE_DELETED" \
             "$OUTPUT_ROOT/$DIR_PEAK_LOGS" \
             "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS"
    
    # Run Group-based CTK Analysis (if groups file provided)
    if [[ -n "$CTK_GROUPS_FILE" ]]; then
        run_group_ctk_analysis "$CTK_GROUPS_FILE" "$OUTPUT_ROOT" "$GENOME_INDEX" \
            "$CIMS_ITERATIONS" "$CIMS_FDR" "$CITS_PVALUE" "$CITS_GAP" \
            "$MOTIF_FLANK" "$RUN_MOTIF" "$RUN_CIMS" "$RUN_CITS"
    fi
             
    console_msg "\n[AGGREGATION]"
    console_msg "  > Collecting files into $OUTPUT_ROOT/..."
    
    # Move Barcode file if exists? Users usually keep input.
    # We leave inputs alone.
    
    # Collect files from analysis directories
    count=0
    for sample in demux_fastq/*.fastq.gz; do
        if [ -f "$sample" ]; then
            sample_name=$(basename "$sample" .fastq.gz)
            
            # Skip unknown samples
            [[ "$sample_name" == "unknown" ]] && continue
            
            analysis_dir="${sample_name}_analysis"
            
            if [[ -d "$analysis_dir" ]]; then
                # log_info "Organizing output for $sample_name..."
                
                # 0. DEMUX FASTQ (Original is in demux_fastq, we can move it)
                 # Wait, we loop over demux_fastq files.
                 # Let's move them separately at the end.
                
                # 1. BAM
                bam_file=$(find "$analysis_dir" -name "*mapped.Aligned.sortedByCoord.out.bam")
                if [[ -n "$bam_file" ]]; then
                    cp "$bam_file" "$OUTPUT_ROOT/$DIR_BAM/${sample_name}.bam"
                    cp "${bam_file}.bai" "$OUTPUT_ROOT/$DIR_BAM/${sample_name}.bam.bai" 2>/dev/null
                fi
                
                # 2. BED (Collapsed)
                bed_file=$(find "$analysis_dir" -name "*_collapsed.bed")
                if [[ -n "$bed_file" ]]; then
                    cp "$bed_file" "$OUTPUT_ROOT/$DIR_BED/${sample_name}.bed"
                    ((count++))
                fi
                
                # 3. Bedgraph & Chrom Sizes
                bg_pos=$(find "$analysis_dir" -name "*_pos.bedgraph")
                bg_neg=$(find "$analysis_dir" -name "*_neg.bedgraph")
                if [[ -n "$bg_pos" ]]; then cp "$bg_pos" "$OUTPUT_ROOT/$DIR_BG/${sample_name}_pos.bedgraph"; fi
                if [[ -n "$bg_neg" ]]; then cp "$bg_neg" "$OUTPUT_ROOT/$DIR_BG/${sample_name}_neg.bedgraph"; fi
                
                # Grab chrom.sizes from ONE sample (overwrite is fine)
                 # It might be in analysis dir if we extracted it, or just use user provided if needed
                 # modules.sh extracts it to 'chrom.sizes' in CWD if not provided?
                 # Actually modules.sh:246 writes to "$chrom_sizes".
                 extracted_chroms="$analysis_dir/chrom.sizes"
                 if [[ -f "$extracted_chroms" ]]; then
                    cp "$extracted_chroms" "$OUTPUT_ROOT/$DIR_BG/chrom.sizes"
                 fi

                # 4. Peaks (Folder)
                peak_dir="${analysis_dir}/${sample_name}_peaks"
                if [[ -d "$peak_dir" ]]; then
                    cp -r "$peak_dir" "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/"
                fi

                # 5. CTK Analysis (CIMS/CITS) - check for new folder names
                for ctk_folder in "CTK_Analysis" "CIMS_Analysis" "CITS_Analysis"; do
                    if [[ -d "$analysis_dir/$ctk_folder" ]]; then
                        mkdir -p "$OUTPUT_ROOT/5_${ctk_folder}"
                        cp -r "$analysis_dir/$ctk_folder/"* "$OUTPUT_ROOT/5_${ctk_folder}/" 2>/dev/null
                        break
                    fi
                done
                    
                # 4.1. Individual Peak Log (Created by updated modules.sh)
                # It was created as ${peak_dir}_homer.log
                peak_log="${peak_dir}_homer.log"
                if [[ -f "$peak_log" ]]; then
                    cp "$peak_log" "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS/${sample_name}_PeakCalling.log"
                fi

                # 5. Reports & Logs
                # FASTP
                cp "$analysis_dir"/*_fastp.html "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
                cp "$analysis_dir"/*_fastp.json "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
                
                # STAR
                cp "$analysis_dir"/*Log.final.out "$OUTPUT_ROOT/$DIR_REPORTS/STAR_LOGS/" 2>/dev/null
                cp "$analysis_dir"/*Log.out "$OUTPUT_ROOT/$DIR_REPORTS/STAR_LOGS/DETAILED_LOGS_CAN_BE_DELETED/" 2>/dev/null
                cp "$analysis_dir"/*Log.progress.out "$OUTPUT_ROOT/$DIR_REPORTS/STAR_LOGS/DETAILED_LOGS_CAN_BE_DELETED/" 2>/dev/null
                cp "$analysis_dir"/*SJ.out.tab "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT/" 2>/dev/null
                
                # ncRNA Mapping outputs
                if [[ -d "$analysis_dir/OTHERS/ncRNA_Mapping" ]]; then
                    mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping"
                    cp "$analysis_dir"/OTHERS/ncRNA_Mapping/*_ncrna_stats.txt "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                    cp "$analysis_dir"/OTHERS/ncRNA_Mapping/*_ncrna.bam* "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                fi
                
                # Pipeline Log (The child specific one)
                # It's named ${sample_name}_analysis.log inside the dir
                cp "$analysis_dir"/*.log "$OUTPUT_ROOT/$DIR_REPORTS/${sample_name}_detailed.log" 2>/dev/null
                
            else
                log_warning "Analysis directory not found for $sample_name"
            fi
        fi
    done
    
    # Move Demux Fastqs
    mv demux_fastq/*.fastq.gz "$OUTPUT_ROOT/$DIR_DEMUX/" 2>/dev/null
    rmdir demux_fastq 2>/dev/null
    
    if [[ "$count" -eq 0 ]]; then
        log_error "No collapsed BED files collected. Skipping Peak Calling."
        exit 1
    fi
    console_msg "  > Found $count BED files."
    log_info "Collected $count BED files."
    
    # Run PEAKittyPeak.sh on Aggregated Data
    # PEAK_SCRIPT requires BED/ directory. We named it 2_COLLAPSED_BED...
    # We should update PEAKittyPeak or symlink.
    # Updating PEAKittyPeak is cleaner, but let's just make a temporary symlink?
    # Or just tell PEAKittyPeak where to look?
    # PEAKittyPeak.sh is simple. Let's just create the "Combined Input" manually here?
    # No, keep modularity. 
    # Let's symlink:
    ln -s "2_COLLAPSED_BED" "$OUTPUT_ROOT/BED"
    
    PEAK_SCRIPT="$(dirname "$0")/PEAKittyPeak.sh"
    
    if [[ -x "$PEAK_SCRIPT" ]]; then
        console_msg "  > Running HOMER Peak Calling (Aggregated)..."
        
        # Run it inside the output root so it sees "BED" symlink and writes "COMBINED_PEAKS" there
        curr_dir=$(pwd)
        cd "$OUTPUT_ROOT" || exit 1
        
        # Define separate log for peak calling details (Specific Name)
        PEAK_LOG="$(pwd)/$DIR_PEAK_LOGS/${INPUT_BASENAME}_PeakCalling.log"
        console_msg "  > Detailed Peak Log: $DIR_PEAK_LOGS/$(basename "$PEAK_LOG")"

        # Call script (using absolute path or relative to old pwd)
        # Call script with -n COMBINED to create COMBINED_peaks folder
        "$PEAK_SCRIPT" -n COMBINED > "$PEAK_LOG" 2>&1
        
        # Remove symlink
        rm BED
        
        # Move generated `COMBINED_peaks` to strictly correct place if needed
        # PEAKittyPeak creates "{BASE_NAME}_peaks" folder in CWD ($OUTPUT_ROOT).
        # We want it in "$DIR_PEAKS/COMBINED_PEAKS".
        if [[ -d "COMBINED_peaks" ]]; then
            mv "COMBINED_peaks"/* "$DIR_PEAKS/COMBINED_PEAKS/" 2>/dev/null
            rmdir "COMBINED_peaks"
        fi
        
        cd "$curr_dir"

        if [ $? -eq 0 ]; then
             console_msg "  > Peak Calling Complete."
        else
             console_msg "  > ${RED}Peak Calling Failed (Check Log)${NC}"
        fi
        
    else
        log_error "PEAKittyPeak.sh not found."
    fi
    
    # ncRNA Filtering Summary (before cleanup so stats files still exist)
    console_msg "\n[ncRNA FILTERING SUMMARY]"
    printf "  %-25s %-15s %-15s %s\n" "Sample" "ncRNA Reads" "Total Reads" "% Filtered"
    console_msg "  ----------------------------------------------------------------"
    
    for analysis_dir in *_analysis; do
        if [[ -d "$analysis_dir" ]]; then
            sample_name="${analysis_dir%_analysis}"
            [[ "$sample_name" == "unknown" ]] && continue
            
            ncrna_stats="${analysis_dir}/OTHERS/ncRNA_Mapping/${sample_name}_ncrna_stats.txt"
            
            if [[ -f "$ncrna_stats" ]]; then
                align_rate=$(grep "overall alignment rate" "$ncrna_stats" | grep -oE "[0-9]+\.[0-9]+%" || echo "N/A")
                total=$(grep "reads; of these:" "$ncrna_stats" | grep -oE "^[0-9]+" || echo "N/A")
                aligned=$(grep "aligned exactly 1 time" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
                multi=$(grep "aligned >1 times" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
                ncrna=$((aligned + multi))
                printf "  %-25s %-15s %-15s %s\n" "$sample_name" "$ncrna" "$total" "$align_rate"
            else
                printf "  %-25s %-15s %-15s %s\n" "$sample_name" "-" "-" "SKIPPED"
            fi
        fi
    done
    console_msg "  ----------------------------------------------------------------"
    
    # Cleanup Analysis Folders and Temp Files
    console_msg "  > Cleaning up temporary analysis directories and files..."
    rm -rf *_analysis barcodes.fasta
    
    # Remove sampled input if it exists (pattern match)
    rm -f *_sampled_*.fastq.gz
    
    # Remove CTK temp directories (CITS.pl_* and CIMS.pl_*)
    rm -rf CITS.pl_* CIMS.pl_* 2>/dev/null
    
    # Remove CTK temp files (*.tmp from perl scripts)
    rm -f *.tmp 2>/dev/null

    # Also remove the main log if it was created in this dir (CLIPittyClip_*.log)
    # We want to MOVE it to REPORTS and rename it.
    if [[ -f "$LOG_FILE" ]]; then
         mv "$LOG_FILE" "$OUTPUT_ROOT/$DIR_REPORTS/${INPUT_BASENAME}_CLIPittyClip.log"
         LOG_FILE="$OUTPUT_ROOT/$DIR_REPORTS/${INPUT_BASENAME}_CLIPittyClip.log"
    fi
    
    # Output Summary
    console_msg "\n[OUTPUT]"
    console_msg "  All results saved to: $OUTPUT_ROOT/"
    console_msg "  ├── 0_DEMUX_FASTQ/"
    console_msg "  ├── 1_BAM/"
    console_msg "  ├── 2_COLLAPSED_BED/"
    console_msg "  ├── 3_BEDGRAPH/"
    console_msg "  ├── 4_PEAKS/"
    console_msg "  ├── 5_OTHERS/"
    console_msg "  └── REPORTS/"

    console_msg "\n[SUCCESS] Pipeline finished."
    
    # Calculate Duration
    PIPELINE_END=$(date +%s)
    DURATION=$((PIPELINE_END - PIPELINE_START))
    # Format HH:MM:SS
    H=$((DURATION/3600))
    M=$(( (DURATION%3600)/60 ))
    S=$((DURATION%60))
    
    console_msg "End Time: $(date '+%Y-%m-%d %H:%M:%S')"
    console_msg "Total Duration: ${H}h ${M}m ${S}s"
    
    # Save console summary log at output folder level
    CONSOLE_LOG="${OUTPUT_ROOT}_summary.log"
    {
        echo "========================================"
        echo "CLIPittyClip Analysis Summary"
        echo "========================================"
        echo "Date: $(date)"
        echo "Input: ${INPUT_FILE}"
        echo "Samples: $total_samples"
        echo "Duration: ${H}h ${M}m ${S}s"
        echo ""
        echo "Output: $OUTPUT_ROOT/"
        echo "  - BAM files:      $OUTPUT_ROOT/$DIR_BAM/"
        echo "  - BED files:      $OUTPUT_ROOT/$DIR_BED/"
        echo "  - Peaks:          $OUTPUT_ROOT/$DIR_PEAKS/"
        echo "  - Detailed logs:  $OUTPUT_ROOT/$DIR_REPORTS/"
        echo "========================================"
    } > "$CONSOLE_LOG"
    console_msg "  > Summary log: $CONSOLE_LOG"
    
    send_notification "CLIPittyClip" "Pipeline execution finished successfully. Duration: ${H}h ${M}m ${S}s"
    exit 0
fi

# Standard Pipeline (Single Sample)
# If CHILD_MODE, header is suppressed.
# Main Banner logic is now handled at the top of the script.
# We just log the start here for file consistency.
log_info "Analysis Started."

log_info "Input: $INPUT_FILE"
log_info "Genome Index: $GENOME_INDEX"
log_info "Threads: $THREADS"

# Set Basename
BASENAME=$(basename "$INPUT_FILE" .fastq.gz)
if [[ "$BASENAME" == "$INPUT_FILE" ]]; then
    BASENAME=$(basename "$INPUT_FILE" .fq.gz)
fi
if [[ -n "$EXP_ID" ]]; then
    BASENAME="$EXP_ID"
fi

# Directory Setup
mkdir -p "${BASENAME}_analysis"
cd "${BASENAME}_analysis" || exit 1
LOG_FILE="${BASENAME}_analysis.log" # redirect log to inside analysis dir
log_info "Working directory: $(pwd)"

# 1. Preprocessing
run_preprocessing "$INPUT_FILE" "$BASENAME" "$UMI_LEN" "$ADAPTER_3" "$THREADS" "$SAMPLE_SIZE" "$DEDUP_MODE"

# 1b. ncRNA Pre-filtering (if enabled and index exists)
CLEANED_FASTQ="${BASENAME}_cleaned.fastq.gz"
if [[ "$SKIP_NCRNA" == "false" ]]; then
    NCRNA_INDEX_DIR=$(check_ncrna_index "$GENOME_INDEX")
    if [[ -n "$NCRNA_INDEX_DIR" ]]; then
        NCRNA_OUTPUT_DIR="OTHERS/ncRNA_Mapping"
        NCRNA_UNMAPPED="${BASENAME}_ncrna_filtered.fastq.gz"
        run_ncrna_filter "$CLEANED_FASTQ" "$NCRNA_UNMAPPED" "$NCRNA_OUTPUT_DIR" "$NCRNA_INDEX_DIR" "$THREADS" "$BASENAME"
        # Use filtered reads for genome mapping
        CLEANED_FASTQ="$NCRNA_UNMAPPED"
    else
        log_warning "ncRNA index not found in $GENOME_INDEX or $GENOME_INDEX/ncRNA. Skipping ncRNA pre-filtering."
        log_warning "To build ncRNA index, see README.md for instructions."
    fi
fi

# 2. Mapping
MAPPED_PREFIX="${BASENAME}_mapped"
if [[ "$ALIGNER" == "bowtie2" ]]; then
    run_mapping_bowtie2 "$CLEANED_FASTQ" "$MAPPED_PREFIX" "$GENOME_INDEX" "$THREADS"
else
    # Default: STAR
    run_mapping_star "$CLEANED_FASTQ" "$MAPPED_PREFIX" "$GENOME_INDEX" "$THREADS" "$MISMATCH_MAX"
fi
# 3. PCR Collapse
# BAM is "${MAPPED_PREFIX}.Aligned.sortedByCoord.out.bam"
BAM_FILE="${MAPPED_PREFIX}.Aligned.sortedByCoord.out.bam"
COLLAPSED_BED="${BASENAME}_collapsed.bed"
MUTATION_FILE="${BASENAME}_mutations.txt"

# Unified preprocessing: Always use parseAlignment.pl for consistent output
# This generates both tags.bed and mutations.txt (future-proofing for CIMS/CITS)
log_info "Unified preprocessing: samtools calmd → parseAlignment.pl → tag2collapse.pl"
check_dependency parseAlignment.pl

run_parse_alignment "${BAM_FILE}" "${BASENAME}_parsed.bed" "${MUTATION_FILE}" "$GENOME_INDEX"

run_collapse_pcr "${BASENAME}_parsed.bed" "${COLLAPSED_BED}" "${UMI_LEN}"

# 4. Coverage Analysis (Bedgraph)
run_coverage "${COLLAPSED_BED}" "${BASENAME}" "$GENOME_INDEX/chrom.sizes" # Pass genome file if available

# 5. Peak Calling (Per-Sample - Optional if using Aggregation, but good for QC)
# 5. Peak Calling (Per-Sample - Optional if using Aggregation, but good for QC)
PEAK_DIR="${BASENAME}_peaks"
run_peak_calling "${COLLAPSED_BED}" "${PEAK_DIR}" "$PEAK_DIST" "$PEAK_SIZE" "$FRAG_LEN"

# 6. CIMS / CITS (CTK Analysis)
if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
    log_info "Running CTK CIMS/CITS Analysis..."
    
    # Determine output folder name based on flags
    if [[ "$RUN_CIMS" == "true" && "$RUN_CITS" == "true" ]]; then
        CTK_OUTPUT="CTK_Analysis"
    elif [[ "$RUN_CIMS" == "true" ]]; then
        CTK_OUTPUT="CIMS_Analysis"
    else
        CTK_OUTPUT="CITS_Analysis"
    fi
    
    # Find reference FASTA for motif analysis
    # Priority: 1) *genome*.fa, 2) *primary*.fa, 3) any .fa excluding *ncrna*
    ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \( -name "*genome*.fa" -o -name "*genome*.fasta" \) 2>/dev/null | head -n 1)
    if [[ -z "$ref_fasta" ]]; then
        ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \( -name "*primary*.fa" -o -name "*primary*.fasta" \) 2>/dev/null | head -n 1)
    fi
    if [[ -z "$ref_fasta" ]]; then
        # Fallback: any .fa/.fasta excluding ncrna
        ref_fasta=$(find "$GENOME_INDEX" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) ! -name "*ncrna*" 2>/dev/null | head -n 1)
    fi
    if [[ -z "$ref_fasta" ]]; then
        log_warning "Reference FASTA not found. Motif analysis may be skipped."
        ref_fasta=""
    else
        log_info "Using reference FASTA for motif analysis: $ref_fasta"
    fi
    
    mkdir -p "$CTK_OUTPUT"
    
    # Run CTK analysis using the STANDARD pipeline outputs (no duplicate preprocessing)
    run_ctk_analysis \
        "${COLLAPSED_BED}" \
        "${MUTATION_FILE}" \
        "$CTK_OUTPUT" \
        "$ref_fasta" \
        "${BASENAME}" \
        "$CIMS_ITERATIONS" \
        "$CIMS_FDR" \
        "$CITS_PVALUE" \
        "$CITS_GAP" \
        "$MOTIF_FLANK" \
        "$RUN_MOTIF" \
        "$RUN_CIMS" \
        "$RUN_CITS"
    
    log_info "CTK Analysis complete. Output: $CTK_OUTPUT"
fi

# Cleanup
if [[ "$KEEP_INTERMEDIATE" != "yes" ]]; then
    log_info "Cleaning up intermediate files..."
    rm -f "${BASENAME}_cleaned.fastq.gz" "${BASENAME}_raw.bed" "${BASENAME}_parsed.bed"
fi

log_info "Analysis Finished Successfully!"
log_info "Results in: $(pwd)"

# Calculate Duration
PIPELINE_END=$(date +%s)
DURATION=$((PIPELINE_END - PIPELINE_START))
H=$((DURATION/3600))
M=$(( (DURATION%3600)/60 ))
S=$((DURATION%60))

log_info "End Time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Total Duration: ${H}h ${M}m ${S}s"

send_notification "CLIPittyClip" "Analysis finished for $BASENAME. Duration: ${H}h ${M}m ${S}s"
