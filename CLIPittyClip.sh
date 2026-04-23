#!/bin/bash

# CLIPittyClip.sh - Modernized CLIP-seq Analysis Pipeline
# Version 3.1.0

# ------------------------------------------------------------------
# Initialization & Setup
# ------------------------------------------------------------------

# Resolve script directory to source libraries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/utils.sh"
source "${SCRIPT_DIR}/lib/modules.sh"
source "${SCRIPT_DIR}/lib/wizard.sh"

# Console output capture will be enabled after argument parsing (only for parent process)

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
BC_LEN=""
SPACER_LEN="0"
ADAPTER_3="GTGTCAGTCACTTCCAGCGG" # L32 default
PEAK_DIST=50
PEAK_SIZE=20
FRAG_LEN=25
PEAK_CALLER="homer" # Peak caller: homer (default) or ctk
ADV_PEAK_CALLER_ARGS=""    # Additional arguments for peak caller (homer or ctk)
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
    echo "  -i, --input-file <path>  Input FASTQ file (gzipped supported)"
    echo "  -d, --input-dir <dir>    Input directory containing .fastq.gz files (Batch Mode)"
    echo ""
    echo "REQUIRED REFERENCE:"
    echo "  -x, --index <path>       Path to genome index directory (STAR or Bowtie2)"
    echo ""
    echo "GENERAL OPTIONS:"
    echo "  -o, --output <str>       Output folder name (default: INPUT_output)"
    echo "  -t, --threads <int>      Number of threads (default: 1)"
    echo "  -m, --mapper <str>       Mapper: 'star' (default) or 'bowtie2'"
    echo "  -v, --verbose            Enable verbose logging"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "PREPROCESSING OPTIONS:"
    echo "  -u, --umi-length <int>   UMI length (e.g., 7 for CoCLIP)"
    echo "  --bc-len <int>           Barcode length to trim (auto-detected if -b is provided)"
    echo "  --spacer-len <int>       Spacer length to trim after barcode (default: 0)"
    echo "  -a, --adapter <str>      3' adapter sequence (default: L32)"
    echo "  --no-dedup               Disable FASTQ deduplication (default: ON)"
    echo "  --eclip                  ENCODE eCLIP mode: UMI in header, uses all eCLIP adapters"
    echo ""
    echo "DEMULTIPLEXING OPTIONS (for -i mode):"
    echo "  -b, --barcodes <path>    Barcode file for demultiplexing"
    echo "  --demux-mismatches <int> Max barcode mismatches (default: 1)"
    echo "  --align-mismatches <int> Max alignment mismatches (default: 2, STAR only)"
    echo ""
    echo "ANALYSIS OPTIONS:"
    echo "  --peak-caller <str>      Peak caller: homer (default) or ctk"
    echo "  --peak-caller-args <str> Additional peak caller arguments (quoted string)"
    echo "  --run-cims-cits          Enable full CTK CIMS+CITS analysis"
    echo "  --run-cims               Enable CIMS analysis (mutation sites only)"
    echo "  --run-cits               Enable CITS analysis (truncation sites only)"
    echo "  --cims-iter <int>        CIMS permutation iterations (default: 5)"
    echo "  --cims-fdr <float>       CIMS FDR threshold (default: 0.05)"
    echo "  --cits-pval <float>      CITS p-value threshold (default: 0.05)"
    echo "  --cits-gap <int>         CITS clustering gap (default: 25, -1=no cluster)"
    echo "  -f, --flank <int>        Flanked BED nucleotides (default: 10)"
    echo "  --no-motif               Skip flanked BED generation"
    echo "  -g, --groups <file>      Groups file for bedgraph/peak grouping"
    echo "  --ctk-group              Enable group CTK analysis (pools samples in groups.txt)"
    echo "  -s, --sample <int>       Test mode: process only first N reads"
    echo "  --filter-ncrna           Enable ncRNA pre-filtering (off by default)"
    echo ""
    echo "OUTPUT OPTIONS:"
    echo "  -k, --keep               Keep intermediate files (in OUTPUT/OTHERS/sample_analysis/)"
    echo "  --notification           Enable system notifications on completion"
    echo "  -w, --wizard             Launch interactive configuration wizard"
    echo ""
    echo "EXAMPLES:"
    echo "  # Standard run (Single File)"
    echo "  $0 -i reads.fastq.gz -x /path/to/star_index -t 8 -u 7"
    echo ""
    echo "  # Demultiplexing run"
    echo "  $0 -i pool.fastq.gz -b barcodes.txt -x /path/to/star_index -t 8"
    echo ""
    echo "  # Directory Batch Mode"
    echo "  $0 -d /path/to/samples/ -x /path/to/star_index -t 8 --run-cims"
    echo ""
    echo "  # Bowtie2 Alignment"
    echo "  $0 -i reads.fastq.gz -x /path/to/bt2_index -t 8 -m bowtie2"
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
FILTER_NCRNA="false"  # ncRNA filtering is OFF by default; enable with --filter-ncrna
ALIGNER="star" # Default aligner
ECLIP_MODE="false"  # ENCODE eCLIP mode (UMI in header, multiple adapters)
FILTER_CHR="true"   # Filter to canonical chromosomes (chr1-22, X, Y, M) - default ON
DEMUX_MISMATCHES="1"   # Default for barcode demultiplexing
ALIGN_MISMATCHES="2"   # Default for STAR --outFilterMismatchNmax

# CTK CIMS/CITS Parameters (with defaults)
CIMS_ITERATIONS="5"
CIMS_FDR="0.05"
CITS_PVALUE="0.05"
CITS_GAP="25"
RUN_MOTIF="yes"
MOTIF_FLANK="10"
CTK_GROUPS_FILE=""  # Optional: group samples for CIMS/CITS aggregation (set by --ctk-group)
CTK_GROUP_MODE="false"  # Explicit flag for group CTK analysis
GROUPS_FILE=""      # Standard Groups File for BedGraph/Matrix aggregation

# Capture Start Time (Seconds) for duration calculation
PIPELINE_START=$(date +%s)

if [[ $# -eq 0 ]]; then
    # Only show header if we are going to exit anyway
    echo "$separator_line"
    echo -e "${BLUE}CLIPittyClip: Modern CLIP-seq Analysis Pipeline${NC}"
    echo "Version 3.1.0"
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
        -i|--input-file) INPUT_FILE="$2"; shift 2 ;;
        -x|--index) GENOME_INDEX="$2"; shift 2 ;;
        -m|--mapper) ALIGNER=$(echo "$2" | tr '[:upper:]' '[:lower:]'); shift 2 ;;
        -o|--output) EXP_ID="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -u|--umi-length) UMI_LEN="$2"; shift 2 ;;
        --bc-len) BC_LEN="$2"; shift 2 ;;
        --spacer-len) SPACER_LEN="$2"; shift 2 ;;
        -a|--adapter) ADAPTER_3="$2"; shift 2 ;;
        -k|--keep) KEEP_INTERMEDIATE="yes"; shift ;;
        --peak-caller) PEAK_CALLER=$(echo "$2" | tr '[:upper:]' '[:lower:]'); shift 2 ;;
        --peak-caller-args) ADV_PEAK_CALLER_ARGS="$2"; shift 2 ;;
        --run-cims) RUN_CIMS=true; RUN_CTK="yes"; shift ;;
        --run-cits) RUN_CITS=true; RUN_CTK="yes"; shift ;;
        --run-cims-cits) RUN_CTK="yes"; RUN_CIMS=true; RUN_CITS=true; shift ;;
        --run-ctk) # Deprecated alias for --run-cims-cits
            echo -e "\033[0;33m[DEPRECATED]\033[0m --run-ctk is deprecated. Use --run-cims-cits instead." >&2
            RUN_CTK="yes"; RUN_CIMS=true; RUN_CITS=true; shift ;;
        --cims-iter) CIMS_ITERATIONS="$2"; shift 2 ;;
        --cims-fdr) CIMS_FDR="$2"; shift 2 ;;
        --cits-pval) CITS_PVALUE="$2"; shift 2 ;;
        --cits-gap) CITS_GAP="$2"; shift 2 ;;
        -f|--flank) MOTIF_FLANK="$2"; shift 2 ;;
        --no-motif) RUN_MOTIF="no"; shift ;;
        -g|--groups) GROUPS_FILE="$2"; shift 2 ;;
        --ctk-group) CTK_GROUP_MODE="true"; shift ;;
        -s|--sample) SAMPLE_SIZE="$2"; shift 2 ;;
        -b|--barcodes) BARCODE_FILE="$2"; DEMUX="yes"; shift 2 ;;
        -d|--input-dir) INPUT_DIR="$2"; shift 2 ;;
        --demux-mismatches) DEMUX_MISMATCHES="$2"; shift 2 ;;
        --align-mismatches) ALIGN_MISMATCHES="$2"; shift 2 ;;
        --no-dedup) DEDUP_MODE="false"; shift ;;
        --filter-ncrna) FILTER_NCRNA="true"; shift ;;
        --eclip) ECLIP_MODE="true"; shift ;;
        --no-chr-filter) FILTER_CHR="false"; shift ;;
        --notification) NOTIFY_MODE="true"; shift ;;
        --child) CHILD_MODE="true"; shift ;;
        -w|--wizard|--advanced) WIZARD_MODE="true"; shift ;;
        -v|--verbose) VERBOSE="true"; shift ;;
        -h|--help) show_usage; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# ------------------------------------------------------------------
# Pre-Wizard Config Checking
# ------------------------------------------------------------------
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
    ADV_PEAK_CALLER_ARGS="$WIZ_HOMER_ARGS"
    
    if [[ -n "$WIZ_PEAK_CALLER" ]]; then
        PEAK_CALLER="$WIZ_PEAK_CALLER"
    fi
    
    if [[ -n "$WIZ_CTK_PEAK_ARGS" ]]; then
        # Prefix the arguments so they append properly to default CTK args
        ADV_PEAK_CALLER_ARGS="$WIZ_CTK_PEAK_ARGS"
    fi
    
    # Important bugfix: fastp arguments from advanced wizard
    # We map WIZ_FASTP_ARGS so the standard CLI loop can process it
    if [[ -n "$WIZ_FASTP_ARGS" ]]; then
        ADV_FASTP_ARGS="$WIZ_FASTP_ARGS"
    fi
fi

# ------------------------------------------------------------------
# Log File Initialization (Post-Parsing / Post-Wizard)
# ------------------------------------------------------------------
# If CHILD_MODE, we suppress main log creation to avoid spamming root dir.
if [[ "$CHILD_MODE" == "true" ]]; then
    LOG_FILE="/dev/null"
    TEMP_CONSOLE_LOG=""  # No console log for child processes
else
    # Enable console output capture for parent process only
    # Derived from user requested fixed name
    TEMP_CONSOLE_LOG="console_output.log"
    
    # Enable tee redirection AFTER wizard to prevent capturing UI ANSI color codes
    > "$TEMP_CONSOLE_LOG"
    exec > >(tee -a "$TEMP_CONSOLE_LOG") 2>&1

    # Log file setup
    # User Request: Fixed name "detailed_output.log"
    LOG_FILE="$(pwd)/detailed_output.log"
    # Overwrite if exists
    > "${LOG_FILE}"
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

# Barcode length validation
if [[ -n "$BARCODE_FILE" ]]; then
    bc_file_len=$(awk '!/^#/{print length($2); exit}' "$BARCODE_FILE")
    if [[ -z "$bc_file_len" ]]; then
        log_error "Failed to read barcode length from $BARCODE_FILE"
        exit 1
    fi
    if [[ -n "$BC_LEN" && "$BC_LEN" != "$bc_file_len" ]]; then
        log_error "Conflict: --bc-len ($BC_LEN) does not match the barcode length in file ($bc_file_len)."
        exit 1
    fi
    BC_LEN="$bc_file_len"
fi
BC_LEN="${BC_LEN:-0}"

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

# Normalize EXP_ID (strip trailing slash if present)
EXP_ID="${EXP_ID%/}"

# Determine input parent directory for default output location
if [[ -n "$INPUT_FILE" ]]; then
    INPUT_PARENT="$(dirname "$INPUT_FILE")"
elif [[ -n "$INPUT_DIR" ]]; then
    INPUT_PARENT="$(dirname "$INPUT_DIR")"
fi

# Thread validation: cap to available cores - 1 (leave 1 for system)
# Cross-platform CPU detection: macOS uses sysctl, Linux uses nproc
if [[ "$(uname)" == "Darwin" ]]; then
    MAX_AVAILABLE_THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    MAX_AVAILABLE_THREADS=$(nproc 2>/dev/null || echo 4)
fi
MAX_SAFE_THREADS=$(( MAX_AVAILABLE_THREADS - 1 ))
[[ $MAX_SAFE_THREADS -lt 1 ]] && MAX_SAFE_THREADS=1

if [[ $THREADS -gt $MAX_SAFE_THREADS ]]; then
    log_warning "Requested $THREADS threads but only $MAX_AVAILABLE_THREADS cores available."
    log_warning "Capping to $MAX_SAFE_THREADS threads (leaving 1 for system stability)."
    THREADS=$MAX_SAFE_THREADS
fi

# Validate groups file (requires --run-cims or --run-cits)
if [[ "$CTK_GROUP_MODE" == "true" ]] && [[ -z "$GROUPS_FILE" ]]; then
    log_error "--ctk-group requires a groups file (-g). Please provide a groups file."
    show_usage
    exit 1
fi

if [[ -n "$GROUPS_FILE" ]]; then
    if [[ "$RUN_CIMS" != "true" ]] && [[ "$RUN_CITS" != "true" ]] && [[ "$CTK_GROUP_MODE" == "true" ]]; then
        log_warning "--ctk-group requires --run-cims and/or --run-cits. Group CTK will be disabled."
        CTK_GROUP_MODE="false"
    fi
    if [[ ! -f "$GROUPS_FILE" ]]; then
        log_error "Groups file not found: $GROUPS_FILE"
        exit 1
    fi
    GROUPS_FILE="$(cd "$(dirname "$GROUPS_FILE")" && pwd)/$(basename "$GROUPS_FILE")"
    # Only set CTK_GROUPS_FILE if group CTK mode is explicitly enabled
    if [[ "$CTK_GROUP_MODE" == "true" ]]; then
        CTK_GROUPS_FILE="$GROUPS_FILE"
        log_info "Group CTK analysis enabled with groups file: $CTK_GROUPS_FILE"
    fi
    log_info "Groups file: $GROUPS_FILE"
fi

# ------------------------------------------------------------------
# Log Configuration (Standard or Advanced Results)
# ------------------------------------------------------------------
log_info "------------------------------------------------------------------"
log_info "Configuration Summary"
log_info "------------------------------------------------------------------"

# Define defaults for display (matching lib/modules.sh architecture)
DEF_FASTP="--length_required 16 --average_qual 30"
DEF_STAR="--outFilterMultimapNmax 10 --outFilterMismatchNmax ${ALIGN_MISMATCHES} --alignEndsType EndToEnd --outSAMattributes ... MD"
DEF_BT2="--md --end-to-end (Standard Sensitivity)"
DEF_HOMER="-style factor -L 2 -localSize 10000 -minDist ${PEAK_DIST:-50}"
DEF_CTK="-big -ss --valley-seeking -minPH 2 -gap ${PEAK_DIST:-50}"

if [[ "$ADVANCED_MODE" == "true" ]]; then
    log_info "Mode:           ADVANCED (Using Overrides)"
    log_info "Fastp Defaults: $DEF_FASTP"
    log_info "Fastp Added:    ${ADV_FASTP_ARGS:-(None)}"

    log_info "Aligner:        $ALIGNER"
    if [[ "$ALIGNER" == "star" ]]; then log_info "Aligner Def:    $DEF_STAR"; else log_info "Aligner Def:    $DEF_BT2"; fi
    log_info "Aligner Added:  ${ADV_ALIGNER_ARGS:-(None)}"

    log_info "Peak Caller:    $PEAK_CALLER"
    if [[ "$PEAK_CALLER" == "ctk" ]]; then
        log_info "CTK Defaults:   $DEF_CTK"
        log_info "CTK Added:      ${ADV_PEAK_CALLER_ARGS:-(None)}"
    else
        log_info "Homer Defaults: $DEF_HOMER"
        log_info "Homer Added:    ${ADV_PEAK_CALLER_ARGS:-(None)}"
    fi
else
    log_info "Mode:           STANDARD"
    log_info "Fastp Config:   $DEF_FASTP"
    log_info "Aligner:        $ALIGNER"
    if [[ "$ALIGNER" == "star" ]]; then log_info "Aligner Config: $DEF_STAR"; else log_info "Aligner Config: $DEF_BT2"; fi
    log_info "Peak Caller:    $PEAK_CALLER"
    if [[ "$PEAK_CALLER" == "ctk" ]]; then
        log_info "CTK Config:     $DEF_CTK"
        log_info "CTK Added:      ${ADV_PEAK_CALLER_ARGS:-(None)}"
    else
        log_info "Homer Config:   $DEF_HOMER"
    fi
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
if [[ "$PEAK_CALLER" == "ctk" ]]; then
    check_dependency tag2peak.pl
fi
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
    
    # eCLIP Mode Banner
    if [[ "$ECLIP_MODE" == "true" ]]; then
        console_msg "╔══════════════════════════════════════════════════════════════════╗"
        console_msg "║                     eCLIP MODE ENABLED                           ║"
        console_msg "╠══════════════════════════════════════════════════════════════════╣"
        console_msg "║  • UMI extraction: DISABLED (UMI already in read header)         ║"
        console_msg "║  • Adapter trimming: Using 9 standard eCLIP adapters             ║"
        console_msg "║  • Header reformat: Converting to CTK-compatible format          ║"
        console_msg "╚══════════════════════════════════════════════════════════════════╝"
        if [[ "$UMI_LEN" -gt 0 ]]; then
            console_msg "[WARNING] eCLIP mode: -u ${UMI_LEN} will be ignored (UMI in header)"
        fi
        if [[ "$ADAPTER_3" != "GTGTCAGTCACTTCCAGCGG" ]]; then
            console_msg "[WARNING] eCLIP mode: -a will be ignored (using eCLIP adapters)"
        fi
        console_msg "Threads: $THREADS | Mode: eCLIP (ENCODE pre-processed)"
    else
        # Display adapter: show "L32" if default, or full custom sequence
        adapter_display="L32"
        if [[ "$ADAPTER_3" != "GTGTCAGTCACTTCCAGCGG" ]]; then
            adapter_display="$ADAPTER_3"
        fi
        console_msg "Threads: $THREADS | UMI: ${UMI_LEN}bp | Adapter: $adapter_display"
    fi
    # Note about dynamic thread scaling for CIMS/CITS
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        console_msg "Note: CIMS/CITS threads may be reduced based on available memory"
    fi
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
    
    # Collect .fastq.gz, .fq.gz, .fastq, and .fq files
    SAMPLE_FILES=($(ls "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz \
                      "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq 2>/dev/null))
    total_samples=${#SAMPLE_FILES[@]}

    console_msg "  > Found $total_samples sample files"
    
    # Build extra flags for child processes
    EXTRA_FLAGS=""
    # Only pass CIMS/CITS to children if NOT using groups file (group CTK runs after batch)
    if [[ -z "$CTK_GROUPS_FILE" ]]; then
        if [[ "$RUN_CIMS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --run-cims"; fi
        if [[ "$RUN_CITS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --run-cits"; fi
    fi
    if [[ "$VERBOSE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --verbose"; fi
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS -k"; fi
    if [[ "$SAMPLE_SIZE" -gt 0 ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --sample $SAMPLE_SIZE"; fi
    EXTRA_FLAGS="$EXTRA_FLAGS -m $ALIGNER"
    if [[ "$ECLIP_MODE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --eclip"; fi
    if [[ "$FILTER_NCRNA" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --filter-ncrna"; fi
    if [[ "$DEDUP_MODE" == "false" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --no-dedup"; fi
    EXTRA_FLAGS="$EXTRA_FLAGS --peak-caller $PEAK_CALLER"
    if [[ -n "$ADV_PEAK_CALLER_ARGS" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --peak-caller-args \"$ADV_PEAK_CALLER_ARGS\""; fi
    if [[ -n "$BC_LEN" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --bc-len $BC_LEN"; fi
    if [[ -n "$SPACER_LEN" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --spacer-len $SPACER_LEN"; fi
    EXTRA_FLAGS="$EXTRA_FLAGS --child"
    
    console_msg "\n[BATCH ANALYSIS]"

    current_sample=0
    for sample in "${SAMPLE_FILES[@]}"; do
        if [ -f "$sample" ]; then
            ((current_sample++))
            sample_name=$(basename "$sample")
            sample_name="${sample_name%.fastq.gz}"
            sample_name="${sample_name%.fq.gz}"
            sample_name="${sample_name%.fastq}"
            sample_name="${sample_name%.fq}"

            # Auto-gzip plain .fastq/.fq files (preserve originals)
            sample_input="$sample"
            gzip_tmp=""
            if [[ "$sample" != *.gz ]]; then
                console_msg "  > Compressing $(basename "$sample") for processing..."
                gzip_tmp="${sample_name}_tmp_input.fastq.gz"
                gzip -c "$sample" > "$gzip_tmp"
                sample_input="$gzip_tmp"
                log_info "Auto-gzipped $sample -> $gzip_tmp (original preserved)"
            fi

            printf "  %2d/%-2d %-20s : " "$current_sample" "$total_samples" "$sample_name"

            echo "[BATCH] Launching analysis for sample: $sample_name" >> "$LOG_FILE"

            self_script="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"

            cmd="bash $self_script \
                -i $sample_input \
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

            # Cleanup gzip temp if created
            if [[ -n "$gzip_tmp" && -f "$gzip_tmp" ]]; then
                rm -f "$gzip_tmp"
                log_info "Removed temp gzip: $gzip_tmp"
            fi
        fi
    done
    
    # Aggregation - same as demux path
    INPUT_BASENAME=$(basename "$INPUT_DIR")
    # Use -o value for output folder if provided, otherwise use input directory name
    # If -o contains a path separator, use as full path; otherwise create next to input
    if [[ -n "$EXP_ID" ]]; then
        if [[ "$EXP_ID" == */* ]]; then
            OUTPUT_ROOT="$EXP_ID"
        else
            OUTPUT_ROOT="${INPUT_PARENT}/${EXP_ID}"
        fi
    else
        OUTPUT_ROOT="${INPUT_PARENT}/${INPUT_BASENAME}_output"
    fi
    
    DIR_BAM="1_BAM"
    DIR_BED="2_COLLAPSED_BED"
    DIR_BG="3_BEDGRAPH"
    DIR_PEAKS="4_PEAKS"
    # OTHERS folder number depends on whether CTK analysis is enabled
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        DIR_OTHERS="6_OTHERS"
        # Determine CTK folder name for PEAKittyPeak.sh
        if [[ "$RUN_CIMS" == "true" ]] && [[ "$RUN_CITS" == "true" ]]; then
            DIR_CTK="5_CTK_Analysis"
        elif [[ "$RUN_CIMS" == "true" ]]; then
            DIR_CTK="5_CIMS_Analysis"
        else
            DIR_CTK="5_CITS_Analysis"
        fi
    else
        DIR_OTHERS="5_OTHERS"
        DIR_CTK="" # No CTK analysis, so no CTK folder
    fi
    DIR_REPORTS="REPORTS"
    DIR_PEAK_LOGS="REPORTS/PEAK"
    DIR_IND_PEAK_LOGS="$DIR_PEAK_LOGS"
    
    # Create base directories (CTK folders created conditionally during aggregation)
    mkdir -p "$OUTPUT_ROOT/$DIR_BAM" \
             "$OUTPUT_ROOT/$DIR_BED" \
             "$OUTPUT_ROOT/$DIR_BG" \
             "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS" \
             "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS" \
             "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/DETAILED_LOGS_CAN_BE_DELETED" \
             "$OUTPUT_ROOT/$DIR_REPORTS/SAMPLES" \
             "$OUTPUT_ROOT/$DIR_PEAK_LOGS" \
             "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS"
    
    # Create aligner-specific folders
    if [[ "$ALIGNER" != "bowtie2" ]]; then
        mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT"
    fi

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
            
            # Use find-based collection (matches demux mode's working approach)
            # Child processes create flat structure in _analysis/, not nested folders
            
            # 1. BAM files
            bam_file=$(find "$sample_out" -name "*mapped.Aligned.sortedByCoord.out.bam" 2>/dev/null | head -n 1)
            if [[ -n "$bam_file" ]]; then
                cp "$bam_file" "$OUTPUT_ROOT/$DIR_BAM/${sample_name}.bam"
                cp "${bam_file}.bai" "$OUTPUT_ROOT/$DIR_BAM/${sample_name}.bam.bai" 2>/dev/null
            fi
            
            # 2. Collapsed BED files
            bed_file=$(find "$sample_out" -name "*_collapsed.bed" 2>/dev/null | head -n 1)
            if [[ -n "$bed_file" ]]; then
                cp "$bed_file" "$OUTPUT_ROOT/$DIR_BED/${sample_name}.bed"
            fi
            
            # 3. Bedgraph files & chrom.sizes
            bg_pos=$(find "$sample_out" -name "*_pos.bedgraph" 2>/dev/null | head -n 1)
            bg_neg=$(find "$sample_out" -name "*_neg.bedgraph" 2>/dev/null | head -n 1)
            if [[ -n "$bg_pos" ]]; then cp "$bg_pos" "$OUTPUT_ROOT/$DIR_BG/${sample_name}_pos.bedgraph"; fi
            if [[ -n "$bg_neg" ]]; then cp "$bg_neg" "$OUTPUT_ROOT/$DIR_BG/${sample_name}_neg.bedgraph"; fi
            
            # Collect scale_factors.tsv (append to master file)
            scale_file=$(find "$sample_out" -name "scale_factors.tsv" 2>/dev/null | head -n 1)
            if [[ -f "$scale_file" ]]; then
                cat "$scale_file" >> "$OUTPUT_ROOT/$DIR_BG/scale_factors.tsv"
            fi
            
            # Grab chrom.sizes
            extracted_chroms=$(find "$sample_out" -name "chrom.sizes" 2>/dev/null | head -n 1)
            if [[ -f "$extracted_chroms" ]]; then
                cp "$extracted_chroms" "$OUTPUT_ROOT/$DIR_BG/chrom.sizes"
            fi
            
            # 4. Peak folder (HOMER produces a directory; CTK produces flat BED + log files)
            peak_dir="${sample_out}/${sample_name}_peaks"
            if [[ -d "$peak_dir" ]]; then
                # HOMER: copy the tag directory
                cp -r "$peak_dir" "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/"
                peak_log="${peak_dir}_homer.log"
            else
                # CTK: copy the flat BED output into a named subdirectory
                local ctk_peak_bed="${sample_out}/${sample_name}_peaks_raw.bed"
                if [[ -f "$ctk_peak_bed" ]]; then
                    mkdir -p "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/${sample_name}_peaks"
                    cp "$ctk_peak_bed" "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/${sample_name}_peaks/"
                fi
                peak_log="${sample_out}/${sample_name}_peaks_ctk.log"
            fi

            # Peak calling log
            if [[ -f "$peak_log" ]]; then
                cp "$peak_log" "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS/${sample_name}_PeakCalling.log"
            fi
            
            # 5. CTK Analysis (CIMS/CITS) - create per-sample subdirectory to avoid overwrites
            if [[ -n "$DIR_CTK" ]]; then
                for ctk_folder in "CTK_Analysis" "CIMS_Analysis" "CITS_Analysis"; do
                    if [[ -d "$sample_out/$ctk_folder" ]]; then
                        # Create sample-specific subdirectory
                        mkdir -p "$OUTPUT_ROOT/$DIR_CTK/${sample_name}"
                        cp -r "$sample_out/$ctk_folder/"* "$OUTPUT_ROOT/$DIR_CTK/${sample_name}/" 2>/dev/null
                        break
                    fi
                done
            fi
            
            # 6. Reports & Logs
            cp "$sample_out"/*_fastp.html "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
            cp "$sample_out"/*_fastp.json "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
            
            # Aligner-specific logs
            if [[ "$ALIGNER" == "bowtie2" ]]; then
                echo -e "Bowtie2 does not generate standalone log files like STAR does.\nAll mapping statistics for this run are recorded in the main detailed_output.log." \
                    > "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/README_Bowtie2_LOGS.txt"
            else
                cp "$sample_out"/*Log.final.out "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/" 2>/dev/null
                cp "$sample_out"/*Log.out "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/DETAILED_LOGS_CAN_BE_DELETED/" 2>/dev/null
            fi
            
            # ncRNA Mapping outputs
            if [[ -d "$sample_out/OTHERS/ncRNA_Mapping" ]]; then
                mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping"
                cp "$sample_out"/OTHERS/ncRNA_Mapping/*_ncrna_stats.txt "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                cp "$sample_out"/OTHERS/ncRNA_Mapping/*_ncrna.bam* "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
            fi
            
            # Pipeline log
            cp "$sample_out"/*.log "$OUTPUT_ROOT/$DIR_REPORTS/SAMPLES/${sample_name}_detailed.log" 2>/dev/null
            
            # Cleanup sample dir or move to OTHERS if -k
            # Cleanup sample dir or move to OTHERS if -k
            if [[ "$KEEP_INTERMEDIATE" != "yes" ]]; then
                rm -rf "$sample_out"
            else
                # Move intermediate files into OUTPUT/OTHERS/ for cleaner organization
                mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/sample_analysis"
                mv "$sample_out" "$OUTPUT_ROOT/$DIR_OTHERS/sample_analysis/" 2>/dev/null
            fi
        fi
    done
    
    # Combined BedGraph Generation (Directory Mode)
    # Only runs when --groups/-g is explicitly provided
    if [[ -n "$GROUPS_FILE" ]]; then
        console_msg "  > Generating Combined Bedgraph by groups..."
        run_combined_bedgraph "$OUTPUT_ROOT" "$GROUPS_FILE" "$OUTPUT_ROOT/$DIR_BG"
    fi
    
    # ncRNA Filtering Summary
    console_msg "\n[ncRNA FILTERING SUMMARY]"
    printf "  %-25s %-15s %-15s %s\n" "Sample" "ncRNA Reads" "Total Reads" "% Filtered"
    console_msg "  ----------------------------------------------------------------"
    
    for sample in "${SAMPLE_FILES[@]}"; do
        if [ -f "$sample" ]; then
            sample_name=$(basename "$sample")
            sample_name="${sample_name%.fastq.gz}"
            sample_name="${sample_name%.fq.gz}"
            # Read from aggregated location (sample_out was deleted after collection)
            ncrna_stats="${OUTPUT_ROOT}/${DIR_OTHERS}/ncRNA_Mapping/${sample_name}_ncrna_stats.txt"
            
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
    
    # Combined Peak Calling
    console_msg "\n[COMBINED PEAK CALLING]"
    BED_DIR="$OUTPUT_ROOT/$DIR_BED"
    
    PEAK_LOG="$OUTPUT_ROOT/$DIR_PEAK_LOGS/Combined_PeakCalling.log"
    
    PEAK_CMD="$SCRIPT_DIR/PEAKittyPeak.sh -i \"$BED_DIR\" --aggregate -n \"Combined\" -p \"$PEAK_DIST\" -z \"$PEAK_SIZE\" -f \"$FRAG_LEN\" --peak-caller \"$PEAK_CALLER\""
    if [[ -n "$ADV_PEAK_CALLER_ARGS" ]]; then PEAK_CMD="$PEAK_CMD --peak-caller-args \"$ADV_PEAK_CALLER_ARGS\""; fi
    if [[ -n "$DIR_CTK" ]]; then
        PEAK_CMD="$PEAK_CMD --ctk-dir \"$OUTPUT_ROOT/$DIR_CTK\""
    fi
    if [[ -n "$CTK_GROUPS_FILE" ]]; then
        PEAK_CMD="$PEAK_CMD --ctk-group \"$CTK_GROUPS_FILE\""
    fi
    if [[ "$RUN_CIMS" == "true" ]]; then
        PEAK_CMD="$PEAK_CMD --cims-fdr \"$CIMS_FDR\""
    fi
    if [[ "$RUN_CITS" == "true" ]]; then
        PEAK_CMD="$PEAK_CMD --cits-pval \"$CITS_PVALUE\""
    fi
    
    # Force hardcoded path to ensure correct filename
    eval "$PEAK_CMD" > "$OUTPUT_ROOT/$DIR_PEAK_LOGS/Combined_PeakCalling.log" 2>&1
    
    # Move combined results (PEAKittyPeak creates ${BASE_NAME}_peaks/, not PEAKS/)
    if [[ -d "Combined_peaks" ]]; then
        mv Combined_peaks/* "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/" 2>/dev/null
        rm -rf Combined_peaks
    fi
    
    console_msg "  > Combined peaks: $OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/"
    
    # Final Summary (print BEFORE moving logs so messages get captured)
    PIPELINE_END=$(date +%s)
    DURATION=$((PIPELINE_END - PIPELINE_START))
    H=$((DURATION/3600))
    M=$(( (DURATION%3600)/60 ))
    S=$((DURATION%60))
    
    console_msg "\n[COMPLETE]"
    console_msg "  > Duration: ${H}h ${M}m ${S}s"
    console_msg "  > Output: $OUTPUT_ROOT/"
    console_msg "  > Console log: $OUTPUT_ROOT/$DIR_REPORTS/${TEMP_CONSOLE_LOG:-CLIPittyClip.log}"
    

    
    # Remove sampled fastq files (created when --sample is used)
    rm -f *_sampled_*.fastq.gz 2>/dev/null
    
    # Remove CITS.pl temporary directories and .tmp files (from failed/empty CTK runs)
    rm -rf CITS.pl_* tag2profile.pl_* 2>/dev/null
    rm -f *.tmp 2>/dev/null
    
    # Add enhanced columns to peak matrix (TC, NC, BC, Cov)
    # This fixes the missing columns issue for eCLIP/Batch runs
    PEAK_MATRIX="$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/COMBINED_peakMatrix.txt"
    PEAKS_BED="$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/peaks_Sorted.bed"
    
    if [[ -f "$PEAK_MATRIX" && -f "$PEAKS_BED" ]]; then
        console_msg "  > Adding enhanced matrix columns..."
        add_matrix_columns "$PEAK_MATRIX" "$PEAKS_BED" \
            "$OUTPUT_ROOT/$DIR_BG" "$OUTPUT_ROOT/$DIR_BG/scale_factors.tsv" \
            "$GROUPS_FILE"
        
        # Promote sorted peak bed to clearly named FINAL target
        mv "$PEAKS_BED" "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed" 2>/dev/null
        
        # Cleanup intermediate peak files to save disk space
        rm -f "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/COMBINED.bed" \
              "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/peaks.bed"
    fi

    send_notification "CLIPittyClip Batch" "Analysis complete for $total_samples samples in $(basename "$INPUT_DIR")"
    
    # Cleanup: Move logs and remove intermediate files (AFTER final messages)
    # Move detailed log (if defined and exists)
    if [[ -f "$LOG_FILE" ]]; then
        # Ensure it goes to REPORTS/detailed_output.log
        mv "$LOG_FILE" "$OUTPUT_ROOT/$DIR_REPORTS/detailed_output.log" 2>/dev/null
        # Update variable so any final messages go to the right place
        LOG_FILE="$OUTPUT_ROOT/$DIR_REPORTS/detailed_output.log"
    else
        mv CLIPittyClip_*.log "$OUTPUT_ROOT/$DIR_REPORTS/" 2>/dev/null
    fi
    
    # Move console log into output
    if [[ -n "$TEMP_CONSOLE_LOG" && -f "$TEMP_CONSOLE_LOG" ]]; then
        mv "$TEMP_CONSOLE_LOG" "$OUTPUT_ROOT/$DIR_REPORTS/console_output.log" 2>/dev/null
    fi
    
    exit 0
fi

# ── Auto-gzip plain .fastq/.fq files for -i mode ────────────────────────────────
GLOBAL_GZIP_TMP=""
if [[ -n "$INPUT_FILE" ]] && [[ "$INPUT_FILE" != *.gz ]]; then
    tmp_base=$(basename "$INPUT_FILE")
    tmp_base="${tmp_base%.fastq}"
    tmp_base="${tmp_base%.fq}"
    
    console_msg "\n[PRE-PROCESSING]"
    console_msg "  > Compressing $(basename "$INPUT_FILE") for processing..."
    GLOBAL_GZIP_TMP="$(pwd)/${tmp_base}_tmp_input.fastq.gz"
    gzip -c "$INPUT_FILE" > "$GLOBAL_GZIP_TMP"
    INPUT_FILE="$GLOBAL_GZIP_TMP"
    log_info "Auto-gzipped original input -> $GLOBAL_GZIP_TMP (original preserved)"
fi

# 0b. Demultiplexing (Recursive Branch)
if [[ "$DEMUX" == "yes" ]]; then
    
    # 1a. Deduplication (if enabled)
    if [[ "$DEDUP_MODE" == "true" ]]; then
        console_msg "\n[DEDUPLICATING]"
        print_section_item "Deduplicating Pooled Reads (fastq2collapse.pl)"
        DEDUP_OUT="pooled_dedup.fastq.gz"
        if run_dedup "$INPUT_FILE" "$DEDUP_OUT"; then
            WORK_INPUT="$DEDUP_OUT"
            print_section_item "Deduplication Complete"
        else
            log_warning "Deduplication failed. Using original file."
            rm -f "$DEDUP_OUT"
            WORK_INPUT="$INPUT_FILE"
        fi
    else
        WORK_INPUT="$INPUT_FILE"
    fi
    
    # 1b. Demultiplexing
    console_msg "\n[DEMULTIPLEXING]"
    print_section_item "Barcodes: $(basename "$BARCODE_FILE")"
    print_section_item "Mismatches Allowed: $DEMUX_MISMATCHES"
    
    # Run demultiplexing (with dedup disabled since we already did it)
    run_demultiplexing "$WORK_INPUT" "$BARCODE_FILE" "$SAMPLE_SIZE" "false"
    
    print_section_item "Demultiplexing Complete"
    
    # Cleanup dedup temp file
    if [[ -f "pooled_dedup.fastq.gz" ]]; then
        rm -f "pooled_dedup.fastq.gz"
        log_info "Cleaned up pooled dedup temp file."
    fi
    
    send_notification "CLIPittyClip: $(basename "$INPUT_FILE")" "Demultiplexing complete"
    
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
    # Only pass CIMS/CITS to children if NOT using group CTK mode
    # When --ctk-group is enabled, skip individual CTK (group CTK runs after batch)
    if [[ "$CTK_GROUP_MODE" != "true" ]]; then
        if [[ "$RUN_CIMS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --run-cims"; fi
        if [[ "$RUN_CITS" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --run-cits"; fi
    else
        log_info "Group CTK mode: skipping individual CIMS/CITS (will run on grouped data)"
    fi
    if [[ "$VERBOSE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --verbose"; fi
    if [[ "$SAMPLE_SIZE" -gt 0 ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --sample $SAMPLE_SIZE"; fi
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS -k"; fi
    
    # Pass Aligner choice
    EXTRA_FLAGS="$EXTRA_FLAGS -m $ALIGNER"
    if [[ "$ECLIP_MODE" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --eclip"; fi
    if [[ "$FILTER_NCRNA" == "true" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --filter-ncrna"; fi
    # Pool was already deduped above; tell children to skip dedup
    EXTRA_FLAGS="$EXTRA_FLAGS --no-dedup"
    
    EXTRA_FLAGS="$EXTRA_FLAGS --peak-caller $PEAK_CALLER"
    if [[ -n "$ADV_PEAK_CALLER_ARGS" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --peak-caller-args \"$ADV_PEAK_CALLER_ARGS\""; fi
    if [[ -n "$BC_LEN" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --bc-len $BC_LEN"; fi
    if [[ -n "$SPACER_LEN" ]]; then EXTRA_FLAGS="$EXTRA_FLAGS --spacer-len $SPACER_LEN"; fi

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
            
            # Run the command with log capture
            # In CHILD_MODE, the script will produce `update_status` outputs (echo -ne ... \r).
            # We capture this into a sample-specific log while still showing it on console.
            $cmd 2>&1 | tee "${sample_name}.log"
            
            # Check exit status of the PIPELINE (not tee) using PIPESTATUS
            pipestatus="${PIPESTATUS[0]}"
            
            # When child finishes, we need to print "Done." and a newline to finalize the line.
            if [ $pipestatus -eq 0 ]; then
                update_status_done
                echo "[BATCH] Sample $sample_name analysis complete." >> "$LOG_FILE"
                
                # Move the captured log into the analysis directory if it exists
                if [[ -d "${sample_name}_analysis" ]]; then
                    mv "${sample_name}.log" "${sample_name}_analysis/${sample_name}.log"
                fi
            else
                echo -e "${RED}FAILED${NC}"
                log_error ">>> Sample $sample_name analysis FAILED."
                # Keep the log file in root for debugging
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
    # Use -o value for output folder if provided, otherwise use input file name
    # If -o contains a path separator, use as full path; otherwise create next to input
    if [[ -n "$EXP_ID" ]]; then
        if [[ "$EXP_ID" == */* ]]; then
            OUTPUT_ROOT="$EXP_ID"
        else
            OUTPUT_ROOT="${INPUT_PARENT}/${EXP_ID}"
        fi
    else
        OUTPUT_ROOT="${INPUT_PARENT}/${INPUT_BASENAME}_output"
    fi
    
    DIR_DEMUX="0_DEMUX_FASTQ"
    DIR_BAM="1_BAM"
    DIR_BED="2_COLLAPSED_BED"
    DIR_BG="3_BEDGRAPH"
    DIR_PEAKS="4_PEAKS"
    # OTHERS folder number depends on whether CTK analysis is enabled
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        DIR_OTHERS="6_OTHERS"
        # Determine CTK folder name for PEAKittyPeak.sh
        if [[ "$RUN_CIMS" == "true" ]] && [[ "$RUN_CITS" == "true" ]]; then
            DIR_CTK="5_CTK_Analysis"
        elif [[ "$RUN_CIMS" == "true" ]]; then
            DIR_CTK="5_CIMS_Analysis"
        else
            DIR_CTK="5_CITS_Analysis"
        fi
    else
        DIR_OTHERS="5_OTHERS"
        DIR_CTK="" # No CTK analysis, so no CTK folder
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
             "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT" \
             "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/DETAILED_LOGS_CAN_BE_DELETED" \
             "$OUTPUT_ROOT/$DIR_REPORTS/SAMPLES" \
             "$OUTPUT_ROOT/$DIR_PEAK_LOGS" \
             "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS"
    
    # Create aligner-specific folders
    if [[ "$ALIGNER" != "bowtie2" ]]; then
        mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT"
    fi

    # Run Group-based CTK Analysis (OPT-IN ONLY with --ctk-group flag)
    if [[ "$CTK_GROUP_MODE" == "true" ]] && [[ -n "$GROUPS_FILE" ]]; then
        run_group_ctk_analysis "$GROUPS_FILE" "$OUTPUT_ROOT" "$GENOME_INDEX" \
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
                
                # Collect scale_factors.tsv (append to master file)
                scale_file=$(find "$analysis_dir" -name "scale_factors.tsv" 2>/dev/null | head -n 1)
                if [[ -f "$scale_file" ]]; then
                    cat "$scale_file" >> "$OUTPUT_ROOT/$DIR_BG/scale_factors.tsv"
                fi
                
                # Grab chrom.sizes from ONE sample (overwrite is fine)
                 # It might be in analysis dir if we extracted it, or just use user provided if needed
                 # modules.sh extracts it to 'chrom.sizes' in CWD if not provided?
                 # Actually modules.sh:246 writes to "$chrom_sizes".
                 extracted_chroms="$analysis_dir/chrom.sizes"
                 if [[ -f "$extracted_chroms" ]]; then
                    cp "$extracted_chroms" "$OUTPUT_ROOT/$DIR_BG/chrom.sizes"
                 fi

                # 4. Peaks (HOMER produces a directory; CTK produces flat BED + log files)
                peak_dir="${analysis_dir}/${sample_name}_peaks"
                if [[ -d "$peak_dir" ]]; then
                    # HOMER: copy the tag directory
                    cp -r "$peak_dir" "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/"
                    peak_log="${peak_dir}_homer.log"
                else
                    # CTK: copy the flat BED output into a named subdirectory
                    local ctk_peak_bed="${analysis_dir}/${sample_name}_peaks_raw.bed"
                    if [[ -f "$ctk_peak_bed" ]]; then
                        mkdir -p "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/${sample_name}_peaks"
                        cp "$ctk_peak_bed" "$OUTPUT_ROOT/$DIR_PEAKS/SAMPLE_PEAKS/${sample_name}_peaks/"
                    fi
                    peak_log="${analysis_dir}/${sample_name}_peaks_ctk.log"
                fi

                # 4.1. Individual Peak Log
                if [[ -f "$peak_log" ]]; then
                    cp "$peak_log" "$OUTPUT_ROOT/$DIR_IND_PEAK_LOGS/${sample_name}_PeakCalling.log"
                fi

                # 5. Reports & Logs
                # FASTP
                cp "$analysis_dir"/*_fastp.html "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
                cp "$analysis_dir"/*_fastp.json "$OUTPUT_ROOT/$DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
                
                # Aligner-specific logs
                if [[ "$ALIGNER" == "bowtie2" ]]; then
                    echo -e "Bowtie2 does not generate standalone log files like STAR does.\nAll mapping statistics for this run are recorded in the main detailed_output.log." \
                        > "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/README_Bowtie2_LOGS.txt"
                else
                    # STAR logs
                    cp "$analysis_dir"/*Log.final.out "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/" 2>/dev/null
                    cp "$analysis_dir"/*Log.out "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/DETAILED_LOGS_CAN_BE_DELETED/" 2>/dev/null
                    cp "$analysis_dir"/*Log.progress.out "$OUTPUT_ROOT/$DIR_REPORTS/ALIGNER_LOGS/DETAILED_LOGS_CAN_BE_DELETED/" 2>/dev/null
                    cp "$analysis_dir"/*SJ.out.tab "$OUTPUT_ROOT/$DIR_OTHERS/STAR_OUTPUT/" 2>/dev/null
                fi
                
                # ncRNA Mapping outputs
                if [[ -d "$analysis_dir/OTHERS/ncRNA_Mapping" ]]; then
                    mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping"
                    cp "$analysis_dir"/OTHERS/ncRNA_Mapping/*_ncrna_stats.txt "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                    cp "$analysis_dir"/OTHERS/ncRNA_Mapping/*_ncrna.bam* "$OUTPUT_ROOT/$DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
                fi
                
                # CTK Analysis outputs
                if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
                    CTK_DEST="$OUTPUT_ROOT/5_CTK_Analysis"
                    mkdir -p "$CTK_DEST"
                    for ctk_dir in CTK_Analysis CIMS_Analysis CITS_Analysis; do
                        if [[ -d "$analysis_dir/$ctk_dir" ]]; then
                            mkdir -p "$CTK_DEST/$sample_name"
                            mv "$analysis_dir/$ctk_dir" "$CTK_DEST/$sample_name/" 2>/dev/null
                        fi
                    done
                fi

                # Pipeline Log (The child specific one)
                # It's named ${sample_name}_analysis.log inside the dir
                cp "$analysis_dir"/*.log "$OUTPUT_ROOT/$DIR_REPORTS/SAMPLES/${sample_name}_detailed.log" 2>/dev/null
                
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
    PEAK_SCRIPT="$(dirname "$0")/PEAKittyPeak.sh"
    
    if [[ -x "$PEAK_SCRIPT" ]]; then
        console_msg "  > Running HOMER Peak Calling (Aggregated)..."
        
        # Run it inside the output root
        # We pass "2_COLLAPSED_BED" which is relative to OUTPUT_ROOT
        curr_dir=$(pwd)
        cd "$OUTPUT_ROOT" || exit 1
        
        # Define separate log for peak calling details
        PEAK_LOG="REPORTS/PEAK/Combined_PeakCalling.log"
        console_msg "  > Detailed Peak Log: $PEAK_LOG"

        # Call script (using absolute path or relative to old pwd)
        # Call script with -n COMBINED to create COMBINED_peaks folder
        # Use -i 2_COLLAPSED_BED (relative to OUTPUT_ROOT) and --aggregate parameter
        PEAK_CMD="bash $PEAK_SCRIPT -i 2_COLLAPSED_BED --aggregate -n COMBINED -p $PEAK_DIST -z $PEAK_SIZE -f $FRAG_LEN"
        
        # Add --ctk-dir if CTK analysis was enabled
        if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
            # Determine CTK folder name
            if [[ "$RUN_CIMS" == "true" ]] && [[ "$RUN_CITS" == "true" ]]; then
                CTK_FOLDER="5_CTK_Analysis"
            elif [[ "$RUN_CIMS" == "true" ]]; then
                CTK_FOLDER="5_CIMS_Analysis"
            else
                CTK_FOLDER="5_CITS_Analysis"
            fi
            PEAK_CMD="$PEAK_CMD --ctk-dir ./$CTK_FOLDER"
            
            # Add CIMS/CITS params
            PEAK_CMD="$PEAK_CMD --cims-fdr $CIMS_FDR --cits-pval $CITS_PVALUE"
            if [[ -n "$GROUPS_FILE" ]]; then
                PEAK_CMD="$PEAK_CMD -g $GROUPS_FILE"
            fi
        fi
        
        # Force hardcoded path to ensure correct filename (bypass potential variable issues)
        eval "$PEAK_CMD" > "REPORTS/PEAK/Combined_PeakCalling.log" 2>&1
        
        # (No Symlink to remove)
        
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
    
    # Combined BedGraph Generation (MUST happen before add_matrix_columns)
    # Uses GROUPS_FILE if present, otherwise creates "All Samples"
    if [[ -n "$GROUPS_FILE" ]]; then
        run_combined_bedgraph "$OUTPUT_ROOT" "$GROUPS_FILE" "$OUTPUT_ROOT/$DIR_BG"
    fi
    
    # Add enhanced columns to peak matrix (after combined bedgraphs are ready)
    PEAK_MATRIX="$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/COMBINED_peakMatrix.txt"
    PEAKS_BED="$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/peaks_Sorted.bed"
    
    if [[ -f "$PEAK_MATRIX" && -f "$PEAKS_BED" ]]; then
        console_msg "  > Adding enhanced matrix columns..."
        add_matrix_columns "$PEAK_MATRIX" "$PEAKS_BED" \
            "$OUTPUT_ROOT/$DIR_BG" "$OUTPUT_ROOT/$DIR_BG/scale_factors.tsv" \
            "$GROUPS_FILE"
        
        # Cleanup intermediate peak files to save disk space
        rm -f "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/COMBINED.bed" \
              "$OUTPUT_ROOT/$DIR_PEAKS/COMBINED_PEAKS/peaks.bed" \
              "$PEAKS_BED"
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
    
    # Mapping Depth Summary (from scale_factors.tsv generated during bedgraph creation)
    SCALE_FILE="$OUTPUT_ROOT/$DIR_BG/scale_factors.tsv"
    if [[ -f "$SCALE_FILE" ]]; then
        console_msg "\n[MAPPING DEPTH SUMMARY]"
        printf "  %-25s %-15s %s\n" "Sample" "Mapped Reads" "Scale Factor"
        console_msg "  ----------------------------------------------------------------"
        while IFS=$'\t' read -r sample reads scale; do
            printf "  %-25s %-15s %s\n" "$sample" "$reads" "$scale"
        done < "$SCALE_FILE"
        console_msg "  ----------------------------------------------------------------"
    fi
    
    # Cleanup Analysis Folders and Temp Files or move to OTHERS if -k
    if [[ "$KEEP_INTERMEDIATE" == "yes" ]]; then
        console_msg "  > Moving sample analysis directories to $OUTPUT_ROOT/$DIR_OTHERS/sample_analysis/..."
        mkdir -p "$OUTPUT_ROOT/$DIR_OTHERS/sample_analysis"
        for analysis_dir in *_analysis; do
            [[ -d "$analysis_dir" ]] && mv "$analysis_dir" "$OUTPUT_ROOT/$DIR_OTHERS/sample_analysis/" 2>/dev/null
        done
    else
        console_msg "  > Cleaning up temporary analysis directories and files..."
        rm -rf *_analysis
    fi
    rm -f barcodes.fasta
    
    # Remove sampled input if it exists (pattern match)
    rm -f *_sampled_*.fastq.gz
    
    # Remove CTK temp directories (CITS.pl_*, CIMS.pl_*, tag2profile.pl_*)
    rm -rf CITS.pl_* CIMS.pl_* tag2profile.pl_* 2>/dev/null
    
    # Remove CTK temp files (*.tmp from perl scripts)
    rm -f *.tmp 2>/dev/null

    # Also remove the main log if it was created in this dir
    # We want to MOVE it to REPORTS/detailed_output.log
    if [[ -f "$LOG_FILE" ]]; then
         mv "$LOG_FILE" "$OUTPUT_ROOT/$DIR_REPORTS/detailed_output.log"
         LOG_FILE="$OUTPUT_ROOT/$DIR_REPORTS/detailed_output.log"
    fi
    
    # Output Summary
    console_msg "\n[OUTPUT]"
    console_msg "  All results saved to: $OUTPUT_ROOT/"
    console_msg "  ├── 0_DEMUX_FASTQ/"
    console_msg "  ├── 1_BAM/"
    console_msg "  ├── 2_COLLAPSED_BED/"
    console_msg "  ├── 3_BEDGRAPH/"
    console_msg "  ├── 4_PEAKS/"
    if [[ -n "$DIR_CTK" ]]; then
        console_msg "  ├── ${DIR_CTK}/"
    fi
    console_msg "  ├── ${DIR_OTHERS}/"
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
    
    # Console log already named and captured (see post-parsing initialization)
    if [[ -n "$TEMP_CONSOLE_LOG" && -f "$TEMP_CONSOLE_LOG" ]]; then
        console_msg "  > Console log: $TEMP_CONSOLE_LOG"
    fi
    
    # Cleanup temp gzip for plain fastq
    if [[ -n "$GLOBAL_GZIP_TMP" ]] && [[ -f "$GLOBAL_GZIP_TMP" ]]; then
        rm -f "$GLOBAL_GZIP_TMP"
        log_info "Removed temp gzip: $GLOBAL_GZIP_TMP"
    fi
    
    send_notification "CLIPittyClip: $(basename "$INPUT_FILE")" "Pipeline execution finished successfully. Duration: ${H}h ${M}m ${S}s"
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

# ── Deduplication ────────────────────────────────────────────────────────────────
# Runs based on DEDUP_MODE only. CHILD_MODE controls console output, not execution.
if [[ "$CHILD_MODE" != "true" ]]; then console_msg "\n[DEDUPLICATING]"; fi
if [[ "$DEDUP_MODE" == "true" ]]; then
    if [[ "$CHILD_MODE" != "true" ]]; then print_section_item "Deduplicating Reads (fastq2collapse.pl)"; fi
    DEDUP_OUT="$(pwd)/${BASENAME}_dedup.fastq.gz"
    if run_dedup "$INPUT_FILE" "$DEDUP_OUT"; then
        INPUT_FILE="$DEDUP_OUT"
        if [[ "$CHILD_MODE" != "true" ]]; then print_section_item "Deduplication Complete"; fi
    else
        log_warning "Deduplication failed. Using original input."
        rm -f "$DEDUP_OUT"
    fi
else
    if [[ "$CHILD_MODE" != "true" ]]; then print_section_item "Deduplication Disabled (--no-dedup)"; fi
fi

# ── Section Headers (non-child mode only) ─────────────────────────────────────
if [[ "$CHILD_MODE" != "true" ]]; then
    console_msg "\n[DEMULTIPLEXING]"
    print_section_item "No Barcode File Provided"
    print_section_item "Skipping Demultiplexing"

    console_msg "\n[ANALYSIS]"
    printf "   1/1  %-20s : " "$BASENAME"
fi

# Directory Setup
mkdir -p "${BASENAME}_analysis"
cd "${BASENAME}_analysis" || exit 1
LOG_FILE="${BASENAME}_analysis.log" # redirect log to inside analysis dir
log_info "Working directory: $(pwd)"

# 1. Preprocessing
run_fastp "$INPUT_FILE" "$BASENAME" "$UMI_LEN" "$ADAPTER_3" "$THREADS" "$SAMPLE_SIZE" "$ECLIP_MODE" "$BC_LEN" "$SPACER_LEN"

# 1b. ncRNA Pre-filtering (if enabled and index exists)
CLEANED_FASTQ="${BASENAME}_cleaned.fastq.gz"
if [[ "$FILTER_NCRNA" == "true" ]]; then
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

# 3a. Canonical chromosome filtering (default ON, disable with --no-chr-filter)
if [[ "$FILTER_CHR" == "true" ]]; then
    FILTERED_BAM="${MAPPED_PREFIX}.filtered.bam"
    filter_canonical_chromosomes "$BAM_FILE" "$FILTERED_BAM"
    BAM_FILE="$FILTERED_BAM"
fi

COLLAPSED_BED="${BASENAME}_collapsed.bed"
MUTATION_FILE="${BASENAME}_mutations.txt"

# Unified preprocessing: Always use parseAlignment.pl for consistent output
# This generates both tags.bed and mutations.txt (future-proofing for CIMS/CITS)
log_info "Unified preprocessing: samtools calmd → parseAlignment.pl → tag2collapse.pl"
check_dependency parseAlignment.pl

run_parse_alignment "${BAM_FILE}" "${BASENAME}_parsed.bed" "${MUTATION_FILE}" "$GENOME_INDEX"

run_collapse_pcr "${BASENAME}_parsed.bed" "${COLLAPSED_BED}" "${UMI_LEN}" "${DEDUP_MODE}"

# 4. Coverage Analysis (Bedgraph)
run_coverage "${COLLAPSED_BED}" "${BASENAME}" "$GENOME_INDEX/chrom.sizes" "${BAM_FILE}" # Pass genome file if available

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

# Finalize status line (prints "Done!" to end the Preprocessing > ... > Peaks > chain)
if [[ "$CHILD_MODE" != "true" ]]; then
    update_status_done
fi

# ── Organize output (non-child single-file mode only) ──────────────────────────
if [[ "$CHILD_MODE" != "true" ]]; then
    # We are in ${BASENAME}_analysis/ right now.
    # Create organized output root one level up alongside (not inside) _analysis/
    SINGLE_OUTPUT_ROOT="$(dirname "$(pwd)")/${BASENAME}_output"
    mkdir -p "$SINGLE_OUTPUT_ROOT"

    SF_DIR_BAM="1_BAM"
    SF_DIR_BED="2_COLLAPSED_BED"
    SF_DIR_BG="3_BEDGRAPH"
    SF_DIR_PEAKS="4_PEAKS"
    if [[ "$RUN_CIMS" == "true" ]] || [[ "$RUN_CITS" == "true" ]]; then
        SF_DIR_OTHERS="6_OTHERS"
    else
        SF_DIR_OTHERS="5_OTHERS"
    fi
    SF_DIR_REPORTS="REPORTS"

    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_BAM"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_BED"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_BG"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_OTHERS"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS"

    # BAM and index
    mv "${BASENAME}"*.bam "$SINGLE_OUTPUT_ROOT/$SF_DIR_BAM/" 2>/dev/null
    mv "${BASENAME}"*.bai "$SINGLE_OUTPUT_ROOT/$SF_DIR_BAM/" 2>/dev/null

    # Collapsed BED
    mv "${COLLAPSED_BED:-${BASENAME}_collapsed.bed}" "$SINGLE_OUTPUT_ROOT/$SF_DIR_BED/" 2>/dev/null
    mv "${MUTATION_FILE:-${BASENAME}_mutations.txt}" "$SINGLE_OUTPUT_ROOT/$SF_DIR_BED/" 2>/dev/null

    # Bedgraphs and scale factors
    mv "${BASENAME}"*.bedgraph "$SINGLE_OUTPUT_ROOT/$SF_DIR_BG/" 2>/dev/null
    mv scale_factors.tsv "$SINGLE_OUTPUT_ROOT/$SF_DIR_BG/" 2>/dev/null

    # Peaks (HOMER directory or CTK flat BED)
    PEAK_DIR_NAME="${BASENAME}_peaks"
    if [[ -d "$PEAK_DIR_NAME" ]]; then
        mv "$PEAK_DIR_NAME" "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/" 2>/dev/null
        mv "${PEAK_DIR_NAME}_homer.log" "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/" 2>/dev/null
        # Duplicate the unified BED format up to root for immediate visibility
        cp "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/$PEAK_DIR_NAME/peaks_Sorted.bed" "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/FINAL_PEAKS.bed" 2>/dev/null
    else
        mv "${BASENAME}_peaks_raw.bed" "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/FINAL_PEAKS.bed" 2>/dev/null
        mv "${BASENAME}_peaks_ctk.log" "$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/" 2>/dev/null
    fi
    
    # -------------------------------------------------------------------------
    # Single-File Peak Matrix Generation (Unification with Batch mode)
    # -------------------------------------------------------------------------
    sample_bed="$SINGLE_OUTPUT_ROOT/$SF_DIR_BED/${COLLAPSED_BED:-${BASENAME}_collapsed.bed}"
    final_peaks="$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/FINAL_PEAKS.bed"
    scale_tsv="$SINGLE_OUTPUT_ROOT/$SF_DIR_BG/scale_factors.tsv"
    
    if [[ -f "$final_peaks" && -f "$sample_bed" ]]; then
        console_msg "  > Generating coverage matrix for Single Sequence..."
        
        # 1. Raw tag coverage
        coverage_file="$SINGLE_OUTPUT_ROOT/$SF_DIR_PEAKS/FINAL_PEAK_MATRIX.txt"
        cp "$final_peaks" "$coverage_file"
        
        bedtools coverage -s -a "$final_peaks" -b "$sample_bed" > "temp_cov.txt"
        awk '{print $7}' "temp_cov.txt" > "col_count.txt"
        paste "$coverage_file" "col_count.txt" > "temp_paste.txt"
        mv "temp_paste.txt" "$coverage_file"
        
        s_name="${BASENAME}"
        HEADER_STR="chr\tstart\tend\tname\tscore\tstrand\tTC_${s_name}"
        echo -e "$HEADER_STR" > "colnames.txt"
        cat "colnames.txt" "$coverage_file" > "temp_final.txt"
        mv "temp_final.txt" "$coverage_file"
        
        rm -f "temp_cov.txt" "col_count.txt" "colnames.txt" 2>/dev/null
        
        # 2. Add Normalized TC columns
        if [[ -f "$scale_tsv" ]]; then
            console_msg "  > Normalizing single-sample tag counts..."
            add_matrix_columns "$coverage_file" "$final_peaks" "$SINGLE_OUTPUT_ROOT/$SF_DIR_BG" "$scale_tsv" ""
        else
            console_msg "  > Warning: scale_factors.tsv missing, skipping Normalization."
        fi
    fi

    # CTK Analysis (CIMS/CITS)
    for ctk_folder in "CTK_Analysis" "CIMS_Analysis" "CITS_Analysis"; do
        if [[ -d "$ctk_folder" ]]; then
            CTK_OUT_DIR=""
            if [[ "$ctk_folder" == "CTK_Analysis" ]]; then CTK_OUT_DIR="5_CTK_Analysis"; fi
            if [[ "$ctk_folder" == "CIMS_Analysis" ]]; then CTK_OUT_DIR="5_CIMS_Analysis"; fi
            if [[ "$ctk_folder" == "CITS_Analysis" ]]; then CTK_OUT_DIR="5_CITS_Analysis"; fi
            mkdir -p "$SINGLE_OUTPUT_ROOT/$CTK_OUT_DIR"
            cp -r "$ctk_folder/"* "$SINGLE_OUTPUT_ROOT/$CTK_OUT_DIR/" 2>/dev/null
            SF_DIR_CTK="$CTK_OUT_DIR"
        fi
    done

    # ncRNA mapping output
    if [[ -d "OTHERS/ncRNA_Mapping" ]]; then
        mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_OTHERS/ncRNA_Mapping"
        mv "OTHERS/ncRNA_Mapping/"* "$SINGLE_OUTPUT_ROOT/$SF_DIR_OTHERS/ncRNA_Mapping/" 2>/dev/null
    fi

    # Reports: fastp, aligner logs, analysis log
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/FASTP_REPORT"
    mkdir -p "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/ALIGNER_LOGS"
    mv "${BASENAME}"*_fastp.html "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
    mv "${BASENAME}"*_fastp.json "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/FASTP_REPORT/" 2>/dev/null
    
    if [[ "$ALIGNER" == "bowtie2" ]]; then
        echo -e "Bowtie2 does not generate standalone log files like STAR does.\nAll mapping statistics for this run are recorded in the main console_output.log." > "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/ALIGNER_LOGS/README_Bowtie2_LOGS.txt"
    else
        mv "${BASENAME}"*.Log.final.out "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/ALIGNER_LOGS/" 2>/dev/null
        mv "${BASENAME}"*.Log.out "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/ALIGNER_LOGS/" 2>/dev/null
    fi

    # Move analysis log into REPORTS
    if [[ -f "$LOG_FILE" ]]; then
        mv "$LOG_FILE" "$SINGLE_OUTPUT_ROOT/$SF_DIR_REPORTS/detailed_output.log"
    fi

    log_info "Output organized: $SINGLE_OUTPUT_ROOT"
fi

# Cleanup
if [[ "$KEEP_INTERMEDIATE" != "yes" ]]; then
    log_info "Cleaning up intermediate files..."
    if [[ -n "$SINGLE_OUTPUT_ROOT" ]]; then
        # Single mode cleans entire scratch directory natively
        cd ..
        rm -rf "${BASENAME}_analysis"
    else
        rm -f "${BASENAME}_cleaned.fastq.gz" "${BASENAME}_raw.bed" "${BASENAME}_parsed.bed"
    fi
else
    # Even if keeping intermediates, we must return to base dir in single-mode
    if [[ -n "$SINGLE_OUTPUT_ROOT" ]]; then
        cd ..
    fi
fi

if [[ -n "$GLOBAL_GZIP_TMP" ]] && [[ -f "$GLOBAL_GZIP_TMP" ]]; then
    rm -f "$GLOBAL_GZIP_TMP"
    log_info "Removed temp gzip: $GLOBAL_GZIP_TMP"
fi

log_info "Analysis Finished Successfully!"
log_info "Results in: ${SINGLE_OUTPUT_ROOT#../}"

# Calculate Duration
PIPELINE_END=$(date +%s)
DURATION=$((PIPELINE_END - PIPELINE_START))
H=$((DURATION/3600))
M=$(( (DURATION%3600)/60 ))
S=$((DURATION%60))

log_info "End Time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Total Duration: ${H}h ${M}m ${S}s"

# ncRNA Filtering Summary (non-child mode, single-file only)
if [[ "$CHILD_MODE" != "true" ]] && [[ "$FILTER_NCRNA" == "true" ]]; then
    ncrna_stats_sf="${SINGLE_OUTPUT_ROOT}/${SF_DIR_OTHERS}/ncRNA_Mapping/${BASENAME}_ncrna_stats.txt"
    if [[ -f "$ncrna_stats_sf" ]]; then
        console_msg "\n[ncRNA FILTERING SUMMARY]"
        printf "  %-25s %-15s %-15s %s\n" "Sample" "ncRNA Reads" "Total Reads" "% Filtered"
        console_msg "  ----------------------------------------------------------------"
        align_rate=$(grep "overall alignment rate" "$ncrna_stats_sf" | grep -oE "[0-9]+\.[0-9]+%" || echo "N/A")
        total=$(grep "reads; of these:" "$ncrna_stats_sf" | grep -oE "^[0-9]+" || echo "N/A")
        aligned=$(grep "aligned exactly 1 time" "$ncrna_stats_sf" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        multi=$(grep "aligned >1 times" "$ncrna_stats_sf" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        ncrna_reads=$(( ${aligned:-0} + ${multi:-0} ))
        printf "  %-25s %-15s %-15s %s\n" "$BASENAME" "$ncrna_reads" "${total:-N/A}" "${align_rate:-N/A}"
        console_msg "  ----------------------------------------------------------------"
    fi
fi

# Mapping Depth Summary (non-child mode, single-file only)
if [[ "$CHILD_MODE" != "true" ]]; then
    scale_sf="${SINGLE_OUTPUT_ROOT}/${SF_DIR_BG}/scale_factors.tsv"
    if [[ -f "$scale_sf" ]]; then
        console_msg "\n[MAPPING DEPTH SUMMARY]"
        printf "  %-25s %-15s %s\n" "Sample" "Mapped Reads" "Scale Factor"
        console_msg "  ----------------------------------------------------------------"
        while IFS=$'\t' read -r sample reads scale; do
            printf "  %-25s %-15s %s\n" "$sample" "$reads" "$scale"
        done < "$scale_sf"
        console_msg "  ----------------------------------------------------------------"
    fi
fi

# Console completion summary (non-child mode only)
if [[ "$CHILD_MODE" != "true" ]]; then
    console_msg "\n[OUTPUT]"
    console_msg "  All results saved to: $SINGLE_OUTPUT_ROOT/"
    console_msg "  ├── 1_BAM/"
    console_msg "  ├── 2_COLLAPSED_BED/"
    console_msg "  ├── 3_BEDGRAPH/"
    console_msg "  ├── 4_PEAKS/"
    if [[ -n "${SF_DIR_CTK:-}" ]]; then
        console_msg "  ├── ${SF_DIR_CTK}/"
    fi
    console_msg "  ├── ${SF_DIR_OTHERS}/"
    console_msg "  └── REPORTS/"

    console_msg "\n[SUCCESS] Pipeline finished."
    console_msg "End Time: $(date '+%Y-%m-%d %H:%M:%S')"
    console_msg "Total Duration: ${H}h ${M}m ${S}s"
    
    # Final step: capture the floating console_output.log 
    if [[ -f "$TEMP_CONSOLE_LOG" ]]; then
        if [[ -n "$SINGLE_OUTPUT_ROOT" ]]; then
            # Strip ../ from the relative path since we are back in root
            mv "$TEMP_CONSOLE_LOG" "${SINGLE_OUTPUT_ROOT#../}/REPORTS/console_output.log" 2>/dev/null
        fi
    fi
fi

send_notification "CLIPittyClip: $BASENAME" "Analysis finished successfully. Duration: ${H}h ${M}m ${S}s"
