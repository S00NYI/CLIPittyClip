#!/bin/bash

# wizard.sh - Interactive configuration wizard for CLIPittyClip Suite
# Part of CLIPittyClip v3.0

# Colors for Wizard
WIZ_CYAN='\033[1;36m'
WIZ_GREEN='\033[1;32m'
WIZ_YELLOW='\033[1;33m'
WIZ_RED='\033[0;31m'
WIZ_BOLD='\033[1m'
WIZ_NC='\033[0m'

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

print_wizard_header() {
    local title="$1"
    clear
    echo -e "${WIZ_CYAN}╔════════════════════════════════════════════════════════════════╗${WIZ_NC}"
    echo -e "${WIZ_CYAN}║${WIZ_NC}                    ${WIZ_BOLD}${title}${WIZ_NC}"
    echo -e "${WIZ_CYAN}╚════════════════════════════════════════════════════════════════╝${WIZ_NC}"
    echo ""
}

print_section() {
    local title="$1"
    echo ""
    echo -e "${WIZ_CYAN}═══════════════════════════════════════════════════════════════${WIZ_NC}"
    echo -e "${WIZ_BOLD}${title}${WIZ_NC}"
    echo -e "${WIZ_CYAN}═══════════════════════════════════════════════════════════════${WIZ_NC}"
}

print_doc_box() {
    local program="$1"
    shift
    local docs=("$@")
    
    echo ""
    echo -e "${WIZ_YELLOW}┌────────────────────────────────────────────────────────────────┐${WIZ_NC}"
    echo -e "${WIZ_YELLOW}│${WIZ_NC}                         ${WIZ_BOLD}DOCUMENTATION${WIZ_NC}                         ${WIZ_YELLOW}│${WIZ_NC}"
    echo -e "${WIZ_YELLOW}├────────────────────────────────────────────────────────────────┤${WIZ_NC}"
    echo -e "${WIZ_YELLOW}│${WIZ_NC} The prompts below show commonly modified options.              ${WIZ_YELLOW}│${WIZ_NC}"
    echo -e "${WIZ_YELLOW}│${WIZ_NC} Additional options can be passed via command line.             ${WIZ_YELLOW}│${WIZ_NC}"
    echo -e "${WIZ_YELLOW}│${WIZ_NC}                                                                ${WIZ_YELLOW}│${WIZ_NC}"
    echo -e "${WIZ_YELLOW}│${WIZ_NC} For full documentation, see:                                   ${WIZ_YELLOW}│${WIZ_NC}"
    for doc in "${docs[@]}"; do
        printf "${WIZ_YELLOW}│${WIZ_NC}  • %-60s ${WIZ_YELLOW}│${WIZ_NC}\n" "$doc"
    done
    echo -e "${WIZ_YELLOW}└────────────────────────────────────────────────────────────────┘${WIZ_NC}"
    echo ""
}

# Prompt for file path with validation
prompt_file_path() {
    local prompt_msg="$1"
    local var_name="$2"
    local result=""
    
    while true; do
        read -p "$prompt_msg" result
        if [[ -z "$result" ]]; then
            echo -e "${WIZ_RED}  ✗ Path cannot be empty. Please try again.${WIZ_NC}"
            continue
        fi
        
        # Expand ~ to home directory
        result="${result/#\~/$HOME}"
        
        if [[ -f "$result" ]]; then
            echo -e "${WIZ_GREEN}  ✓ File found: $result${WIZ_NC}"
            eval "$var_name=\"$result\""
            break
        else
            echo -e "${WIZ_RED}  ✗ File not found: $result${WIZ_NC}"
            read -p "  Try again? [Y/n]: " retry
            if [[ "$retry" =~ ^[Nn]$ ]]; then
                return 1
            fi
        fi
    done
    return 0
}

# Prompt for directory path with validation
prompt_dir_path() {
    local prompt_msg="$1"
    local var_name="$2"
    local check_type="${3:-any}"  # "star", "bowtie2", "bed", or "any"
    local result=""
    
    while true; do
        read -p "$prompt_msg" result
        if [[ -z "$result" ]]; then
            echo -e "${WIZ_RED}  ✗ Path cannot be empty. Please try again.${WIZ_NC}"
            continue
        fi
        
        # Expand ~ to home directory
        result="${result/#\~/$HOME}"
        
        if [[ -d "$result" ]]; then
            echo -e "${WIZ_GREEN}  ✓ Directory found: $result${WIZ_NC}"
            
            # Additional checks based on type
            if [[ "$check_type" == "bed" ]]; then
                local bed_count=$(ls "$result"/BED/*.bed 2>/dev/null | wc -l | tr -d ' ')
                if [[ "$bed_count" -eq 0 ]]; then
                    echo -e "${WIZ_RED}  ✗ No BED/ folder or .bed files found.${WIZ_NC}"
                    read -p "  Try again? [Y/n]: " retry
                    if [[ "$retry" =~ ^[Nn]$ ]]; then return 1; fi
                    continue
                fi
                echo -e "${WIZ_GREEN}  ✓ Found $bed_count .bed files in BED/ folder${WIZ_NC}"
            fi
            
            eval "$var_name=\"$result\""
            break
        else
            echo -e "${WIZ_RED}  ✗ Directory not found: $result${WIZ_NC}"
            read -p "  Try again? [Y/n]: " retry
            if [[ "$retry" =~ ^[Nn]$ ]]; then
                return 1
            fi
        fi
    done
    return 0
}

# Prompt for a value with default
prompt_value() {
    local prompt_msg="$1"
    local default="$2"
    local var_name="$3"
    local validation="${4:-any}"  # "int", "float", or "any"
    local result=""
    
    while true; do
        read -p "$prompt_msg (default: $default): " result
        if [[ -z "$result" ]]; then
            result="$default"
            break
        fi
        
        if [[ "$validation" == "int" ]]; then
            if [[ "$result" =~ ^[0-9]+$ ]]; then break; fi
            echo -e "${WIZ_RED}  Must be a positive integer.${WIZ_NC}"
        elif [[ "$validation" == "float" ]]; then
            if [[ "$result" =~ ^[0-9]*\.?[0-9]+$ ]]; then break; fi
            echo -e "${WIZ_RED}  Must be a number.${WIZ_NC}"
        else
            break
        fi
    done
    
    eval "$var_name=\"$result\""
}

# Prompt for yes/no with default
prompt_yesno() {
    local prompt_msg="$1"
    local default="$2"  # "y" or "n"
    local var_name="$3"
    local result=""
    
    local prompt_suffix="[y/N]"
    if [[ "$default" == "y" ]]; then prompt_suffix="[Y/n]"; fi
    
    while true; do
        read -p "$prompt_msg $prompt_suffix: " result
        if [[ -z "$result" ]]; then result="$default"; fi
        result=$(echo "$result" | tr '[:upper:]' '[:lower:]')
        if [[ "$result" =~ ^[yn]$ ]]; then
            eval "$var_name=\"$result\""
            break
        fi
        echo -e "${WIZ_RED}  Please enter 'y' or 'n'.${WIZ_NC}"
    done
}

# Prompt for selection from options
prompt_select() {
    local prompt_msg="$1"
    local var_name="$2"
    shift 2
    local options=("$@")
    local result=""
    
    echo "$prompt_msg"
    local i=1
    for opt in "${options[@]}"; do
        echo "  [$i] $opt"
        ((i++))
    done
    
    while true; do
        read -p "Selection: " result
        if [[ "$result" =~ ^[0-9]+$ ]] && [[ "$result" -ge 1 ]] && [[ "$result" -le "${#options[@]}" ]]; then
            eval "$var_name=\"$result\""
            break
        fi
        echo -e "${WIZ_RED}  Invalid selection. Please enter a number 1-${#options[@]}.${WIZ_NC}"
    done
}

# Print configuration summary
print_summary_box() {
    echo ""
    echo -e "${WIZ_CYAN}═══════════════════════════════════════════════════════════════${WIZ_NC}"
    echo -e "${WIZ_BOLD}CONFIGURATION SUMMARY${WIZ_NC}"
    echo -e "${WIZ_CYAN}═══════════════════════════════════════════════════════════════${WIZ_NC}"
}

# ═══════════════════════════════════════════════════════════════════════════════
# CLIPittyClip WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_clipittyclip() {
    print_wizard_header "CLIPittyClip Wizard"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Input Type
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 1: Input Type"
    echo ""
    prompt_select "Is your sample pooled/multiplexed or already demultiplexed?" WIZ_INPUT_TYPE \
        "Pooled/Multiplexed (requires demultiplexing)" \
        "Already demultiplexed"
    
    if [[ "$WIZ_INPUT_TYPE" == "1" ]]; then
        # Pooled sample
        echo ""
        prompt_file_path "  Enter path to pooled FASTQ file: " WIZ_INPUT_FILE || return 1
        prompt_file_path "  Enter path to barcode file: " WIZ_BARCODE_FILE || return 1
        WIZ_MODE="pooled"
    else
        # Demultiplexed
        echo ""
        prompt_select "Do you have a single file or multiple files?" WIZ_FILE_COUNT \
            "Single FASTQ file" \
            "Multiple files in a directory"
        
        if [[ "$WIZ_FILE_COUNT" == "1" ]]; then
            prompt_file_path "  Enter path to FASTQ file: " WIZ_INPUT_FILE || return 1
            WIZ_MODE="single"
        else
            prompt_dir_path "  Enter path to directory containing FASTQ files: " WIZ_INPUT_DIR || return 1
            WIZ_MODE="directory"
        fi
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Genome Index
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Genome Index"
    echo ""
    prompt_dir_path "  Enter path to genome index directory: " WIZ_GENOME_INDEX || return 1
    
    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Default or Advanced
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 3: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings" \
        "Modify program settings (advanced)"
    
    # Initialize defaults
    WIZ_ALIGNER="star"
    WIZ_THREADS="1"
    WIZ_UMI_LEN="7"
    WIZ_ADAPTER="L32"
    WIZ_CIMS="n"
    WIZ_CITS="n"
    WIZ_PEAK_DIST="50"
    WIZ_PEAK_SIZE="20"
    WIZ_FRAG_LEN="25"
    WIZ_HOMER_ARGS=""
    
    if [[ "$WIZ_SETTINGS_MODE" == "2" ]]; then
        # ─────────────────────────────────────────────────────────────────────
        # ADVANCED MODE
        # ─────────────────────────────────────────────────────────────────────
        print_doc_box "CLIPittyClip" \
            "CLIPittyClip.sh --help" \
            "STAR: https://github.com/alexdobin/STAR" \
            "Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2" \
            "CTK: https://zhanglab.c2b2.columbia.edu/" \
            "HOMER: http://homer.ucsd.edu/homer/"
        
        # Alignment Settings
        print_section "ALIGNMENT SETTINGS"
        echo ""
        prompt_select "Select aligner:" WIZ_ALIGNER_SEL "STAR (default)" "Bowtie2"
        if [[ "$WIZ_ALIGNER_SEL" == "2" ]]; then WIZ_ALIGNER="bowtie2"; fi
        
        prompt_value "  Enter number of threads" "1" WIZ_THREADS "int"
        
        echo ""
        if [[ "$WIZ_ALIGNER" == "star" ]]; then
            echo -e "  ${WIZ_YELLOW}Current STAR defaults:${WIZ_NC}"
            echo "    --outFilterMultimapNmax 10"
            echo "    --outFilterMismatchNmax 2"
            echo "    --alignEndsType EndToEnd"
            echo ""
            echo -e "  ${WIZ_GREEN}Common STAR options:${WIZ_NC}"
            echo "    --outFilterMultimapNmax <int>    Max mapped locations (default: 10)"
            echo "    --outFilterMismatchNmax <int>    Max mismatches (default: 2)"
            echo "    --alignEndsType <str>            EndToEnd or Local"
            echo "    --chimSegmentMin <int>           Min chimeric segment length"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See STAR manual for full options."
            echo "  Enter additional STAR arguments (optional):"
            read -p "  > " WIZ_ALIGNER_ARGS
        else
            echo -e "  ${WIZ_YELLOW}Current Bowtie2 defaults:${WIZ_NC}"
            echo "    --end-to-end (Standard sensitivity)"
            echo "    --md"
            echo ""
            echo -e "  ${WIZ_GREEN}Common Bowtie2 options:${WIZ_NC}"
            echo "    --local                          Local alignment mode"
            echo "    --very-sensitive                 More accurate alignment"
            echo "    -N <int>                         Max mismatches in seed (0 or 1)"
            echo "    -L <int>                         Seed length"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See Bowtie2 manual for full options."
            echo "  Enter additional Bowtie2 arguments (optional):"
            read -p "  > " WIZ_ALIGNER_ARGS
        fi
        
        # Preprocessing Settings
        print_section "PREPROCESSING SETTINGS (fastp)"
        echo ""
        prompt_value "  Enter UMI length" "7" WIZ_UMI_LEN "int"
        echo "  Adapter options: L32 (default), L19, or custom sequence"
        prompt_value "  Enter 3' adapter" "L32" WIZ_ADAPTER
        
        echo ""
        echo -e "  ${WIZ_YELLOW}Current fastp defaults:${WIZ_NC}"
        echo "    --length_required 16"
        echo "    --average_qual 30"
        echo ""
        echo -e "  ${WIZ_GREEN}Common fastp options:${WIZ_NC}"
        echo "    -q, --qualified_quality_phred <int>   Quality threshold (default: 15)"
        echo "    -l, --length_required <int>           Min read length (default: 15)"
        echo "    --average_qual <int>                  Mean quality requirement"
        echo "    --trim_front1 <int>                   Trim bases from front of read"
        echo "    --trim_tail1 <int>                    Trim bases from tail of read"
        echo ""
        echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See fastp documentation for full options."
        echo "  Enter additional fastp arguments (optional):"
        read -p "  > " WIZ_FASTP_ARGS
        
        # CIMS/CITS Settings
        print_section "CIMS/CITS SETTINGS (CTK)"
        echo ""
        prompt_yesno "  Enable CIMS analysis?" "n" WIZ_CIMS
        prompt_yesno "  Enable CITS analysis?" "n" WIZ_CITS
        
        if [[ "$WIZ_CIMS" == "y" ]] || [[ "$WIZ_CITS" == "y" ]]; then
            echo ""
            echo -e "  ${WIZ_YELLOW}Current CTK defaults:${WIZ_NC}"
            if [[ "$WIZ_CIMS" == "y" ]]; then
                echo "    CIMS: iterations=10, FDR=1 (all sites)"
            fi
            if [[ "$WIZ_CITS" == "y" ]]; then
                echo "    CITS: p-value=1 (all sites), gap=25"
            fi
            echo ""
            echo -e "  ${WIZ_GREEN}Common CTK options (pass via CLI):${WIZ_NC}"
            echo "    --cims-iter <int>      Permutation iterations (default: 10)"
            echo "    --cims-fdr <float>     FDR cutoff (e.g., 0.001)"
            echo "    --cits-pval <float>    P-value cutoff (e.g., 0.001)"
            echo "    --cits-gap <int>       Clustering gap (-1 = no clustering)"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See CTK documentation for full options."
            echo "  Enter additional CTK arguments (optional):"
            read -p "  > " WIZ_CTK_ARGS
        fi
        
        # Peak Calling Settings
        print_section "PEAK CALLING SETTINGS"
        echo ""
        echo "  [1] HOMER findPeaks (Default)"
        echo "  [2] CTK tag2peak.pl"
        echo ""
        prompt_value "  Select peak caller [1-2]" "1" caller_sel "int"
        
        if [[ "$caller_sel" == "2" ]]; then
            WIZ_PEAK_CALLER="ctk"
            echo -e "  ${WIZ_YELLOW}Current CTK tag2peak defaults:${WIZ_NC}"
            echo "    -big -ss --valley-seeking -minPH 2"
            echo ""
            prompt_value "  Enter gap (-gap)" "50" WIZ_PEAK_DIST "int"
            echo ""
            echo "  Enter additional CTK tag2peak arguments (optional):"
            read -p "  > " WIZ_CTK_PEAK_ARGS
        else
            WIZ_PEAK_CALLER="homer"
            echo -e "  ${WIZ_YELLOW}Current HOMER defaults:${WIZ_NC}"
            echo "    -style factor"
            echo "    -L 2"
            echo "    -localSize 10000"
            echo "    -minDist 50"
            echo "    -size 20"
            echo "    -fragLength 25"
            echo ""
            prompt_value "  Enter min distance between peaks" "50" WIZ_PEAK_DIST "int"
            prompt_value "  Enter peak size" "20" WIZ_PEAK_SIZE "int"
            prompt_value "  Enter fragment length" "25" WIZ_FRAG_LEN "int"
            echo ""
            echo -e "  ${WIZ_GREEN}Common HOMER findPeaks options:${WIZ_NC}"
            echo "    -style <str>           factor (TF), histone (broad), groseq, tss"
            echo "    -F <float>             Fold enrichment over control (default: 4.0)"
            echo "    -L <float>             Local filtering enrichment (default: 2.0)"
            echo "    -localSize <int>       Region size for local background"
            echo "    -strand <str>          separate or both"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See HOMER documentation for full options."
            echo "  Enter additional HOMER findPeaks arguments (optional):"
            read -p "  > " WIZ_HOMER_ARGS
        fi
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # CONFIGURATION SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    print_summary_box
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────┐"
    if [[ "$WIZ_MODE" == "pooled" ]]; then
        printf "  │ %-20s %-38s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
        printf "  │ %-20s %-38s │\n" "Barcodes:" "$(basename "$WIZ_BARCODE_FILE")"
    elif [[ "$WIZ_MODE" == "single" ]]; then
        printf "  │ %-20s %-38s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
    else
        printf "  │ %-20s %-38s │\n" "Input Dir:" "$WIZ_INPUT_DIR"
    fi
    printf "  │ %-20s %-38s │\n" "Genome Index:" "$(basename "$WIZ_GENOME_INDEX")"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Aligner:" "$(echo "$WIZ_ALIGNER" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Threads:" "$WIZ_THREADS"
    if [[ -n "$WIZ_ALIGNER_ARGS" ]]; then
        printf "  │ %-20s %-38s │\n" "Aligner Args:" "$WIZ_ALIGNER_ARGS"
    fi
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "UMI Length:" "$WIZ_UMI_LEN"
    printf "  │ %-20s %-38s │\n" "Adapter:" "$WIZ_ADAPTER"
    if [[ -n "$WIZ_FASTP_ARGS" ]]; then
        printf "  │ %-20s %-38s │\n" "fastp Args:" "$WIZ_FASTP_ARGS"
    fi
    echo "  ├─────────────────────────────────────────────────────────────┤"
    local cims_status="Disabled"; [[ "$WIZ_CIMS" == "y" ]] && cims_status="Enabled"
    local cits_status="Disabled"; [[ "$WIZ_CITS" == "y" ]] && cits_status="Enabled"
    printf "  │ %-20s %-38s │\n" "CIMS:" "$cims_status"
    printf "  │ %-20s %-38s │\n" "CITS:" "$cits_status"
    if [[ -n "$WIZ_CTK_ARGS" ]]; then
        printf "  │ %-20s %-38s │\n" "CTK Args:" "$WIZ_CTK_ARGS"
    fi
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Peak Caller:" "$(echo "${WIZ_PEAK_CALLER:-homer}" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Peak Distance/Gap:" "$WIZ_PEAK_DIST"
    if [[ "${WIZ_PEAK_CALLER:-homer}" == "homer" ]]; then
        printf "  │ %-20s %-38s │\n" "Peak Size:" "$WIZ_PEAK_SIZE"
        printf "  │ %-20s %-38s │\n" "Fragment Length:" "$WIZ_FRAG_LEN"
        if [[ -n "$WIZ_HOMER_ARGS" ]]; then
            printf "  │ %-20s %-38s │\n" "HOMER Args:" "$WIZ_HOMER_ARGS"
        fi
    else
        if [[ -n "$WIZ_CTK_PEAK_ARGS" ]]; then
            printf "  │ %-20s %-38s │\n" "CTK Peak Args:" "$WIZ_CTK_PEAK_ARGS"
        fi
    fi
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""
    
    prompt_yesno "  Start analysis with these settings?" "y" WIZ_CONFIRM
    
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi
    
    # Export variables for main script
    export WIZ_MODE WIZ_INPUT_FILE WIZ_INPUT_DIR WIZ_BARCODE_FILE WIZ_GENOME_INDEX
    export WIZ_ALIGNER WIZ_THREADS WIZ_UMI_LEN WIZ_ADAPTER
    export WIZ_ALIGNER_ARGS WIZ_FASTP_ARGS WIZ_CTK_ARGS WIZ_HOMER_ARGS
    export WIZ_CIMS WIZ_CITS
    export WIZ_PEAK_CALLER WIZ_PEAK_DIST WIZ_PEAK_SIZE WIZ_FRAG_LEN WIZ_CTK_PEAK_ARGS
    
    echo -e "${WIZ_GREEN}Starting analysis...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAPittyMap WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_mapittymap() {
    print_wizard_header "MAPittyMap Wizard"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Input File
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 1: Input File"
    echo ""
    prompt_file_path "  Enter path to input FASTQ file: " WIZ_INPUT_FILE || return 1
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Genome Index
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Genome Index"
    echo ""
    prompt_dir_path "  Enter path to genome index directory: " WIZ_GENOME_INDEX || return 1
    
    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Default or Advanced
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 3: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings" \
        "Modify program settings (advanced)"
    
    # Initialize defaults
    WIZ_ALIGNER="star"
    WIZ_THREADS="1"
    WIZ_ALIGN_MISMATCHES="2"
    WIZ_OUTPUT_NAME=""
    
    if [[ "$WIZ_SETTINGS_MODE" == "2" ]]; then
        # ─────────────────────────────────────────────────────────────────────
        # ADVANCED MODE
        # ─────────────────────────────────────────────────────────────────────
        print_doc_box "MAPittyMap" \
            "MAPittyMap.sh --help" \
            "STAR: https://github.com/alexdobin/STAR" \
            "Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2"
        
        print_section "ALIGNER SETTINGS"
        echo ""
        prompt_select "Select aligner:" WIZ_ALIGNER_SEL "STAR (default)" "Bowtie2"
        if [[ "$WIZ_ALIGNER_SEL" == "2" ]]; then WIZ_ALIGNER="bowtie2"; fi
        
        prompt_value "  Enter number of threads" "1" WIZ_THREADS "int"
        prompt_value "  Enter max alignment mismatches" "2" WIZ_ALIGN_MISMATCHES "int"
        
        local default_name=$(basename "$WIZ_INPUT_FILE" .fastq.gz)
        prompt_value "  Enter output name" "$default_name" WIZ_OUTPUT_NAME
        
        echo ""
        if [[ "$WIZ_ALIGNER" == "star" ]]; then
            echo -e "  ${WIZ_YELLOW}Current STAR defaults:${WIZ_NC}"
            echo "    --outFilterMultimapNmax 10"
            echo "    --outFilterMismatchNmax $WIZ_ALIGN_MISMATCHES"
            echo "    --alignEndsType EndToEnd"
            echo ""
            echo -e "  ${WIZ_GREEN}Common STAR options:${WIZ_NC}"
            echo "    --outFilterMultimapNmax <int>    Max mapped locations (default: 10)"
            echo "    --outFilterMismatchNmax <int>    Max mismatches"
            echo "    --alignEndsType <str>            EndToEnd or Local"
            echo "    --chimSegmentMin <int>           Min chimeric segment length"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See STAR manual for full options."
            echo "  Enter additional STAR arguments (optional):"
            read -p "  > " WIZ_ALIGNER_ARGS
        else
            echo -e "  ${WIZ_YELLOW}Current Bowtie2 defaults:${WIZ_NC}"
            echo "    --end-to-end (Standard sensitivity)"
            echo "    --md"
            echo ""
            echo -e "  ${WIZ_GREEN}Common Bowtie2 options:${WIZ_NC}"
            echo "    --local                          Local alignment mode"
            echo "    --very-sensitive                 More accurate alignment"
            echo "    -N <int>                         Max mismatches in seed (0 or 1)"
            echo "    -L <int>                         Seed length"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See Bowtie2 manual for full options."
            echo "  Enter additional Bowtie2 arguments (optional):"
            read -p "  > " WIZ_ALIGNER_ARGS
        fi
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # CONFIGURATION SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    print_summary_box
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────┐"
    printf "  │ %-20s %-38s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
    printf "  │ %-20s %-38s │\n" "Genome Index:" "$(basename "$WIZ_GENOME_INDEX")"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Aligner:" "$(echo "$WIZ_ALIGNER" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Threads:" "$WIZ_THREADS"
    printf "  │ %-20s %-38s │\n" "Max Mismatches:" "$WIZ_ALIGN_MISMATCHES"
    if [[ -n "$WIZ_OUTPUT_NAME" ]]; then
        printf "  │ %-20s %-38s │\n" "Output Name:" "$WIZ_OUTPUT_NAME"
    fi
    if [[ -n "$WIZ_ALIGNER_ARGS" ]]; then
        printf "  │ %-20s %-38s │\n" "Aligner Args:" "$WIZ_ALIGNER_ARGS"
    fi
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""
    
    prompt_yesno "  Start mapping with these settings?" "y" WIZ_CONFIRM
    
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi
    
    # Export variables for main script
    export WIZ_INPUT_FILE WIZ_GENOME_INDEX
    export WIZ_ALIGNER WIZ_THREADS WIZ_ALIGN_MISMATCHES WIZ_OUTPUT_NAME WIZ_ALIGNER_ARGS
    
    echo -e "${WIZ_GREEN}Starting mapping...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# PEAKittyPeak WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_peakittypeak() {
    print_wizard_header "PEAKittyPeak Wizard"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Working Directory
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 1: Input Directory"
    echo ""
    echo "  This tool requires a BED/ folder containing collapsed BED files."
    echo ""
    read -p "  Enter path to working directory (or press Enter for current): " WIZ_WORK_DIR
    if [[ -z "$WIZ_WORK_DIR" ]]; then WIZ_WORK_DIR="."; fi
    
    # Expand ~ and validate
    WIZ_WORK_DIR="${WIZ_WORK_DIR/#\~/$HOME}"
    
    if [[ ! -d "$WIZ_WORK_DIR" ]]; then
        echo -e "${WIZ_RED}  ✗ Directory not found: $WIZ_WORK_DIR${WIZ_NC}"
        return 1
    fi
    
    # Check for BED folder
    WIZ_BED_COUNT=$(ls "$WIZ_WORK_DIR"/BED/*.bed 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$WIZ_BED_COUNT" -eq 0 ]]; then
        echo -e "${WIZ_RED}  ✗ No BED/ folder or .bed files found in $WIZ_WORK_DIR${WIZ_NC}"
        return 1
    fi
    
    echo -e "${WIZ_GREEN}  ✓ Found $WIZ_BED_COUNT .bed files in BED/ folder${WIZ_NC}"
    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Default or Advanced
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings" \
        "Modify program settings (advanced)"
    
    # Initialize defaults
    WIZ_PEAK_DIST="50"
    WIZ_PEAK_SIZE="20"
    WIZ_FRAG_LEN="25"
    WIZ_OUTPUT_NAME="Combined"
    WIZ_HOMER_ARGS=""
    
    if [[ "$WIZ_SETTINGS_MODE" == "2" ]]; then
        # ─────────────────────────────────────────────────────────────────────
        # ADVANCED MODE
        # ─────────────────────────────────────────────────────────────────────
        print_doc_box "PEAKittyPeak" \
            "PEAKittyPeak.sh --help" \
            "HOMER findPeaks: http://homer.ucsd.edu/homer/ngs/peaks.html" \
            "HOMER Manual: http://homer.ucsd.edu/homer/"
        
        print_section "PEAK PARAMETERS"
        echo ""
        echo "  [1] HOMER findPeaks (Default)"
        echo "  [2] CTK tag2peak.pl"
        echo ""
        prompt_value "  Select peak caller [1-2]" "1" caller_sel "int"
        
        if [[ "$caller_sel" == "2" ]]; then
            WIZ_PEAK_CALLER="ctk"
            echo -e "  ${WIZ_YELLOW}Current CTK tag2peak defaults:${WIZ_NC}"
            echo "    -big -ss --valley-seeking -minPH 2"
            echo ""
            prompt_value "  Enter gap (-gap)" "50" WIZ_PEAK_DIST "int"
            prompt_value "  Enter output name" "Combined" WIZ_OUTPUT_NAME
            echo ""
            echo "  Enter additional CTK tag2peak arguments (optional):"
            read -p "  > " WIZ_CTK_PEAK_ARGS
        else
            WIZ_PEAK_CALLER="homer"
            echo -e "  ${WIZ_YELLOW}Current HOMER defaults:${WIZ_NC}"
            echo "    -style factor"
            echo "    -L 2"
            echo "    -localSize 10000"
            echo "    -minDist 50"
            echo "    -size 20"
            echo "    -fragLength 25"
            echo ""
            prompt_value "  Enter min distance between peaks" "50" WIZ_PEAK_DIST "int"
            prompt_value "  Enter peak size" "20" WIZ_PEAK_SIZE "int"
            prompt_value "  Enter fragment length" "25" WIZ_FRAG_LEN "int"
            prompt_value "  Enter output name" "Combined" WIZ_OUTPUT_NAME
            
            echo ""
            echo -e "  ${WIZ_GREEN}Common HOMER findPeaks options:${WIZ_NC}"
            echo "    -style <str>           factor (TF), histone (broad), groseq, tss"
            echo "    -F <float>             Fold enrichment over control (default: 4.0)"
            echo "    -L <float>             Local filtering enrichment (default: 2.0)"
            echo "    -localSize <int>       Region size for local background"
            echo "    -strand <str>          separate or both"
            echo ""
            echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} These are examples. See HOMER documentation for full options."
            echo "  Enter additional HOMER findPeaks arguments (optional):"
            read -p "  > " WIZ_HOMER_ARGS
        fi
    fi
    
    # ─────────────────────────────────────────────────────────────────────────
    # CONFIGURATION SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    print_summary_box
    echo ""
    echo "  ┌─────────────────────────────────────────────────────────────┐"
    printf "  │ %-20s %-38s │\n" "Working Dir:" "$WIZ_WORK_DIR"
    printf "  │ %-20s %-38s │\n" "BED Files:" "$WIZ_BED_COUNT files found"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Peak Caller:" "$(echo "${WIZ_PEAK_CALLER:-homer}" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Peak Distance/Gap:" "$WIZ_PEAK_DIST"
    if [[ "${WIZ_PEAK_CALLER:-homer}" == "homer" ]]; then
        printf "  │ %-20s %-38s │\n" "Peak Size:" "$WIZ_PEAK_SIZE"
        printf "  │ %-20s %-38s │\n" "Fragment Length:" "$WIZ_FRAG_LEN"
        if [[ -n "$WIZ_HOMER_ARGS" ]]; then
            printf "  │ %-20s %-38s │\n" "HOMER Args:" "$WIZ_HOMER_ARGS"
        fi
    else
        if [[ -n "$WIZ_CTK_PEAK_ARGS" ]]; then
            printf "  │ %-20s %-38s │\n" "CTK Peak Args:" "$WIZ_CTK_PEAK_ARGS"
        fi
    fi
    printf "  │ %-20s %-38s │\n" "Output Name:" "$WIZ_OUTPUT_NAME"
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""
    
    prompt_yesno "  Start peak calling with these settings?" "y" WIZ_CONFIRM
    
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi
    
    # Export variables for main script
    export WIZ_WORK_DIR WIZ_BED_COUNT WIZ_PEAK_CALLER WIZ_CTK_PEAK_ARGS
    export WIZ_PEAK_DIST WIZ_PEAK_SIZE WIZ_FRAG_LEN WIZ_OUTPUT_NAME WIZ_HOMER_ARGS
    
    echo -e "${WIZ_GREEN}Starting peak calling...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# LEGACY WIZARD FUNCTIONS (for backward compatibility)
# ═══════════════════════════════════════════════════════════════════════════════

# Keep old function names for backward compatibility
print_wiz_header() {
    print_wizard_header "$1"
}

run_wizard_mapping() {
    # Legacy function - redirect to new wizard
    run_wizard_mapittymap
}

run_wizard_homer() {
    # Legacy function - minimal implementation for standalone HOMER config
    prompt_value "Peak distance" "50" PEAK_DIST "int"
    prompt_value "Peak size" "20" PEAK_SIZE "int"
    prompt_value "Fragment length" "25" FRAG_LEN "int"
    read -p "Additional HOMER arguments: " ADV_HOMER_ARGS
}
