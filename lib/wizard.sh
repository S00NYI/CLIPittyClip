#!/bin/bash

# wizard.sh - Interactive configuration for CLIPittyClip
# Part of CLIPittyClip v3.0

# Colors for Wizard
WIZ_CYAN='\033[1;36m'
WIZ_GREEN='\033[1;32m'
WIZ_YELLOW='\033[1;33m'
WIZ_RED='\033[0;31m'
WIZ_NC='\033[0m'

print_wiz_header() {
    clear
    echo -e "${WIZ_CYAN}==================================================================${WIZ_NC}"
    echo -e "${WIZ_CYAN}        CLIPittyClip: Advanced Configuration Wizard${WIZ_NC}"
    echo -e "${WIZ_CYAN}==================================================================${WIZ_NC}"
    echo "You have enabled advanced mode. You can override default parameters"
    echo "for each tool in the pipeline."
    echo ""
    echo -e "${WIZ_YELLOW}[INSTRUCTIONS]${WIZ_NC}"
    echo "- Press [ENTER] to keep the default settings."
    echo "- Format: --option1 value --option2 value --flag"
    echo "- Note: Not all options require values."
    echo ""
}

validate_option() {
    local tool="$1"
    local arg_string="$2"
    local valid=true
    
    if [[ -z "$arg_string" ]]; then return 0; fi
    
    # NEW CHECK: Input must look like flags (start with -)
    if [[ "$arg_string" != -* ]]; then
          echo -ne "\n"
          echo -e "  ${WIZ_RED}[ERROR] Invalid input. Options must start with '-' (e.g., -F or --option).${WIZ_NC}"
          return 1
    fi
    
    echo -ne "  ${WIZ_YELLOW}[VALIDATION] Checking options...${WIZ_NC}"
    
    # Get help text once
    local help_text=""
    if [[ "$tool" == "fastp" ]]; then help_text=$(fastp --help 2>&1); fi
    if [[ "$tool" == "STAR" ]]; then help_text=$(STAR --help 2>&1); fi
    if [[ "$tool" == "bowtie2" ]]; then help_text=$(bowtie2 --help 2>&1); fi
    if [[ "$tool" == "findPeaks" ]]; then help_text=$(findPeaks 2>&1); fi # findPeaks prints help on stderr
    
    # Parse args
    for word in $arg_string; do
        if [[ "$word" == -* ]]; then
             # It's a flag
             if [[ "$help_text" != *"$word"* ]]; then
                 echo ""
                 echo -e "  ${WIZ_RED}[ERROR] Option '$word' not found in $tool help text.${WIZ_NC}"
                 valid=false
             fi
        fi
    done
    
    if [[ "$valid" == "true" ]]; then
        echo -e "${WIZ_GREEN} OK.${WIZ_NC}"
        return 0
    else
        return 1
    fi
}

run_wizard_fastp() {
    local def_fastp="--length_required 16 --average_qual 30"
    
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo -e "${WIZ_CYAN}1. Preprocessing (fastp)${WIZ_NC}"
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo "Description: Adapters, Quality Trimming, UMI extraction."
    echo "Current Default: $def_fastp"
    echo ""
    echo -e "${WIZ_GREEN}[COMMON OPTIONS]${WIZ_NC}"
    echo "  -q, --qualified_quality_phred <int>   Quality threshold (Default: 15)"
    echo "  -l, --length_required <int>           Min read length (Default: 15)"
    echo "  --average_qual <int>                  Mean quality requirement (Default: 0)"
    echo "  --trim_front1 <int>                   Trim bases from front of read 1"
    echo "  --trim_tail1 <int>                    Trim bases from tail of read 1"
    echo "  --detect_adapter_for_pe               Auto-detect PE adapters"
    echo "  --adapter_sequence <str>              Manually specify adapter"
    echo -e "  --dedup                               ${WIZ_RED}[DISABLED] handled by seqkit${WIZ_NC}"
    echo ""
    
    while true; do
        read -p "Do you want to modify fastp parameters? [y/N]: " mod_fastp
        if [[ -z "$mod_fastp" ]]; then mod_fastp="n"; fi
        if [[ "$mod_fastp" =~ ^[YyNn]$ ]]; then break; fi
        echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
    done
    
    if [[ "$mod_fastp" =~ ^[Yy]$ ]]; then
        while true; do
            read -p "Enter additional fastp arguments: " ADV_FASTP_ARGS
            if validate_option "fastp" "$ADV_FASTP_ARGS"; then
                break
            else
                read -p "Retry? (Press Enter to skip, or y to retry): " retry
                if [[ "$retry" != "y" ]]; then ADV_FASTP_ARGS=""; break; fi
            fi
        done
        echo -e "${WIZ_GREEN}[UPDATED]${WIZ_NC} fastp will run with: $def_fastp $ADV_FASTP_ARGS"
    fi
}

run_wizard_mapping() {
    echo ""
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo -e "${WIZ_CYAN}2. Mapping (Select Strategy)${WIZ_NC}"
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    
    local chosen_aligner="star"
    
    while true; do
        echo "Which aligner do you want to use?"
        echo "  1) STAR (Default)"
        echo "  2) Bowtie2"
        read -p "Selection [1]: " align_sel
        
        # Default to 1
        if [[ -z "$align_sel" ]]; then align_sel="1"; fi
        
        if [[ "$align_sel" == "1" ]]; then
            chosen_aligner="star"
            break
        elif [[ "$align_sel" == "2" ]]; then
            chosen_aligner="bowtie2"
            break
        else
             echo -e "${WIZ_RED}[ERROR] Invalid selection. Please enter 1 or 2.${WIZ_NC}"
        fi
    done
    
    export ALIGNER="$chosen_aligner" # Export for main script
    
    # Defaults for Aligners
    local def_star="--outFilterMultimapNmax 10 --outFilterMismatchNmax ${MISMATCHES:-2} --alignEndsType EndToEnd --outSAMattributes ... MD"
    local def_bt2="--md --end-to-end (Standard Sensitivity)"
    
    echo ""
    # Capitalize for display
    local display_aligner="STAR"
    if [[ "$chosen_aligner" == "bowtie2" ]]; then display_aligner="Bowtie2"; fi
    
    echo -e "${WIZ_CYAN}2a. Mapping ($display_aligner)${WIZ_NC}"
    echo "Description: Genome Alignment."
    
    if [[ "$chosen_aligner" == "star" ]]; then
        echo "Current Default: $def_star"
        echo ""
        echo -e "${WIZ_GREEN}[STAR OPTIONS]${WIZ_NC}"
        echo "  --outFilterMultimapNmax <int>         Max mapped places allowed (Default: 10)"
        echo "  --outFilterMismatchNmax <int>         Max mismatches allowed (Default: 10)"
        echo "  --alignEndsType <str>                 EndToEnd (Default) or Local"
        echo "  --outFilterScoreMinOverLread <float>  Normalized score overhead (0.66)"
        echo "  --chimSegmentMin <int>                Min chimeric segment length (0 to disable)"
        
        while true; do
            read -p "Do you want to modify STAR parameters? [y/N]: " mod_star
            if [[ -z "$mod_star" ]]; then mod_star="n"; fi
            if [[ "$mod_star" =~ ^[YyNn]$ ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
        done
        
        if [[ "$mod_star" =~ ^[Yy]$ ]]; then
             while true; do
                read -p "Enter additional STAR arguments: " ADV_ALIGNER_ARGS
                if validate_option "STAR" "$ADV_ALIGNER_ARGS"; then break; else
                    read -p "Retry? (Press Enter to skip, or y to retry): " retry
                    if [[ "$retry" != "y" ]]; then ADV_ALIGNER_ARGS=""; break; fi
                fi
             done
        fi
    else
        # Bowtie2
        echo "Current Default: $def_bt2"
        echo ""
        echo -e "${WIZ_GREEN}[BOWTIE2 OPTIONS]${WIZ_NC}"
        echo "  --local                               Local alignment"
        echo "  --end-to-end                          End-to-End alignment (Default)"
        echo "  -N <int>                              Max mismatches in seed (0 or 1)"
        echo "  -L <int>                              Seed length"
        
        while true; do
            read -p "Do you want to modify Bowtie2 parameters? [y/N]: " mod_bt2
            if [[ -z "$mod_bt2" ]]; then mod_bt2="n"; fi
            if [[ "$mod_bt2" =~ ^[YyNn]$ ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
        done
        
        if [[ "$mod_bt2" =~ ^[Yy]$ ]]; then
             while true; do
                read -p "Enter additional Bowtie2 arguments: " ADV_ALIGNER_ARGS
                if validate_option "bowtie2" "$ADV_ALIGNER_ARGS"; then break; else
                    read -p "Retry? (Press Enter to skip, or y to retry): " retry
                    if [[ "$retry" != "y" ]]; then ADV_ALIGNER_ARGS=""; break; fi
                fi
             done
        fi
    fi
}

run_wizard_homer() {
    local def_homer="-style factor -L 2 -localSize 10000 -minDist ${PEAK_DIST:-50}"
    
    echo ""
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo -e "${WIZ_CYAN}3. Peak Calling (HOMER)${WIZ_NC}"
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo "Description: Peak Identification."
    echo "Current Default: $def_homer"
    echo ""
    echo -e "${WIZ_GREEN}[HOMER OPTIONS]${WIZ_NC}"
    echo "  -style <str>                          factor (TF), histone (Broad), groseq, tss"
    echo "  -F <float>                            Fold enrichment over control (Default: 4.0)"
    echo "  -L <float>                            Local filtering enrichment (Default: 4.0)"
    echo "  -localSize <int>                      Region size for local background (Default: 10000)"
    echo "  -size <int>                           Peak size (Default: auto/varies)"
    echo "  -minDist <int>                        Min distance between peaks"
    
    while true; do
        read -p "Do you want to modify HOMER parameters? [y/N]: " mod_homer
        if [[ -z "$mod_homer" ]]; then mod_homer="n"; fi
        if [[ "$mod_homer" =~ ^[YyNn]$ ]]; then break; fi
        echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
    done
    
    if [[ "$mod_homer" =~ ^[Yy]$ ]]; then
         while true; do
            read -p "Enter additional HOMER arguments: " ADV_HOMER_ARGS
            if validate_option "findPeaks" "$ADV_HOMER_ARGS"; then break; else
                read -p "Retry? (Press Enter to skip, or y to retry): " retry
                if [[ "$retry" != "y" ]]; then ADV_HOMER_ARGS=""; break; fi
            fi
         done
    fi
}

# CTK CIMS/CITS Analysis Configuration
run_wizard_ctk() {
    echo ""
    echo -e "${WIZ_CYAN}══════════════════════════════════════════════════════════════${WIZ_NC}"
    echo -e "${WIZ_CYAN}  CTK CIMS/CITS ANALYSIS${WIZ_NC}"
    echo -e "${WIZ_CYAN}══════════════════════════════════════════════════════════════${WIZ_NC}"
    echo ""
    echo "CIMS/CITS analysis identifies crosslink sites at single-nucleotide resolution."
    echo "  - CIMS: Detects crosslink-induced mutations (deletions/substitutions)"
    echo "  - CITS: Detects crosslink-induced truncations"
    echo ""
    
    # Enable CTK Analysis?
    while true; do
        read -p "Run CIMS/CITS analysis? [cims/cits/both/no] (default: no): " ctk_choice
        if [[ -z "$ctk_choice" ]]; then ctk_choice="no"; fi
        ctk_choice=$(echo "$ctk_choice" | tr '[:upper:]' '[:lower:]')
        if [[ "$ctk_choice" =~ ^(cims|cits|both|no)$ ]]; then break; fi
        echo -e "${WIZ_RED}[ERROR] Invalid input. Enter 'cims', 'cits', 'both', or 'no'.${WIZ_NC}"
    done
    
    # Set flags based on choice
    case "$ctk_choice" in
        cims)
            RUN_CIMS="true"
            RUN_CITS="false"
            ;;
        cits)
            RUN_CIMS="false"
            RUN_CITS="true"
            ;;
        both)
            RUN_CIMS="true"
            RUN_CITS="true"
            ;;
        no)
            RUN_CIMS="false"
            RUN_CITS="false"
            CIMS_ITERATIONS="10"
            CIMS_FDR="1"
            CITS_PVALUE="1"
            CITS_GAP="25"
            RUN_MOTIF="no"
            MOTIF_FLANK="10"
            return
            ;;
    esac
    
    # CIMS Parameters (only if CIMS enabled)
    if [[ "$RUN_CIMS" == "true" ]]; then
        echo ""
        echo -e "${WIZ_YELLOW}CIMS Parameters:${WIZ_NC}"
        
        # CIMS Iterations
        while true; do
            read -p "  Permutation iterations [default: 10]: " cims_iter
            if [[ -z "$cims_iter" ]]; then cims_iter="10"; break; fi
            if [[ "$cims_iter" =~ ^[0-9]+$ ]] && [[ "$cims_iter" -gt 0 ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Must be a positive integer.${WIZ_NC}"
        done
        CIMS_ITERATIONS="$cims_iter"
        
        # CIMS FDR
        while true; do
            read -p "  FDR threshold [default: 1 (all sites)]: " cims_fdr
            if [[ -z "$cims_fdr" ]]; then cims_fdr="1"; break; fi
            if [[ "$cims_fdr" =~ ^[0-9]*\.?[0-9]+$ ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Must be a decimal number (e.g., 1 or 0.001).${WIZ_NC}"
        done
        CIMS_FDR="$cims_fdr"
    fi
    
    # CITS Parameters (only if CITS enabled)
    if [[ "$RUN_CITS" == "true" ]]; then
        echo ""
        echo -e "${WIZ_YELLOW}CITS Parameters:${WIZ_NC}"
        
        # CITS P-value
        while true; do
            read -p "  P-value threshold [default: 1 (all sites)]: " cits_pval
            if [[ -z "$cits_pval" ]]; then cits_pval="1"; break; fi
            if [[ "$cits_pval" =~ ^[0-9]*\.?[0-9]+$ ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Must be a decimal number (e.g., 1 or 0.001).${WIZ_NC}"
        done
        CITS_PVALUE="$cits_pval"
        
        # CITS Gap
        while true; do
            read -p "  Clustering gap (-1 = no cluster) [default: 25]: " cits_gap
            if [[ -z "$cits_gap" ]]; then cits_gap="25"; break; fi
            if [[ "$cits_gap" =~ ^-?[0-9]+$ ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Must be an integer.${WIZ_NC}"
        done
        CITS_GAP="$cits_gap"
    fi
    
    echo ""
    echo -e "${WIZ_YELLOW}Motif Enrichment:${WIZ_NC}"
    
    # Run Motif Analysis?
    while true; do
        read -p "  Run HOMER motif enrichment? [Y/n]: " run_motif
        if [[ -z "$run_motif" ]]; then run_motif="y"; fi
        if [[ "$run_motif" =~ ^[YyNn]$ ]]; then break; fi
        echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
    done
    
    if [[ "$run_motif" =~ ^[Yy]$ ]]; then
        RUN_MOTIF="yes"
        
        # Motif Flank
        while true; do
            read -p "  Flanking nucleotides (±) [default: 10]: " motif_flank
            if [[ -z "$motif_flank" ]]; then motif_flank="10"; break; fi
            if [[ "$motif_flank" =~ ^[0-9]+$ ]] && [[ "$motif_flank" -gt 0 ]]; then break; fi
            echo -e "${WIZ_RED}[ERROR] Must be a positive integer.${WIZ_NC}"
        done
        MOTIF_FLANK="$motif_flank"
    else
        RUN_MOTIF="no"
        MOTIF_FLANK="10"
    fi
    
    echo ""
    echo -e "${WIZ_GREEN}CTK configuration complete.${WIZ_NC}"
}

# The main wizard entry point (sequential)
run_wizard() {
    print_wiz_header
    run_wizard_fastp
    run_wizard_mapping
    run_wizard_homer
    run_wizard_ctk
    
    # Summary
    echo ""
    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    echo -e "${WIZ_CYAN}[SUMMARY]${WIZ_NC}"
    
    echo "1. Preprocessing (Fastp)"
    echo -e "   Default:        --length_required 16 --average_qual 30"
    echo -e "   Added/Modified: ${ADV_FASTP_ARGS:-(None)}"
    
    local display_aligner="STAR"
    if [[ "$ALIGNER" == "bowtie2" ]]; then display_aligner="Bowtie2"; fi
    
    echo "2. Aligner ($display_aligner)"
    if [[ "$ALIGNER" == "star" ]]; then
        echo -e "   Default:        --outFilterMultimapNmax 10 --outFilterMismatchNmax ${MISMATCHES:-2} --alignEndsType EndToEnd --outSAMattributes ... MD"
    else
        echo -e "   Default:        --md --end-to-end (Standard Sensitivity)"
    fi
    echo -e "   Added/Modified: ${ADV_ALIGNER_ARGS:-(None)}"
    
    echo "3. Peak Calling (HOMER)"
    echo -e "   Default:        -style factor -L 2 -localSize 10000 -minDist ${PEAK_DIST:-50}"
    echo -e "   Added/Modified: ${ADV_HOMER_ARGS:-(None)}"
    
    echo "4. CTK CIMS/CITS Analysis"
    if [[ "$RUN_CIMS" == "true" && "$RUN_CITS" == "true" ]]; then
        echo -e "   Mode:           Both (CTK_Analysis)"
        echo -e "   CIMS:           iterations=${CIMS_ITERATIONS}, FDR=${CIMS_FDR}"
        echo -e "   CITS:           p-value=${CITS_PVALUE}, gap=${CITS_GAP}"
    elif [[ "$RUN_CIMS" == "true" ]]; then
        echo -e "   Mode:           CIMS only (CIMS_Analysis)"
        echo -e "   CIMS:           iterations=${CIMS_ITERATIONS}, FDR=${CIMS_FDR}"
    elif [[ "$RUN_CITS" == "true" ]]; then
        echo -e "   Mode:           CITS only (CITS_Analysis)"
        echo -e "   CITS:           p-value=${CITS_PVALUE}, gap=${CITS_GAP}"
    else
        echo -e "   Enabled:        No"
    fi
    if [[ "$RUN_CIMS" == "true" || "$RUN_CITS" == "true" ]]; then
        if [[ "$RUN_MOTIF" == "yes" ]]; then
            echo -e "   Motif:          HOMER (±${MOTIF_FLANK}nt flanks)"
        else
            echo -e "   Motif:          Disabled"
        fi
    fi

    echo -e "${WIZ_CYAN}------------------------------------------------------------------${WIZ_NC}"
    
    while true; do
        read -p "Proceed with these settings? [Y/n]: " confirm
        if [[ -z "$confirm" ]]; then confirm="y"; fi
        if [[ "$confirm" =~ ^[YyNn]$ ]]; then break; fi
        echo -e "${WIZ_RED}[ERROR] Invalid input. Please enter 'y' or 'n'.${WIZ_NC}"
    done
    
    if [[ "$confirm" =~ ^[Nn]$ ]]; then
        echo "Configuration aborted."
        exit 1
    fi
    
    # Persist all configuration
    echo "ADV_FASTP_ARGS=\"$ADV_FASTP_ARGS\"" > analysis_config.env
    echo "ALIGNER=\"$ALIGNER\"" >> analysis_config.env
    echo "ADV_ALIGNER_ARGS=\"$ADV_ALIGNER_ARGS\"" >> analysis_config.env
    echo "ADV_HOMER_ARGS=\"$ADV_HOMER_ARGS\"" >> analysis_config.env
    
    # CTK Configuration
    echo "RUN_CIMS=\"$RUN_CIMS\"" >> analysis_config.env
    echo "RUN_CITS=\"$RUN_CITS\"" >> analysis_config.env
    echo "CIMS_ITERATIONS=\"$CIMS_ITERATIONS\"" >> analysis_config.env
    echo "CIMS_FDR=\"$CIMS_FDR\"" >> analysis_config.env
    echo "CITS_PVALUE=\"$CITS_PVALUE\"" >> analysis_config.env
    echo "CITS_GAP=\"$CITS_GAP\"" >> analysis_config.env
    echo "RUN_MOTIF=\"$RUN_MOTIF\"" >> analysis_config.env
    echo "MOTIF_FLANK=\"$MOTIF_FLANK\"" >> analysis_config.env
    
    echo -e "${WIZ_GREEN}Configuration saved to 'analysis_config.env'. Starting Analysis...${WIZ_NC}"
    echo "=================================================================="
}

