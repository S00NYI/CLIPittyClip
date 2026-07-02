#!/bin/bash

# wizard.sh - Interactive configuration wizard for CLIPittyClip Suite
# Part of CLIPittyClip v3.5

# Colors for Wizard
WIZ_CYAN='\033[1;36m'
WIZ_GREEN='\033[1;32m'
WIZ_YELLOW='\033[1;33m'
WIZ_RED='\033[0;31m'
WIZ_BOLD='\033[1m'
WIZ_NC='\033[0m'

WIZ_INTRO_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/wizard_intro.py"

_play_wizard_intro() {
    [[ "$WIZ_INTRO_PLAYED" == "1" ]] && return 0
    WIZ_INTRO_PLAYED=1
    python3 "$WIZ_INTRO_PY" || true
}

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
        read -p "$prompt_msg" result || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
            read -p "  Try again? [Y/n]: " retry || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
        read -p "$prompt_msg" result || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
                    read -p "  Try again? [Y/n]: " retry || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
                    if [[ "$retry" =~ ^[Nn]$ ]]; then return 1; fi
                    continue
                fi
                echo -e "${WIZ_GREEN}  ✓ Found $bed_count .bed files in BED/ folder${WIZ_NC}"
            fi
            
            eval "$var_name=\"$result\""
            break
        else
            echo -e "${WIZ_RED}  ✗ Directory not found: $result${WIZ_NC}"
            read -p "  Try again? [Y/n]: " retry || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
        read -p "$prompt_msg (default: $default): " result || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
        read -p "$prompt_msg $prompt_suffix: " result || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
        read -p "Selection: " result || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
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
# TOOL TRIAGE DISPATCHER
# ═══════════════════════════════════════════════════════════════════════════════
# Top-level entry point: ask which tool to configure, dispatch to its wizard.
# Invalid input is re-prompted by prompt_select. Ctrl-C exits normally.

run_wizard_dispatcher() {
    _play_wizard_intro
    print_wizard_header "CLIPittyClip Suite Wizard"

    echo "  Which tool would you like to configure?"
    echo ""
    prompt_select "Select tool:" WIZ_TOOL_CHOICE \
        "CLIPittyClip   — full CLIP-seq pipeline (preprocess → map → peak call)" \
        "PREPittyPrep   — preprocessing only (groomed *_prepped.fastq.gz)" \
        "MAPittyMap     — standalone mapping module" \
        "PEAKittyPeak   — standalone peak calling on existing BED files"

    case "$WIZ_TOOL_CHOICE" in
        1) run_wizard_clipittyclip ;;
        2) run_wizard_prepittyprep ;;
        3) run_wizard_mapittymap ;;
        4) run_wizard_peakittypeak ;;
        *)
            echo -e "${WIZ_RED}Internal error: unexpected tool choice '${WIZ_TOOL_CHOICE}'${WIZ_NC}" >&2
            return 1
            ;;
    esac
}

# ═══════════════════════════════════════════════════════════════════════════════
# CLIPittyClip KNOB REGISTRY (single source of truth for review/edit/diff/cmdline)
# ═══════════════════════════════════════════════════════════════════════════════
# Record format (pipe-separated):
#   category|key|label|var_name|type|default|help_text|cli_flag
# type: int | float | str | yn | select:opt1,opt2
# cli_flag: long form (e.g. --umi-length, --no-dedup); empty = no CLI emission

_cc_knob_registry() {
    cat <<'EOF'
1|threads|Threads|WIZ_THREADS|int|1|Number of parallel threads|-t
1|umi_len|UMI Length|WIZ_UMI_LEN|int|0|UMI length in bp (0 if none)|-u
1|adapter|3' Adapter|WIZ_ADAPTER|str|L32|3' adapter sequence (L32 default)|-a
1|bc_len|Barcode Length|WIZ_BC_LEN|str||Barcode length (blank=auto-detect)|--bc-len
1|spacer_len|Spacer Length|WIZ_SPACER_LEN|int|0|Spacer length trimmed after barcode|--spacer-len
1|bc_first|Barcode First|WIZ_BC_FIRST|yn|false|Layout [BC][UMI] not [UMI][BC]|--bc-first
1|min_qual|fastp Min Quality|WIZ_MIN_QUAL|int|30|fastp average quality threshold|--min-qual
1|no_dedup|Disable Dedup|WIZ_NO_DEDUP|yn|false|Skip PCR deduplication|--no-dedup
1|align_mismatches|Align Mismatches|WIZ_ALIGN_MISMATCHES|int|2|Max STAR alignment mismatches|--align-mismatches
1|demux_mismatches|Demux Mismatches|WIZ_DEMUX_MISMATCHES|int|1|Max barcode mismatches during demux|--demux-mismatches
1|aligner|Aligner|WIZ_ALIGNER|select:star,bowtie2|star|Read aligner|-m
1|aligner_args|Aligner Args|WIZ_ALIGNER_ARGS|str||Additional aligner pass-through args|
1|fastp_args|fastp Args|WIZ_FASTP_ARGS|str||Additional fastp pass-through args|
2|filter_repeat|Repeat Filter|WIZ_FILTER_REPEAT|yn|false|Enable repeat-element pre-filter|--filter-repeat
2|no_chr_filter|Disable Chr Filter|WIZ_NO_CHR_FILTER|yn|false|Skip canonical-chromosome filter|--no-chr-filter
3|peak_caller|Peak Caller|WIZ_PEAK_CALLER|select:homer,ctk|homer|Peak caller (HOMER or CTK)|--peak-caller
3|peak_dist|Peak Distance|WIZ_PEAK_DIST|int|50|Min distance between peaks (or -gap for CTK)|
3|peak_size|Peak Size|WIZ_PEAK_SIZE|int|20|Peak size (HOMER only; rolled into HOMER args)|
3|frag_len|Fragment Length|WIZ_FRAG_LEN|int|25|Fragment length (HOMER only; rolled into HOMER args)|
3|peak_caller_args|Peak Caller Args|WIZ_PEAK_CALLER_ARGS|str||Additional peak-caller args|--peak-caller-args
4|cims_iter|CIMS Iterations|WIZ_CIMS_ITER|int|5|CIMS permutation iterations|--cims-iter
4|cims_fdr|CIMS FDR|WIZ_CIMS_FDR|float|0.05|CIMS FDR threshold|--cims-fdr
4|cits_pval|CITS P-value|WIZ_CITS_PVAL|float|0.05|CITS p-value threshold|--cits-pval
4|cits_gap|CITS Gap|WIZ_CITS_GAP|int|25|CITS clustering gap (-1=no clustering)|--cits-gap
5|genome_fasta|Genome FASTA|WIZ_GENOME_FASTA|str||Path to genome FASTA (Clink uses)|--genome-fasta
5|clink_umi_len|Clink UMI Length|WIZ_CLINK_UMI_LEN|str||UMI length for umi_tools (blank=auto)|--clink-umi-len
5|clink_fdr|Clink FDR|WIZ_CLINK_FDR|float|0.05|Clink FDR threshold|--clink-fdr
5|clink_min_cov|Clink Min Coverage|WIZ_CLINK_MIN_COV|int|5|Minimum coverage to test|--clink-min-cov
5|clink_multi_map|Clink Multi-Map|WIZ_CLINK_MULTI_MAP|yn|false|EM-rescue of multi-mapped reads|--clink-multi-map
5|xl_bigwig|XL BigWig|WIZ_XL_BIGWIG|yn|false|Per-sample crosslink bigWigs|--xl-bigwig
6|no_motif|Skip Motif|WIZ_NO_MOTIF|yn|false|Skip flanked BED generation|--no-motif
6|flank|Flank|WIZ_FLANK|int|10|Flanked BED nucleotides|-f
7|groups_file|Groups File|WIZ_GROUPS_FILE|str||Path to groups file|-g
7|ctk_group|CTK Group|WIZ_CTK_GROUP|yn|false|Pool samples for CTK aggregation|--ctk-group
7|group_xlsite|Group Xlsite|WIZ_GROUP_XLSITE|yn|false|Group crosslink-site analysis|--group-xlsite
8|output|Output Name|WIZ_OUTPUT|str||Output folder name (blank=INPUT_output)|-o
8|keep|Keep Intermediate|WIZ_KEEP|yn|false|Keep intermediate files|-k
8|sample_size|Sample Size|WIZ_SAMPLE_SIZE|str||Test mode: process N reads only|-s
8|notification|Notifications|WIZ_NOTIFICATION|yn|false|System notifications on completion|--notification
EOF
}

_cc_category_name() {
    case "$1" in
        1) echo "Preprocessing" ;;
        2) echo "Filters" ;;
        3) echo "Peak Calling" ;;
        4) echo "CTK CIMS/CITS" ;;
        5) echo "Clink" ;;
        6) echo "Motif / Flanked BED" ;;
        7) echo "Grouping" ;;
        8) echo "Output / Runtime" ;;
    esac
}

# Return 0 if category is active for current WIZ_* state
_cc_should_show_category() {
    case "$1" in
        4) [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]] ;;
        5) [[ "$WIZ_RUN_CLINK" == "true" ]] ;;
        7) [[ "$WIZ_MODE" == "pooled" || "$WIZ_MODE" == "directory" ]] ;;
        *) return 0 ;;
    esac
}

# Return 0 if knob is active for current WIZ_* state
_cc_should_show_knob() {
    case "$1" in
        bc_len|spacer_len|bc_first|demux_mismatches)
            [[ "$WIZ_MODE" == "pooled" || "$WIZ_MODE" == "directory" ]] ;;
        peak_size|frag_len)
            [[ "$WIZ_PEAK_CALLER" == "homer" ]] ;;
        *) return 0 ;;
    esac
}

# Lookup full record for a knob key. Empty output if not found.
_cc_knob_record() {
    _cc_knob_registry | awk -F'|' -v k="$1" '$2==k {print; exit}'
}

# Field accessor — args: key, field_index (1-based)
_cc_field() {
    local rec; rec=$(_cc_knob_record "$1")
    [[ -z "$rec" ]] && return 1
    echo "$rec" | awk -F'|' -v f="$2" '{print $f}'
}

# Get current value of a knob (indirect expansion)
_cc_current() {
    local var; var=$(_cc_field "$1" 4)
    [[ -z "$var" ]] && return 1
    echo "${!var}"
}

# Display-format a knob value: yn→yes/no, blank→placeholder
_cc_display_value() {
    local key="$1"
    local type; type=$(_cc_field "$key" 5)
    local val; val=$(_cc_current "$key")
    case "$type" in
        yn) [[ "$val" == "true" ]] && echo "yes" || echo "no" ;;
        *) [[ -z "$val" ]] && echo "(blank)" || echo "$val" ;;
    esac
}

# True if current value matches default
_cc_is_default() {
    local def; def=$(_cc_field "$1" 6)
    local val; val=$(_cc_current "$1")
    [[ "$val" == "$def" ]]
}

# Print one-line help
_cc_print_help() {
    local label help
    label=$(_cc_field "$1" 3)
    help=$(_cc_field "$1" 7)
    echo -e "  ${WIZ_YELLOW}[$label]${WIZ_NC} $help"
}

# Render category submenu — args: category number
# Side effect: sets global _cc_pos_map = "key1 key2 key3..." (position → key lookup)
_cc_render_category() {
    local cat="$1"
    clear
    local name; name=$(_cc_category_name "$cat")
    echo -e "${WIZ_CYAN}─── ${WIZ_BOLD}Category $cat: $name${WIZ_NC}${WIZ_CYAN} ──────────────────────────${WIZ_NC}"
    echo ""

    _cc_pos_map=""
    local n=0
    while IFS='|' read -r c key label var type def help flag; do
        [[ "$c" != "$cat" ]] && continue
        _cc_should_show_knob "$key" || continue
        n=$((n+1))
        _cc_pos_map="$_cc_pos_map $key"
        local val; val=$(_cc_display_value "$key")
        local marker=" "
        _cc_is_default "$key" || marker="*"
        printf "  ${marker}[%2d] %-26s %s\n" "$n" "$label" "$val"
    done < <(_cc_knob_registry)

    echo ""
    echo "  [?]<n>  Show help for option n"
    echo "  [ b ]   Back to categories"
    echo "  [Enter] Accept and return"
    echo ""
    echo -e "  ${WIZ_YELLOW}* = changed from default${WIZ_NC}"
}

# Return knob key at position N in category (uses _cc_pos_map set by _cc_render_category)
_cc_knob_at_position() {
    local n="$1"
    local i=0
    for k in $_cc_pos_map; do
        i=$((i+1))
        [[ "$i" == "$n" ]] && { echo "$k"; return 0; }
    done
    return 1
}

# Spot-edit one knob (type-aware dispatch)
_cc_edit_knob() {
    local key="$1"
    local label var type def
    label=$(_cc_field "$key" 3)
    var=$(_cc_field "$key" 4)
    type=$(_cc_field "$key" 5)
    def=$(_cc_field "$key" 6)
    echo ""
    case "$type" in
        int)   prompt_value "  $label" "${!var:-$def}" "$var" "int" ;;
        float) prompt_value "  $label" "${!var:-$def}" "$var" "float" ;;
        yn)
            local cur_yn="n"; [[ "${!var}" == "true" ]] && cur_yn="y"
            local _new_yn
            prompt_yesno "  $label?" "$cur_yn" _new_yn
            if [[ "$_new_yn" == "y" ]]; then eval "$var=true"; else eval "$var=false"; fi
            ;;
        select:*)
            local opts_csv="${type#select:}"
            local opts_arr; IFS=',' read -r -a opts_arr <<< "$opts_csv"
            local _sel
            prompt_select "  Select $label:" _sel "${opts_arr[@]}"
            eval "$var=\"${opts_arr[$((_sel-1))]}\""
            ;;
        str|*)
            echo "  $label (current: ${!var:-blank}, blank input keeps current):"
            local _v
            read -p "  > " _v || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
            [[ -n "$_v" ]] && eval "$var=\"$_v\""
            ;;
    esac
}

# Category submenu loop
_cc_category_submenu() {
    local cat="$1"
    while true; do
        _cc_render_category "$cat"
        echo ""
        local choice
        read -p "> " choice || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
        case "$choice" in
            ""|b|B) return 0 ;;
            \?*)
                local num="${choice#\?}"
                local key; key=$(_cc_knob_at_position "$num") || { echo -e "${WIZ_RED}  Invalid number.${WIZ_NC}"; sleep 1; continue; }
                _cc_print_help "$key"
                read -p "  Press Enter to continue..." _ || exit 1
                ;;
            [0-9]*)
                local key; key=$(_cc_knob_at_position "$choice") || { echo -e "${WIZ_RED}  Invalid number.${WIZ_NC}"; sleep 1; continue; }
                _cc_edit_knob "$key"
                ;;
            *) echo -e "${WIZ_RED}  Invalid input. Use a number, ?<n>, b, or Enter.${WIZ_NC}"; sleep 1 ;;
        esac
    done
}

# Render top-level review menu (categories)
# Re-prompt the STEP 4 analysis-plan track selection. Current values pre-fill defaults.
_cc_select_tracks() {
    print_section "Analysis Plan — re-pick tracks"
    echo ""
    echo "  HOMER peaks runs by default. Pick additional tracks:"
    echo ""
    local _yn
    local _cims_def="n"; [[ "$WIZ_RUN_CIMS" == "true" ]] && _cims_def="y"
    local _cits_def="n"; [[ "$WIZ_RUN_CITS" == "true" ]] && _cits_def="y"
    local _clink_def="n"; [[ "$WIZ_RUN_CLINK" == "true" ]] && _clink_def="y"
    prompt_yesno "  Enable CTK CIMS analysis (mutation sites)?" "$_cims_def" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CIMS="true" || WIZ_RUN_CIMS="false"
    prompt_yesno "  Enable CTK CITS analysis (truncation sites)?" "$_cits_def" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CITS="true" || WIZ_RUN_CITS="false"
    prompt_yesno "  Enable Clink pipeline (Python pileup → CITS/CIMS)?" "$_clink_def" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CLINK="true" || WIZ_RUN_CLINK="false"
}

_cc_render_top_menu() {
    clear
    echo -e "${WIZ_CYAN}╔════════════════════════════════════════════════════════════════╗${WIZ_NC}"
    echo -e "${WIZ_CYAN}║${WIZ_NC}              ${WIZ_BOLD}STEP 5.5: Review Defaults${WIZ_NC}                       ${WIZ_CYAN}║${WIZ_NC}"
    echo -e "${WIZ_CYAN}╚════════════════════════════════════════════════════════════════╝${WIZ_NC}"
    echo ""
    local _tracks="HOMER"
    [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]] && _tracks="$_tracks + CTK"
    [[ "$WIZ_RUN_CLINK" == "true" ]] && _tracks="$_tracks + Clink"
    echo -e "  Selected analysis tracks: ${WIZ_GREEN}${_tracks}${WIZ_NC}"
    echo ""
    echo "  Which category would you like to review or edit?"
    echo ""
    local cat
    for cat in 1 2 3 4 5 6 7 8; do
        local name; name=$(_cc_category_name "$cat")
        if _cc_should_show_category "$cat"; then
            printf "  [%d] %s\n" "$cat" "$name"
        else
            printf "  ${WIZ_YELLOW}[%d] %s (not applicable)${WIZ_NC}\n" "$cat" "$name"
        fi
    done
    echo ""
    echo "  [ t ]   Change analysis plan (re-pick tracks)"
    echo "  [ a ]   Walk through ALL options step-by-step (advanced)"
    echo "  [Enter] Accept and continue to confirmation"
}

# Top-level review loop. Sets _cc_advanced_fallback=1 if user opts into advanced.
_cc_review_loop() {
    _cc_advanced_fallback=0
    while true; do
        _cc_render_top_menu
        echo ""
        local choice
        read -p "> " choice || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
        case "$choice" in
            "") return 0 ;;
            a|A)
                echo ""
                echo "  Switch to walking through ALL options step-by-step?"
                echo "  [y] Yes, switch to advanced"
                echo "  [n] No, cancel"
                echo "  [b] Back to review menu"
                local conf
                read -p "  > " conf || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
                case "$conf" in
                    y|Y) _cc_advanced_fallback=1; return 0 ;;
                    *) continue ;;
                esac
                ;;
            t|T) _cc_select_tracks ;;
            [1-8])
                _cc_should_show_category "$choice" || { echo -e "${WIZ_RED}  Category not active for your settings.${WIZ_NC}"; sleep 1; continue; }
                _cc_category_submenu "$choice"
                ;;
            *) echo -e "${WIZ_RED}  Invalid input. Use 1-8, t, a, or Enter.${WIZ_NC}"; sleep 1 ;;
        esac
    done
}

# Diff view: print only knobs changed from default
_cc_print_diff() {
    local any=0
    echo "  ┌─── Changed from defaults ──────────────────────────────────┐"
    while IFS='|' read -r c key label var type def help flag; do
        _cc_should_show_category "$c" || continue
        _cc_should_show_knob "$key" || continue
        _cc_is_default "$key" && continue
        any=1
        local val; val=$(_cc_display_value "$key")
        local def_disp="$def"; [[ -z "$def_disp" ]] && def_disp="(blank)"
        [[ "$type" == "yn" ]] && { [[ "$def" == "true" ]] && def_disp="yes" || def_disp="no"; }
        printf "  │ %-22s %-16s (default: %s)\n" "$label:" "$val" "$def_disp"
    done < <(_cc_knob_registry)
    [[ "$any" -eq 0 ]] && echo "  │   (all defaults — nothing changed)"
    echo "  └────────────────────────────────────────────────────────────┘"
}

# Build the equivalent CLIPittyClip.sh command line for the current WIZ_* state.
# Emits non-default flags only; omits knobs whose cli_flag is empty.
_cc_build_cmdline() {
    local cmd="CLIPittyClip.sh"
    # Input + index + barcodes/dir
    if [[ "$WIZ_MODE" == "pooled" ]]; then
        cmd="$cmd -i \"$WIZ_INPUT_FILE\" -b \"$WIZ_BARCODE_FILE\""
    elif [[ "$WIZ_MODE" == "single" ]]; then
        cmd="$cmd -i \"$WIZ_INPUT_FILE\""
    elif [[ "$WIZ_MODE" == "directory" ]]; then
        cmd="$cmd -d \"$WIZ_INPUT_DIR\""
    fi
    cmd="$cmd -x \"$WIZ_GENOME_INDEX\""
    [[ -n "$WIZ_ECLIP_MODE" ]]      && cmd="$cmd --eclip $WIZ_ECLIP_MODE"
    [[ "$WIZ_PARCLIP_MODE" == "true" ]] && cmd="$cmd --parclip"
    # Track flags
    if [[ "$WIZ_RUN_CIMS" == "true" && "$WIZ_RUN_CITS" == "true" ]]; then
        cmd="$cmd --run-cims-cits"
    else
        [[ "$WIZ_RUN_CIMS" == "true" ]] && cmd="$cmd --run-cims"
        [[ "$WIZ_RUN_CITS" == "true" ]] && cmd="$cmd --run-cits"
    fi
    [[ "$WIZ_RUN_CLINK" == "true" ]] && cmd="$cmd --run-clink"
    # Per-knob non-default flags
    while IFS='|' read -r c key label var type def help flag; do
        [[ -z "$flag" ]] && continue
        _cc_should_show_category "$c" || continue
        _cc_should_show_knob "$key" || continue
        _cc_is_default "$key" && continue
        local val; val="${!var}"
        if [[ "$type" == "yn" ]]; then
            [[ "$val" == "true" ]] && cmd="$cmd $flag"
        else
            [[ -n "$val" ]] && cmd="$cmd $flag \"$val\""
        fi
    done < <(_cc_knob_registry)
    echo "$cmd"
}

# ═══════════════════════════════════════════════════════════════════════════════
# CLIPittyClip WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_clipittyclip() {
    _play_wizard_intro
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
        echo ""
        prompt_file_path "  Enter path to pooled FASTQ file: " WIZ_INPUT_FILE || return 1
        prompt_file_path "  Enter path to barcode file: " WIZ_BARCODE_FILE || return 1
        WIZ_MODE="pooled"
    else
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

    # Optional output folder name (blank = INPUT_output)
    echo ""
    echo "  Output folder name (leave blank to use default: INPUT_output):"
    read -p "  > " WIZ_OUTPUT || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Protocol
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Protocol"
    echo ""
    prompt_select "Which CLIP protocol was used?" WIZ_PROTOCOL \
        "Standard CLIP (iCLIP, CoCLIP, etc.)" \
        "eCLIP — paired-end (post-eclipdemux R2, UMI in header)" \
        "eCLIP — single-end / seCLIP (raw R1, UMI in sequence)" \
        "PAR-CLIP (4SU; T>C signal)"

    WIZ_ECLIP_MODE=""
    WIZ_PARCLIP_MODE="false"
    case "$WIZ_PROTOCOL" in
        2) WIZ_ECLIP_MODE="pe" ;;
        3) WIZ_ECLIP_MODE="se" ;;
        4) WIZ_PARCLIP_MODE="true"
           echo ""
           echo -e "  ${WIZ_YELLOW}Note:${WIZ_NC} PAR-CLIP signal is T>C substitutions."
           echo -e "  Recommend enabling Clink in the next step." ;;
    esac

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Genome Index
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 3: Genome Index"
    echo ""
    prompt_dir_path "  Enter path to genome index directory: " WIZ_GENOME_INDEX || return 1
    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4: Analysis Plan (which tracks to run)
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 4: Analysis Plan"
    WIZ_RUN_CIMS="false"
    WIZ_RUN_CITS="false"
    WIZ_RUN_CLINK="false"
    local _yn
    _cc_select_tracks

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Settings tier
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 5: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings (you can review/spot-edit them next)" \
        "Modify program settings step-by-step (advanced)"

    # Defaults (must match main-branch CLIPittyClip.sh AND knob registry)
    WIZ_ALIGNER="star"
    WIZ_THREADS="1"
    WIZ_UMI_LEN="0"
    WIZ_ADAPTER="L32"
    WIZ_BC_LEN=""
    WIZ_SPACER_LEN="0"
    WIZ_BC_FIRST="false"
    WIZ_MIN_QUAL="30"
    WIZ_NO_DEDUP="false"
    WIZ_ALIGN_MISMATCHES="2"
    WIZ_DEMUX_MISMATCHES="1"
    WIZ_FILTER_REPEAT="false"
    WIZ_NO_CHR_FILTER="false"
    WIZ_GENOME_FASTA=""
    WIZ_PEAK_CALLER="homer"
    WIZ_PEAK_CALLER_ARGS=""
    WIZ_PEAK_DIST="50"
    WIZ_PEAK_SIZE="20"
    WIZ_FRAG_LEN="25"
    WIZ_CIMS_ITER="5"
    WIZ_CIMS_FDR="0.05"
    WIZ_CITS_PVAL="0.05"
    WIZ_CITS_GAP="25"
    WIZ_CLINK_UMI_LEN=""
    WIZ_CLINK_FDR="0.05"
    WIZ_CLINK_MIN_COV="5"
    WIZ_CLINK_MULTI_MAP="false"
    WIZ_NO_MOTIF="false"
    WIZ_FLANK="10"
    WIZ_GROUPS_FILE=""
    WIZ_CTK_GROUP="false"
    WIZ_GROUP_XLSITE="false"
    WIZ_KEEP="false"
    WIZ_SAMPLE_SIZE=""
    WIZ_NOTIFICATION="false"
    WIZ_XL_BIGWIG="false"
    WIZ_ALIGNER_ARGS=""
    WIZ_FASTP_ARGS=""

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5.5: Categorized default review (only when defaults picked)
    # Sets _cc_advanced_fallback=1 if user escapes to advanced via [a].
    # ─────────────────────────────────────────────────────────────────────────
    if [[ "$WIZ_SETTINGS_MODE" == "1" ]]; then
        _cc_review_loop
        [[ "$_cc_advanced_fallback" == "1" ]] && WIZ_SETTINGS_MODE="2"
    fi

    # ─────────────────────────────────────────────────────────────────────────
    # Advanced walkthrough (picked at STEP 5, or via [a] escape from review)
    # ─────────────────────────────────────────────────────────────────────────
    if [[ "$WIZ_SETTINGS_MODE" == "2" ]]; then
        print_doc_box "CLIPittyClip" \
            "CLIPittyClip.sh --help" \
            "STAR: https://github.com/alexdobin/STAR" \
            "Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2" \
            "CTK: https://zhanglab.c2b2.columbia.edu/" \
            "HOMER: http://homer.ucsd.edu/homer/"

        # Genome FASTA (only relevant if Clink picked)
        if [[ "$WIZ_RUN_CLINK" == "true" ]]; then
            echo ""
            echo "  Clink track can use a genome FASTA. Leave blank to skip."
            read -p "  Genome FASTA path (optional): " WIZ_GENOME_FASTA || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
            WIZ_GENOME_FASTA="${WIZ_GENOME_FASTA/#\~/$HOME}"
        fi

        # ─── PREPROCESSING ─────────────────────────────────────────────────
        print_section "PREPROCESSING SETTINGS"
        echo ""
        prompt_value "  Enter number of threads" "1" WIZ_THREADS "int"
        prompt_value "  Enter UMI length (0 if none)" "0" WIZ_UMI_LEN "int"
        echo "  Adapter options: L32 (default), L19, or custom sequence"
        prompt_value "  Enter 3' adapter" "L32" WIZ_ADAPTER

        if [[ "$WIZ_MODE" == "pooled" || "$WIZ_MODE" == "directory" ]]; then
            echo ""
            echo "  Barcode/spacer layout (used when demultiplexing):"
            prompt_value "  Barcode length (blank = auto-detect from barcode file)" "" WIZ_BC_LEN
            prompt_value "  Spacer length to trim after barcode" "0" WIZ_SPACER_LEN "int"
            prompt_yesno "  Barcode precedes UMI? (e.g., BrdU-CLIP, irCLIP2)" "n" _yn
            [[ "$_yn" == "y" ]] && WIZ_BC_FIRST="true"
            prompt_value "  Max barcode mismatches during demux" "1" WIZ_DEMUX_MISMATCHES "int"
        fi

        prompt_value "  fastp average quality threshold" "30" WIZ_MIN_QUAL "int"
        prompt_yesno "  Disable PCR deduplication?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_NO_DEDUP="true"
        prompt_value "  Max alignment mismatches (STAR)" "2" WIZ_ALIGN_MISMATCHES "int"

        echo ""
        local _sel
        prompt_select "Select aligner:" _sel "STAR (default)" "Bowtie2"
        [[ "$_sel" == "2" ]] && WIZ_ALIGNER="bowtie2"
        echo "  Enter additional aligner arguments (optional, blank to skip):"
        read -p "  > " WIZ_ALIGNER_ARGS || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

        echo "  Enter additional fastp arguments (optional, blank to skip):"
        read -p "  > " WIZ_FASTP_ARGS || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

        # ─── FILTERS ───────────────────────────────────────────────────────
        print_section "FILTERS"
        echo ""
        prompt_yesno "  Enable repeat-element pre-filter? (default: off)" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_FILTER_REPEAT="true"
        prompt_yesno "  Disable canonical-chromosome filter? (keep on by default)" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_NO_CHR_FILTER="true"

        # ─── PEAK CALLING ──────────────────────────────────────────────────
        print_section "PEAK CALLING"
        echo ""
        prompt_select "Select peak caller:" _sel \
            "HOMER findPeaks (default)" \
            "CTK tag2peak.pl"
        if [[ "$_sel" == "2" ]]; then
            WIZ_PEAK_CALLER="ctk"
            echo -e "  ${WIZ_YELLOW}CTK tag2peak defaults:${WIZ_NC} -big -ss --valley-seeking"
            prompt_value "  Enter gap (-gap)" "50" WIZ_PEAK_DIST "int"
        else
            echo -e "  ${WIZ_YELLOW}HOMER defaults:${WIZ_NC} -style factor · -L 2 · -localSize 1000 · -minDist 50 · -size 20 · -fragLength 25"
            prompt_value "  Enter min distance between peaks" "50" WIZ_PEAK_DIST "int"
            prompt_value "  Enter peak size" "20" WIZ_PEAK_SIZE "int"
            prompt_value "  Enter fragment length" "25" WIZ_FRAG_LEN "int"
        fi
        echo "  Enter additional peak-caller arguments (optional, blank to skip):"
        read -p "  > " WIZ_PEAK_CALLER_ARGS || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

        # ─── CTK CIMS/CITS (only if track selected at STEP 4) ──────────────
        if [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]]; then
            print_section "CTK CIMS/CITS"
            echo -e "  ${WIZ_YELLOW}Defaults:${WIZ_NC} cims-iter=5, cims-fdr=0.05, cits-pval=0.05, cits-gap=25"
            prompt_value "  CIMS permutation iterations" "5" WIZ_CIMS_ITER "int"
            prompt_value "  CIMS FDR threshold" "0.05" WIZ_CIMS_FDR "float"
            prompt_value "  CITS p-value threshold" "0.05" WIZ_CITS_PVAL "float"
            prompt_value "  CITS clustering gap (-1 = no clustering)" "25" WIZ_CITS_GAP "int"
        fi

        # ─── Clink (only if track selected at STEP 4) ──────────────────────
        if [[ "$WIZ_RUN_CLINK" == "true" ]]; then
            print_section "CLINK PIPELINE"
            echo -e "  ${WIZ_YELLOW}Defaults:${WIZ_NC} fdr=0.05, min-cov=5"
            prompt_value "  Clink UMI length (blank = auto-detect)" "" WIZ_CLINK_UMI_LEN
            prompt_value "  Clink FDR threshold" "0.05" WIZ_CLINK_FDR "float"
            prompt_value "  Clink minimum coverage" "5" WIZ_CLINK_MIN_COV "int"
            prompt_yesno "  Rescue multi-mapped reads via EM (--clink-multi-map)?" "n" _yn
            [[ "$_yn" == "y" ]] && WIZ_CLINK_MULTI_MAP="true"
        fi

        # ─── Motif / flanked BED ───────────────────────────────────────────
        print_section "MOTIF / FLANKED BED"
        echo ""
        prompt_yesno "  Skip flanked BED generation (--no-motif)?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_NO_MOTIF="true"
        if [[ "$WIZ_NO_MOTIF" == "false" ]]; then
            prompt_value "  Flanked BED nucleotides" "10" WIZ_FLANK "int"
        fi

        # ─── GROUPING (conditional on multi-sample input) ──────────────────
        if [[ "$WIZ_MODE" == "directory" || "$WIZ_MODE" == "pooled" ]]; then
            print_section "GROUP ANALYSIS"
            echo ""
            prompt_yesno "  Provide a groups file for aggregation?" "n" _yn
            if [[ "$_yn" == "y" ]]; then
                prompt_file_path "  Enter path to groups file: " WIZ_GROUPS_FILE || return 1
                prompt_yesno "  Enable --ctk-group (pool samples for CTK)?" "n" _yn
                [[ "$_yn" == "y" ]] && WIZ_CTK_GROUP="true"
                prompt_yesno "  Enable --group-xlsite (group crosslink-site analysis)?" "n" _yn
                [[ "$_yn" == "y" ]] && WIZ_GROUP_XLSITE="true"
            fi
        fi

        # ─── OUTPUT / RUNTIME ──────────────────────────────────────────────
        # (Output folder name was collected in STEP 1; review/edit via category 8.)
        print_section "OUTPUT / RUNTIME"
        echo ""
        prompt_yesno "  Keep intermediate files (-k)?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_KEEP="true"
        prompt_value "  Test mode: process N reads only (blank = all)" "" WIZ_SAMPLE_SIZE
        prompt_yesno "  Enable system notifications on completion?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_NOTIFICATION="true"

        if [[ "$WIZ_RUN_CLINK" == "true" ]]; then
            prompt_yesno "  Generate per-sample crosslink bigWigs (--xl-bigwig)?" "n" _yn
            [[ "$_yn" == "y" ]] && WIZ_XL_BIGWIG="true"
        fi
    fi

    # ─────────────────────────────────────────────────────────────────────────
    # FINAL CONFIRMATION: summary + diff view + equivalent command line + confirm
    # ─────────────────────────────────────────────────────────────────────────
    clear
    echo -e "${WIZ_CYAN}╔════════════════════════════════════════════════════════════════╗${WIZ_NC}"
    echo -e "${WIZ_CYAN}║${WIZ_NC}                  ${WIZ_BOLD}FINAL CONFIRMATION${WIZ_NC}                          ${WIZ_CYAN}║${WIZ_NC}"
    echo -e "${WIZ_CYAN}║${WIZ_NC}        ${WIZ_YELLOW}Review the settings below before starting.${WIZ_NC}              ${WIZ_CYAN}║${WIZ_NC}"
    echo -e "${WIZ_CYAN}╚════════════════════════════════════════════════════════════════╝${WIZ_NC}"
    echo ""
    # Input/protocol/genome — always shown
    echo "  ┌─── Inputs ─────────────────────────────────────────────────┐"
    if [[ "$WIZ_MODE" == "pooled" ]]; then
        printf "  │ %-18s %-40s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
        printf "  │ %-18s %-40s │\n" "Barcodes:" "$(basename "$WIZ_BARCODE_FILE")"
    elif [[ "$WIZ_MODE" == "single" ]]; then
        printf "  │ %-18s %-40s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
    else
        printf "  │ %-18s %-40s │\n" "Input Dir:" "$WIZ_INPUT_DIR"
    fi
    local _proto="Standard CLIP"
    [[ "$WIZ_ECLIP_MODE" == "pe" ]] && _proto="eCLIP-pe"
    [[ "$WIZ_ECLIP_MODE" == "se" ]] && _proto="eCLIP-se"
    [[ "$WIZ_PARCLIP_MODE" == "true" ]] && _proto="PAR-CLIP"
    printf "  │ %-18s %-40s │\n" "Protocol:" "$_proto"
    printf "  │ %-18s %-40s │\n" "Genome Index:" "$(basename "$WIZ_GENOME_INDEX")"
    local _tracks="HOMER"
    [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]] && _tracks="$_tracks + CTK"
    [[ "$WIZ_RUN_CLINK" == "true" ]] && _tracks="$_tracks + Clink"
    printf "  │ %-18s %-40s │\n" "Tracks:" "$_tracks"
    echo "  └────────────────────────────────────────────────────────────┘"
    echo ""

    # Diff view: only knobs changed from default
    _cc_print_diff
    echo ""

    # Equivalent CLI command — copy/paste to skip wizard next time
    local _cmd; _cmd=$(_cc_build_cmdline)
    CLIPITTY_EQUIV_CMD="$_cmd"
    echo -e "  ${WIZ_YELLOW}Equivalent command (copy to skip the wizard next time):${WIZ_NC}"
    echo "    $_cmd"
    echo ""

    prompt_yesno "  Start analysis with these settings?" "y" WIZ_CONFIRM
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi

    # Export all WIZ_* vars for host script
    export WIZ_MODE WIZ_INPUT_FILE WIZ_INPUT_DIR WIZ_BARCODE_FILE WIZ_GENOME_INDEX WIZ_GENOME_FASTA
    export WIZ_ECLIP_MODE WIZ_PARCLIP_MODE
    export WIZ_ALIGNER WIZ_THREADS WIZ_UMI_LEN WIZ_ADAPTER
    export WIZ_BC_LEN WIZ_SPACER_LEN WIZ_BC_FIRST WIZ_MIN_QUAL
    export WIZ_NO_DEDUP WIZ_ALIGN_MISMATCHES WIZ_DEMUX_MISMATCHES
    export WIZ_FILTER_REPEAT WIZ_NO_CHR_FILTER
    export WIZ_ALIGNER_ARGS WIZ_FASTP_ARGS WIZ_PEAK_CALLER_ARGS
    export WIZ_PEAK_CALLER WIZ_PEAK_DIST WIZ_PEAK_SIZE WIZ_FRAG_LEN
    export WIZ_RUN_CIMS WIZ_RUN_CITS WIZ_CIMS_ITER WIZ_CIMS_FDR WIZ_CITS_PVAL WIZ_CITS_GAP
    export WIZ_RUN_CLINK WIZ_CLINK_UMI_LEN WIZ_CLINK_FDR WIZ_CLINK_MIN_COV WIZ_CLINK_MULTI_MAP
    export WIZ_NO_MOTIF WIZ_FLANK
    export WIZ_GROUPS_FILE WIZ_CTK_GROUP WIZ_GROUP_XLSITE
    export WIZ_OUTPUT WIZ_KEEP WIZ_SAMPLE_SIZE WIZ_NOTIFICATION WIZ_XL_BIGWIG
    export CLIPITTY_EQUIV_CMD

    echo -e "${WIZ_GREEN}Starting analysis...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# PREPittyPrep WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_prepittyprep() {
    _play_wizard_intro
    print_wizard_header "PREPittyPrep Wizard"
    local _yn

    # ─── STEP 1: Input ─────────────────────────────────────────────────────
    print_section "STEP 1: Input"
    echo ""
    prompt_select "What is your input?" WIZ_INPUT_TYPE \
        "Single FASTQ file" \
        "Directory of FASTQ files (batch mode)"
    if [[ "$WIZ_INPUT_TYPE" == "1" ]]; then
        prompt_file_path "  Enter path to FASTQ file: " WIZ_INPUT_FILE || return 1
        WIZ_MODE="single"
    else
        prompt_dir_path "  Enter directory containing FASTQ files: " WIZ_INPUT_DIR || return 1
        WIZ_MODE="directory"
    fi

    WIZ_BARCODE_FILE=""
    if [[ "$WIZ_MODE" == "single" ]]; then
        echo ""
        prompt_yesno "  Is this a pooled library (need to demultiplex by barcode)?" "n" _yn
        if [[ "$_yn" == "y" ]]; then
            prompt_file_path "  Enter barcode file: " WIZ_BARCODE_FILE || return 1
        fi
    fi

    # ─── STEP 2: Mode (standard vs GEO) ────────────────────────────────────
    print_section "STEP 2: Mode"
    echo ""
    echo "  Standard prep: dedup → demux → fastp → *_prepped.fastq.gz (ready to map)"
    echo "  GEO deposit:   raw barcode-split FASTQ only (no dedup, no fastp)"
    echo ""
    prompt_select "Which mode?" WIZ_MODE_SEL \
        "Standard prep" \
        "GEO deposit (requires barcode file + single -i)"

    WIZ_GEO_MODE="false"
    if [[ "$WIZ_MODE_SEL" == "2" ]]; then
        if [[ -z "$WIZ_BARCODE_FILE" || "$WIZ_MODE" != "single" ]]; then
            echo -e "${WIZ_RED}  ✗ GEO mode requires a single -i input AND a barcode file.${WIZ_NC}" >&2
            return 1
        fi
        WIZ_GEO_MODE="true"
    fi

    # Defaults (match PREPittyPrep.sh)
    WIZ_THREADS="1"
    WIZ_UMI_LEN="0"
    WIZ_ADAPTER="AGATCGGAAGAGC"
    WIZ_BC_LEN=""
    WIZ_SPACER_LEN="0"
    WIZ_NO_DEDUP="false"
    WIZ_ECLIP="false"
    WIZ_DEMUX_MISMATCHES="1"
    WIZ_FILTER_NCRNA="false"
    WIZ_GENOME_INDEX=""
    WIZ_OUTPUT=""
    WIZ_KEEP="false"
    WIZ_SAMPLE_SIZE="0"

    if [[ "$WIZ_GEO_MODE" == "true" ]]; then
        echo ""
        echo -e "  ${WIZ_YELLOW}GEO mode: skipping preprocessing & ncRNA prompts.${WIZ_NC}"
    else
        # ─── STEP 3: Preprocessing ─────────────────────────────────────────
        print_section "STEP 3: Preprocessing"
        echo ""
        prompt_value "  Enter UMI length (0 if none)" "0" WIZ_UMI_LEN "int"
        echo "  Adapter default: AGATCGGAAGAGC (Illumina universal)"
        prompt_value "  Enter 3' adapter" "AGATCGGAAGAGC" WIZ_ADAPTER

        if [[ -n "$WIZ_BARCODE_FILE" ]]; then
            echo ""
            echo "  Barcode/spacer layout (used during demultiplexing):"
            prompt_value "  Barcode length (blank = auto-detect)" "" WIZ_BC_LEN
            prompt_value "  Spacer length to trim after barcode" "0" WIZ_SPACER_LEN "int"
            prompt_value "  Max barcode mismatches during demux" "1" WIZ_DEMUX_MISMATCHES "int"
        fi

        echo ""
        prompt_yesno "  Disable PCR deduplication?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_NO_DEDUP="true"
        prompt_yesno "  eCLIP mode (UMI in read header)?" "n" _yn
        [[ "$_yn" == "y" ]] && WIZ_ECLIP="true"

        # ─── STEP 4: ncRNA filter (optional) ───────────────────────────────
        print_section "STEP 4: ncRNA Filter"
        echo ""
        prompt_yesno "  Enable ncRNA pre-filtering? (requires genome index)" "n" _yn
        if [[ "$_yn" == "y" ]]; then
            WIZ_FILTER_NCRNA="true"
            prompt_dir_path "  Enter path to genome index: " WIZ_GENOME_INDEX || return 1
        fi
    fi

    # ─── STEP 5: Output / runtime ──────────────────────────────────────────
    print_section "STEP 5: Output / Runtime"
    echo ""
    prompt_value "  Output folder name (blank = INPUT_output)" "" WIZ_OUTPUT
    prompt_value "  Number of threads" "1" WIZ_THREADS "int"
    prompt_yesno "  Keep intermediate files (-k)?" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_KEEP="true"
    prompt_value "  Test mode: process N reads only (0 = all)" "0" WIZ_SAMPLE_SIZE "int"

    # ─── Summary + equivalent command ──────────────────────────────────────
    print_summary_box
    echo ""
    echo "  ┌─── PREPittyPrep configuration ─────────────────────────────┐"
    if [[ "$WIZ_MODE" == "single" ]]; then
        printf "  │ %-18s %-40s │\n" "Input:" "$(basename "$WIZ_INPUT_FILE")"
    else
        printf "  │ %-18s %-40s │\n" "Input Dir:" "$WIZ_INPUT_DIR"
    fi
    [[ -n "$WIZ_BARCODE_FILE" ]] && printf "  │ %-18s %-40s │\n" "Barcodes:" "$(basename "$WIZ_BARCODE_FILE")"
    local _m="Standard prep"; [[ "$WIZ_GEO_MODE" == "true" ]] && _m="GEO deposit"
    printf "  │ %-18s %-40s │\n" "Mode:" "$_m"
    if [[ "$WIZ_GEO_MODE" != "true" ]]; then
        printf "  │ %-18s %-40s │\n" "UMI Length:" "$WIZ_UMI_LEN"
        printf "  │ %-18s %-40s │\n" "Adapter:" "$WIZ_ADAPTER"
        local _d="Enabled"; [[ "$WIZ_NO_DEDUP" == "true" ]] && _d="Disabled"
        printf "  │ %-18s %-40s │\n" "Dedup:" "$_d"
        [[ "$WIZ_ECLIP" == "true" ]] && printf "  │ %-18s %-40s │\n" "eCLIP mode:" "yes"
        local _f="Disabled"; [[ "$WIZ_FILTER_NCRNA" == "true" ]] && _f="Enabled"
        printf "  │ %-18s %-40s │\n" "ncRNA Filter:" "$_f"
    fi
    printf "  │ %-18s %-40s │\n" "Threads:" "$WIZ_THREADS"
    [[ "$WIZ_KEEP" == "true" ]] && printf "  │ %-18s %-40s │\n" "Keep:" "yes"
    echo "  └────────────────────────────────────────────────────────────┘"
    echo ""

    # Build equivalent CLI
    local _cmd="PREPittyPrep.sh"
    if [[ "$WIZ_MODE" == "single" ]]; then
        _cmd="$_cmd -i \"$WIZ_INPUT_FILE\""
    else
        _cmd="$_cmd -d \"$WIZ_INPUT_DIR\""
    fi
    [[ -n "$WIZ_BARCODE_FILE" ]] && _cmd="$_cmd -b \"$WIZ_BARCODE_FILE\""
    [[ "$WIZ_GEO_MODE" == "true" ]] && _cmd="$_cmd --geo"
    [[ "$WIZ_UMI_LEN" != "0" ]] && _cmd="$_cmd -u $WIZ_UMI_LEN"
    [[ "$WIZ_ADAPTER" != "AGATCGGAAGAGC" ]] && _cmd="$_cmd -a \"$WIZ_ADAPTER\""
    [[ -n "$WIZ_BC_LEN" ]] && _cmd="$_cmd --bc-len $WIZ_BC_LEN"
    [[ "$WIZ_SPACER_LEN" != "0" ]] && _cmd="$_cmd --spacer-len $WIZ_SPACER_LEN"
    [[ "$WIZ_NO_DEDUP" == "true" ]] && _cmd="$_cmd --no-dedup"
    [[ "$WIZ_ECLIP" == "true" ]] && _cmd="$_cmd --eclip"
    [[ "$WIZ_DEMUX_MISMATCHES" != "1" ]] && _cmd="$_cmd --demux-mismatches $WIZ_DEMUX_MISMATCHES"
    [[ "$WIZ_FILTER_NCRNA" == "true" ]] && _cmd="$_cmd --filter-ncrna -x \"$WIZ_GENOME_INDEX\""
    [[ -n "$WIZ_OUTPUT" ]] && _cmd="$_cmd -o \"$WIZ_OUTPUT\""
    [[ "$WIZ_KEEP" == "true" ]] && _cmd="$_cmd -k"
    [[ "$WIZ_THREADS" != "1" ]] && _cmd="$_cmd -t $WIZ_THREADS"
    [[ "$WIZ_SAMPLE_SIZE" != "0" ]] && _cmd="$_cmd -s $WIZ_SAMPLE_SIZE"
    CLIPITTY_EQUIV_CMD="$_cmd"

    echo -e "  ${WIZ_YELLOW}Equivalent command (copy to skip the wizard next time):${WIZ_NC}"
    echo "    $_cmd"
    echo ""

    prompt_yesno "  Start preprocessing with these settings?" "y" WIZ_CONFIRM
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi

    export WIZ_MODE WIZ_INPUT_FILE WIZ_INPUT_DIR WIZ_BARCODE_FILE WIZ_GEO_MODE
    export WIZ_THREADS WIZ_UMI_LEN WIZ_ADAPTER WIZ_BC_LEN WIZ_SPACER_LEN
    export WIZ_NO_DEDUP WIZ_ECLIP WIZ_DEMUX_MISMATCHES
    export WIZ_FILTER_NCRNA WIZ_GENOME_INDEX
    export WIZ_OUTPUT WIZ_KEEP WIZ_SAMPLE_SIZE
    export CLIPITTY_EQUIV_CMD

    echo -e "${WIZ_GREEN}Starting preprocessing...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAPittyMap WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_mapittymap() {
    _play_wizard_intro
    print_wizard_header "MAPittyMap Wizard"
    local _yn

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Input File
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 1: Input File"
    echo ""
    prompt_file_path "  Enter path to input FASTQ file: " WIZ_INPUT_FILE || return 1

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Genome Index (+ optional --genome-fasta)
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Genome Index"
    echo ""
    prompt_dir_path "  Enter path to genome index directory: " WIZ_GENOME_INDEX || return 1

    WIZ_GENOME_FASTA=""
    echo ""
    echo "  Optional: genome FASTA path (recommended for CIMS; leave blank to skip)."
    read -p "  > " WIZ_GENOME_FASTA || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
    WIZ_GENOME_FASTA="${WIZ_GENOME_FASTA/#\~/$HOME}"

    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Filters
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 3: Filters"
    echo ""
    WIZ_FILTER_REPEAT="false"
    WIZ_NO_CHR_FILTER="false"
    prompt_yesno "  Enable repeat-element pre-filter? (default: off)" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_FILTER_REPEAT="true"
    prompt_yesno "  Disable canonical-chromosome filter? (kept ON by default)" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_NO_CHR_FILTER="true"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4: Analysis Plan (crosslink-site tracks)
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 4: Analysis Plan"
    echo ""
    echo "  Mapping always produces BAM + collapsed BED + bedgraph."
    echo "  Optional crosslink-site tracks (each is independent):"
    echo ""
    WIZ_RUN_CIMS="false"
    WIZ_RUN_CITS="false"
    WIZ_RUN_CLINK="false"
    prompt_yesno "  Enable CTK CIMS analysis (mutation sites)?" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CIMS="true"
    prompt_yesno "  Enable CTK CITS analysis (truncation sites)?" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CITS="true"
    prompt_yesno "  Enable Clink pipeline (Python pileup → CITS/CIMS)?" "n" _yn
    [[ "$_yn" == "y" ]] && WIZ_RUN_CLINK="true"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Default or Advanced
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 5: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings" \
        "Modify program settings (advanced)"

    # Initialize defaults
    WIZ_ALIGNER="star"
    WIZ_THREADS="1"
    WIZ_ALIGN_MISMATCHES="2"
    WIZ_OUTPUT_NAME=""
    WIZ_ALIGNER_ARGS=""
    WIZ_UMI_LEN="0"
    WIZ_CIMS_ITER="5"
    WIZ_CIMS_FDR="0.05"
    WIZ_CITS_PVAL="0.05"
    WIZ_CITS_GAP="25"
    WIZ_CLINK_UMI_LEN=""
    WIZ_CLINK_FDR="0.05"
    WIZ_CLINK_MIN_COV="5"
    WIZ_CLINK_MULTI_MAP="false"
    WIZ_NO_MOTIF="false"
    WIZ_FLANK="10"

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
        prompt_value "  Enter UMI length (0 if none)" "0" WIZ_UMI_LEN "int"

        local default_name=$(basename "$WIZ_INPUT_FILE" .fastq.gz)
        prompt_value "  Enter output name" "$default_name" WIZ_OUTPUT_NAME

        echo ""
        echo "  Enter additional aligner arguments (optional, blank to skip):"
        read -p "  > " WIZ_ALIGNER_ARGS || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

        # Track-conditional knob prompts
        if [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]]; then
            print_section "CTK CIMS/CITS"
            echo -e "  ${WIZ_YELLOW}Defaults:${WIZ_NC} cims-iter=5, cims-fdr=0.05, cits-pval=0.05, cits-gap=25"
            prompt_value "  CIMS permutation iterations" "5" WIZ_CIMS_ITER "int"
            prompt_value "  CIMS FDR threshold" "0.05" WIZ_CIMS_FDR "float"
            prompt_value "  CITS p-value threshold" "0.05" WIZ_CITS_PVAL "float"
            prompt_value "  CITS clustering gap (-1 = no clustering)" "25" WIZ_CITS_GAP "int"
        fi

        if [[ "$WIZ_RUN_CLINK" == "true" ]]; then
            print_section "CLINK PIPELINE"
            echo -e "  ${WIZ_YELLOW}Defaults:${WIZ_NC} fdr=0.05, min-cov=5"
            prompt_value "  Clink UMI length (blank = auto-detect)" "" WIZ_CLINK_UMI_LEN
            prompt_value "  Clink FDR threshold" "0.05" WIZ_CLINK_FDR "float"
            prompt_value "  Clink minimum coverage" "5" WIZ_CLINK_MIN_COV "int"
            prompt_yesno "  Rescue multi-mapped reads via EM (--clink-multi-map)?" "n" _yn
            [[ "$_yn" == "y" ]] && WIZ_CLINK_MULTI_MAP="true"
        fi

        if [[ "$WIZ_RUN_CIMS" == "true" || "$WIZ_RUN_CITS" == "true" ]]; then
            print_section "MOTIF / FLANKED BED"
            echo ""
            prompt_yesno "  Skip flanked BED generation (--no-motif)?" "n" _yn
            [[ "$_yn" == "y" ]] && WIZ_NO_MOTIF="true"
            if [[ "$WIZ_NO_MOTIF" == "false" ]]; then
                prompt_value "  Flanked BED nucleotides" "10" WIZ_FLANK "int"
            fi
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
    [[ -n "$WIZ_GENOME_FASTA" ]] && printf "  │ %-20s %-38s │\n" "Genome FASTA:" "$(basename "$WIZ_GENOME_FASTA")"
    local _rf="Disabled"; [[ "$WIZ_FILTER_REPEAT" == "true" ]] && _rf="Enabled"
    local _cf="Enabled";  [[ "$WIZ_NO_CHR_FILTER" == "true" ]] && _cf="Disabled"
    printf "  │ %-20s %-38s │\n" "Repeat Filter:" "$_rf"
    printf "  │ %-20s %-38s │\n" "Chr Filter:" "$_cf"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Aligner:" "$(echo "$WIZ_ALIGNER" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Threads:" "$WIZ_THREADS"
    printf "  │ %-20s %-38s │\n" "Max Mismatches:" "$WIZ_ALIGN_MISMATCHES"
    [[ "$WIZ_UMI_LEN" != "0" ]] && printf "  │ %-20s %-38s │\n" "UMI Length:" "$WIZ_UMI_LEN"
    [[ -n "$WIZ_OUTPUT_NAME" ]] && printf "  │ %-20s %-38s │\n" "Output Name:" "$WIZ_OUTPUT_NAME"
    [[ -n "$WIZ_ALIGNER_ARGS" ]] && printf "  │ %-20s %-38s │\n" "Aligner Args:" "$WIZ_ALIGNER_ARGS"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    local _cims="Disabled"; [[ "$WIZ_RUN_CIMS" == "true" ]] && _cims="Enabled"
    local _cits="Disabled"; [[ "$WIZ_RUN_CITS" == "true" ]] && _cits="Enabled"
    local _clink="Disabled"; [[ "$WIZ_RUN_CLINK" == "true" ]] && _clink="Enabled"
    printf "  │ %-20s %-38s │\n" "CTK CIMS:" "$_cims"
    printf "  │ %-20s %-38s │\n" "CTK CITS:" "$_cits"
    printf "  │ %-20s %-38s │\n" "Clink:" "$_clink"
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""

    # Equivalent command line
    local _cmd="MAPittyMap.sh -i \"$WIZ_INPUT_FILE\" -x \"$WIZ_GENOME_INDEX\""
    [[ -n "$WIZ_GENOME_FASTA" ]] && _cmd="$_cmd --genome-fasta \"$WIZ_GENOME_FASTA\""
    [[ "$WIZ_ALIGNER" != "star" ]] && _cmd="$_cmd --aligner $WIZ_ALIGNER"
    [[ "$WIZ_THREADS" != "1" ]] && _cmd="$_cmd -t $WIZ_THREADS"
    [[ "$WIZ_ALIGN_MISMATCHES" != "2" ]] && _cmd="$_cmd -m $WIZ_ALIGN_MISMATCHES"
    [[ "$WIZ_UMI_LEN" != "0" ]] && _cmd="$_cmd -u $WIZ_UMI_LEN"
    [[ -n "$WIZ_OUTPUT_NAME" ]] && _cmd="$_cmd -o \"$WIZ_OUTPUT_NAME\""
    [[ "$WIZ_FILTER_REPEAT" == "true" ]] && _cmd="$_cmd --filter-repeat"
    [[ "$WIZ_NO_CHR_FILTER" == "true" ]] && _cmd="$_cmd --no-chr-filter"
    if [[ "$WIZ_RUN_CIMS" == "true" && "$WIZ_RUN_CITS" == "true" ]]; then
        _cmd="$_cmd --run-cims-cits"
    else
        [[ "$WIZ_RUN_CIMS" == "true" ]] && _cmd="$_cmd --run-cims"
        [[ "$WIZ_RUN_CITS" == "true" ]] && _cmd="$_cmd --run-cits"
    fi
    [[ "$WIZ_CIMS_ITER" != "5" ]] && _cmd="$_cmd --cims-iter $WIZ_CIMS_ITER"
    [[ "$WIZ_CIMS_FDR" != "0.05" ]] && _cmd="$_cmd --cims-fdr $WIZ_CIMS_FDR"
    [[ "$WIZ_CITS_PVAL" != "0.05" ]] && _cmd="$_cmd --cits-pval $WIZ_CITS_PVAL"
    [[ "$WIZ_CITS_GAP" != "25" ]] && _cmd="$_cmd --cits-gap $WIZ_CITS_GAP"
    [[ "$WIZ_RUN_CLINK" == "true" ]] && _cmd="$_cmd --run-clink"
    [[ -n "$WIZ_CLINK_UMI_LEN" ]] && _cmd="$_cmd --clink-umi-len $WIZ_CLINK_UMI_LEN"
    [[ "$WIZ_CLINK_FDR" != "0.05" ]] && _cmd="$_cmd --clink-fdr $WIZ_CLINK_FDR"
    [[ "$WIZ_CLINK_MIN_COV" != "5" ]] && _cmd="$_cmd --clink-min-cov $WIZ_CLINK_MIN_COV"
    [[ "$WIZ_CLINK_MULTI_MAP" == "true" ]] && _cmd="$_cmd --clink-multi-map"
    [[ "$WIZ_NO_MOTIF" == "true" ]] && _cmd="$_cmd --no-motif"
    [[ "$WIZ_FLANK" != "10" ]] && _cmd="$_cmd -f $WIZ_FLANK"
    CLIPITTY_EQUIV_CMD="$_cmd"
    echo -e "  ${WIZ_YELLOW}Equivalent command (copy to skip the wizard next time):${WIZ_NC}"
    echo "    $_cmd"
    echo ""

    prompt_yesno "  Start mapping with these settings?" "y" WIZ_CONFIRM

    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi

    # Export variables for main script
    export WIZ_INPUT_FILE WIZ_GENOME_INDEX WIZ_GENOME_FASTA WIZ_FILTER_REPEAT WIZ_NO_CHR_FILTER
    export WIZ_ALIGNER WIZ_THREADS WIZ_ALIGN_MISMATCHES WIZ_OUTPUT_NAME WIZ_ALIGNER_ARGS WIZ_UMI_LEN
    export WIZ_RUN_CIMS WIZ_RUN_CITS WIZ_CIMS_ITER WIZ_CIMS_FDR WIZ_CITS_PVAL WIZ_CITS_GAP
    export WIZ_RUN_CLINK WIZ_CLINK_UMI_LEN WIZ_CLINK_FDR WIZ_CLINK_MIN_COV WIZ_CLINK_MULTI_MAP
    export WIZ_NO_MOTIF WIZ_FLANK
    export CLIPITTY_EQUIV_CMD

    echo -e "${WIZ_GREEN}Starting mapping...${WIZ_NC}"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════════
# PEAKittyPeak WIZARD
# ═══════════════════════════════════════════════════════════════════════════════

run_wizard_peakittypeak() {
    _play_wizard_intro
    print_wizard_header "PEAKittyPeak Wizard"
    local _yn

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Input Directory
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 1: Input Directory"
    echo ""
    echo "  This tool requires a BED/ folder containing collapsed BED files."
    echo ""
    read -p "  Enter path to working directory (or press Enter for current): " WIZ_WORK_DIR || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }
    [[ -z "$WIZ_WORK_DIR" ]] && WIZ_WORK_DIR="."
    WIZ_WORK_DIR="${WIZ_WORK_DIR/#\~/$HOME}"

    if [[ ! -d "$WIZ_WORK_DIR" ]]; then
        echo -e "${WIZ_RED}  ✗ Directory not found: $WIZ_WORK_DIR${WIZ_NC}"
        return 1
    fi
    WIZ_BED_COUNT=$(ls "$WIZ_WORK_DIR"/BED/*.bed 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$WIZ_BED_COUNT" -eq 0 ]]; then
        echo -e "${WIZ_RED}  ✗ No BED/ folder or .bed files found in $WIZ_WORK_DIR${WIZ_NC}"
        return 1
    fi
    echo -e "${WIZ_GREEN}  ✓ Found $WIZ_BED_COUNT .bed files in BED/ folder${WIZ_NC}"
    echo ""
    echo -e "${WIZ_GREEN}✓ Required inputs are ready!${WIZ_NC}"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Aggregation Mode
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 2: Aggregation Mode"
    echo ""
    echo "  Aggregate mode: pool all BED files into one peak-calling run."
    echo "  Individual mode: call peaks on each BED file separately."
    echo ""
    prompt_select "Which mode?" WIZ_AGG_SEL \
        "Individual (default) — per-sample peak calling" \
        "Aggregate — pool all BED files"
    WIZ_AGGREGATE="false"
    [[ "$WIZ_AGG_SEL" == "2" ]] && WIZ_AGGREGATE="true"

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Settings tier
    # ─────────────────────────────────────────────────────────────────────────
    print_section "STEP 3: Settings"
    echo ""
    prompt_select "How would you like to proceed?" WIZ_SETTINGS_MODE \
        "Use default settings" \
        "Modify program settings (advanced)"

    # Defaults (match PEAKittyPeak.sh)
    WIZ_PEAK_CALLER="homer"
    WIZ_PEAK_DIST="50"
    WIZ_PEAK_SIZE="20"
    WIZ_FRAG_LEN="25"
    WIZ_OUTPUT_NAME="Combined"
    WIZ_PEAK_CALLER_ARGS=""
    WIZ_GROUPS_FILE=""
    WIZ_GROUP_PEAKS_FILE=""
    WIZ_CTK_DIR=""
    WIZ_CIMS_FDR="0.05"
    WIZ_CITS_PVAL="0.05"

    if [[ "$WIZ_SETTINGS_MODE" == "2" ]]; then
        print_doc_box "PEAKittyPeak" \
            "PEAKittyPeak.sh --help" \
            "HOMER findPeaks: http://homer.ucsd.edu/homer/ngs/peaks.html"

        # ─── Peak Caller ───────────────────────────────────────────────────
        print_section "PEAK CALLER"
        echo ""
        local _sel
        prompt_select "Select peak caller:" _sel \
            "HOMER findPeaks (default)" \
            "CTK tag2peak.pl"
        if [[ "$_sel" == "2" ]]; then
            WIZ_PEAK_CALLER="ctk"
            echo -e "  ${WIZ_YELLOW}CTK tag2peak defaults:${WIZ_NC} -big -ss --valley-seeking"
            prompt_value "  Enter gap (-gap)" "50" WIZ_PEAK_DIST "int"
        else
            echo -e "  ${WIZ_YELLOW}HOMER defaults:${WIZ_NC} -style factor · -L 2 · -localSize 1000 · -minDist 50 · -size 20 · -fragLength 25"
            prompt_value "  Enter min distance between peaks" "50" WIZ_PEAK_DIST "int"
            prompt_value "  Enter peak size" "20" WIZ_PEAK_SIZE "int"
            prompt_value "  Enter fragment length" "25" WIZ_FRAG_LEN "int"
        fi
        prompt_value "  Output base name" "Combined" WIZ_OUTPUT_NAME
        echo "  Enter additional peak-caller arguments (optional, blank to skip):"
        read -p "  > " WIZ_PEAK_CALLER_ARGS || { echo -e "${WIZ_RED}  ✗ Input stream closed — aborting wizard.${WIZ_NC}" >&2; exit 1; }

        # ─── Grouping ──────────────────────────────────────────────────────
        print_section "GROUPING"
        echo ""
        prompt_yesno "  Provide a groups file for aggregation metrics?" "n" _yn
        if [[ "$_yn" == "y" ]]; then
            prompt_file_path "  Enter path to groups file (-g): " WIZ_GROUPS_FILE || return 1
        fi
        prompt_yesno "  Use --group-peaks (peak-level grouping mode)?" "n" _yn
        if [[ "$_yn" == "y" ]]; then
            prompt_file_path "  Enter path to group-peaks file: " WIZ_GROUP_PEAKS_FILE || return 1
        fi

        # ─── CTK Annotation ────────────────────────────────────────────────
        print_section "CTK CIMS/CITS ANNOTATION (optional)"
        echo ""
        echo "  If you have a CTK analysis output directory, the peak matrix"
        echo "  can be annotated with CIMS/CITS site counts."
        prompt_yesno "  Annotate peaks with CTK CIMS/CITS results?" "n" _yn
        if [[ "$_yn" == "y" ]]; then
            prompt_dir_path "  Enter path to CTK output directory (--ctk-dir): " WIZ_CTK_DIR || return 1
            prompt_value "  CIMS FDR threshold" "0.05" WIZ_CIMS_FDR "float"
            prompt_value "  CITS p-value threshold" "0.05" WIZ_CITS_PVAL "float"
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
    local _agg="Individual"; [[ "$WIZ_AGGREGATE" == "true" ]] && _agg="Aggregate"
    printf "  │ %-20s %-38s │\n" "Mode:" "$_agg"
    echo "  ├─────────────────────────────────────────────────────────────┤"
    printf "  │ %-20s %-38s │\n" "Peak Caller:" "$(echo "$WIZ_PEAK_CALLER" | tr '[:lower:]' '[:upper:]')"
    printf "  │ %-20s %-38s │\n" "Peak Distance/Gap:" "$WIZ_PEAK_DIST"
    if [[ "$WIZ_PEAK_CALLER" == "homer" ]]; then
        printf "  │ %-20s %-38s │\n" "Peak Size:" "$WIZ_PEAK_SIZE"
        printf "  │ %-20s %-38s │\n" "Fragment Length:" "$WIZ_FRAG_LEN"
    fi
    printf "  │ %-20s %-38s │\n" "Output Name:" "$WIZ_OUTPUT_NAME"
    [[ -n "$WIZ_PEAK_CALLER_ARGS" ]] && printf "  │ %-20s %-38s │\n" "Caller Args:" "$WIZ_PEAK_CALLER_ARGS"
    if [[ -n "$WIZ_GROUPS_FILE" || -n "$WIZ_GROUP_PEAKS_FILE" ]]; then
        echo "  ├─────────────────────────────────────────────────────────────┤"
        [[ -n "$WIZ_GROUPS_FILE" ]] && printf "  │ %-20s %-38s │\n" "Groups File:" "$(basename "$WIZ_GROUPS_FILE")"
        [[ -n "$WIZ_GROUP_PEAKS_FILE" ]] && printf "  │ %-20s %-38s │\n" "Group Peaks File:" "$(basename "$WIZ_GROUP_PEAKS_FILE")"
    fi
    if [[ -n "$WIZ_CTK_DIR" ]]; then
        echo "  ├─────────────────────────────────────────────────────────────┤"
        printf "  │ %-20s %-38s │\n" "CTK Dir:" "$(basename "$WIZ_CTK_DIR")"
        printf "  │ %-20s %-38s │\n" "CIMS FDR:" "$WIZ_CIMS_FDR"
        printf "  │ %-20s %-38s │\n" "CITS p-value:" "$WIZ_CITS_PVAL"
    fi
    echo "  └─────────────────────────────────────────────────────────────┘"
    echo ""

    # Equivalent command line
    local _cmd="PEAKittyPeak.sh -i \"$WIZ_WORK_DIR/BED\""
    [[ "$WIZ_AGGREGATE" == "true" ]] && _cmd="$_cmd --aggregate"
    [[ "$WIZ_PEAK_CALLER" != "homer" ]] && _cmd="$_cmd --peak-caller $WIZ_PEAK_CALLER"
    [[ "$WIZ_PEAK_DIST" != "50" ]] && _cmd="$_cmd -p $WIZ_PEAK_DIST"
    [[ "$WIZ_PEAK_CALLER" == "homer" && "$WIZ_PEAK_SIZE" != "20" ]] && _cmd="$_cmd -z $WIZ_PEAK_SIZE"
    [[ "$WIZ_PEAK_CALLER" == "homer" && "$WIZ_FRAG_LEN" != "25" ]] && _cmd="$_cmd -f $WIZ_FRAG_LEN"
    [[ "$WIZ_OUTPUT_NAME" != "Combined" ]] && _cmd="$_cmd -n \"$WIZ_OUTPUT_NAME\""
    [[ -n "$WIZ_PEAK_CALLER_ARGS" ]] && _cmd="$_cmd --peak-caller-args \"$WIZ_PEAK_CALLER_ARGS\""
    [[ -n "$WIZ_GROUPS_FILE" ]] && _cmd="$_cmd -g \"$WIZ_GROUPS_FILE\""
    [[ -n "$WIZ_GROUP_PEAKS_FILE" ]] && _cmd="$_cmd --group-peaks \"$WIZ_GROUP_PEAKS_FILE\""
    [[ -n "$WIZ_CTK_DIR" ]] && _cmd="$_cmd --ctk-dir \"$WIZ_CTK_DIR\""
    [[ "$WIZ_CIMS_FDR" != "0.05" ]] && _cmd="$_cmd --cims-fdr $WIZ_CIMS_FDR"
    [[ "$WIZ_CITS_PVAL" != "0.05" ]] && _cmd="$_cmd --cits-pval $WIZ_CITS_PVAL"
    CLIPITTY_EQUIV_CMD="$_cmd"
    echo -e "  ${WIZ_YELLOW}Equivalent command (copy to skip the wizard next time):${WIZ_NC}"
    echo "    $_cmd"
    echo ""

    prompt_yesno "  Start peak calling with these settings?" "y" WIZ_CONFIRM
    if [[ "$WIZ_CONFIRM" == "n" ]]; then
        echo "Configuration aborted."
        return 1
    fi

    export WIZ_WORK_DIR WIZ_BED_COUNT WIZ_AGGREGATE
    export WIZ_PEAK_CALLER WIZ_PEAK_CALLER_ARGS
    export WIZ_PEAK_DIST WIZ_PEAK_SIZE WIZ_FRAG_LEN WIZ_OUTPUT_NAME
    export WIZ_GROUPS_FILE WIZ_GROUP_PEAKS_FILE
    export WIZ_CTK_DIR WIZ_CIMS_FDR WIZ_CITS_PVAL
    export CLIPITTY_EQUIV_CMD

    echo -e "${WIZ_GREEN}Starting peak calling...${WIZ_NC}"
    return 0
}
