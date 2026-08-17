#!/bin/bash

# lib/dedup.sh - FASTQ Deduplication Module for CLIPittyClip
#
# Standalone module: source this file independently of modules.sh for testing.
# Requires: lib/utils.sh (log_info, log_warning, log_error must be defined)
#
# Public API:
#   run_dedup <input.fastq[.gz]> <output.fastq> [compress_output=false]
#     Emits plain .fastq by default; pass "true" as $3 to gzip the output.
#
# Internal helpers (also used by run_fastp in modules.sh):
#   _fastq_collapse_core <input.fastq> <output.fastq>
#   detect_eclip_umi_length <input.fastq.gz> [user_umi_len]
#   reformat_eclip_umi_to_sequence <input.fastq.gz> <output.fastq.gz> <umi_len>
#   strip_eclip_barcode <input.fastq.gz> <output.fastq.gz> <umi_len>

# ═══════════════════════════════════════════════════════════════════════════
# Core Deduplication Engine
# ═══════════════════════════════════════════════════════════════════════════

# _fastq_collapse_core — dedup a plain (uncompressed) FASTQ
#
# Engine: hash-based (fastq_collapse_hash.py) by default — fastest, RAM
# scales ~169MB per 1M reads (measured). Pass LOW_MEMORY="true" (set via
# CLIPittyClip's --low-memory flag) to force the sort-based engine instead
# (fastq_collapse_sort.sh) — flat ~30-40MB regardless of size, ~2x faster
# than CTK's fastq2collapse.pl at every scale tested, but slower than hash.
# Above 30M reads with the hash engine still selected, a warning is printed
# suggesting --low-memory (hash RAM use at that size is untested/uncapped).
# Falls back to the sort engine unconditionally if python3/hash script are
# unavailable, and to an inline awk|sort|uniq pipeline only if
# fastq_collapse_sort.sh itself is missing.
#
# Args:    $1 = input.fastq (plain, not gzipped)
#          $2 = output.fastq (plain, not gzipped)
# Returns: 0 on success, non-zero on failure
_fastq_collapse_core() {
    local input="$1"
    local output="$2"
    local hash_script sort_script
    hash_script="$(dirname "${BASH_SOURCE[0]}")/fastq_collapse_hash.py"
    sort_script="$(dirname "${BASH_SOURCE[0]}")/fastq_collapse_sort.sh"

    local n_reads
    n_reads=$(( $(wc -l < "$input") / 4 ))
    local warn_threshold=30000000

    # Large-library warning.
    #
    # This fires per sample, so in batch mode it used to print a multi-line
    # block for every large library regardless of -v, burying the actual
    # progress output. Now:
    #   -v given     → full explanation inline, as before
    #   -v not given → a single yellow [!] marker on the progress line, so the
    #                  condition is still visible without dominating the output
    #   either way   → recorded to the log, and to a run-scoped tally that the
    #                  final summary reports (see print_dedup_warning_summary)
    if [[ "$n_reads" -gt "$warn_threshold" && "${LOW_MEMORY:-false}" != "true" ]]; then
        log_warning "Library has ${n_reads} reads (>${warn_threshold}); hash-based dedup RAM use is unbounded at this scale. Consider re-running with --low-memory."
        if [[ -n "${DEDUP_WARN_FILE:-}" ]]; then
            printf '%s\t%s\n' "${BASENAME:-$(basename "$input")}" "$n_reads" >> "$DEDUP_WARN_FILE"
        fi
        if [[ "${VERBOSE:-false}" == "true" ]]; then
            echo -e "\n      [WARNING: ${n_reads} reads > 30M — consider --low-memory to use the flat-RAM sort engine]" >&2
        else
            echo -ne " \033[0;33m[!]\033[0m" >&2
        fi
    fi

    # Engine-selection detail is diagnostic, not progress — verbose only.
    _dedup_note () { [[ "${VERBOSE:-false}" == "true" ]] && echo -e "\n      $1" >&2; return 0; }

    if [[ "${LOW_MEMORY:-false}" != "true" ]] && [[ -f "$hash_script" ]] && command -v python3 &>/dev/null; then
        log_info "Dedup engine: hash-based (fastq_collapse_hash.py) — ${n_reads} reads"
        _dedup_note "[dedup: hash-based engine, ${n_reads} reads]"
        python3 "$hash_script" "$input" "$output" 2>> "${LOG_FILE:-/dev/null}"
    elif [[ -x "$sort_script" ]]; then
        log_info "Dedup engine: sort-based (fastq_collapse_sort.sh) — ${n_reads} reads"
        _dedup_note "[dedup: sort-based engine, ${n_reads} reads]"
        "$sort_script" "$input" "$output" 2>> "${LOG_FILE:-/dev/null}"
    else
        log_info "Dedup engine: sort-based inline fallback (awk | sort | uniq) — ${n_reads} reads"
        _dedup_note "[dedup: sort-based inline fallback, ${n_reads} reads]"
        # Replicates fastq2collapse.pl column order exactly:
        #   paste order: ID | QUAL | SEQ  →  sort -k3 (by SEQ)
        #   uniq -f2 -c  →  fields: $1=count $2=ID $3=QUAL $4=SEQ
        #   output: $2"#"$1 = @id#COUNT
        awk 'NR%4==1{id=$1} NR%4==2{seq=$0} NR%4==0{print id"\t"$0"\t"seq}' "$input" \
            | sort -k3 \
            | uniq -f2 -c \
            | awk '{print $2"#"$1"\n"$4"\n+\n"$3}' > "$output"
    fi
}

# ═══════════════════════════════════════════════════════════════════════════
# Public: run_dedup
# ═══════════════════════════════════════════════════════════════════════════

# run_dedup — deduplicate a (gzipped or plain) FASTQ
# Accepts .fastq or .fastq.gz input (auto-decompresses if needed).
# Default output: plain .fastq (no gzip overhead on intermediate files).
# Pass compress_output="true" as $3 only when a gzipped handoff product is needed.
# Returns 0 on success, 1 on failure (caller handles fallback + messaging).
#
# Args: $1 = input file        (.fastq.gz or .fastq)
#       $2 = output file       (.fastq by default; .fastq.gz when compress_output=true)
#       $3 = compress_output   (optional; "true" to gzip the output; default: "false")
run_dedup() {
    local input_file="$1"
    local output_file="$2"
    local compress_output="${3:-false}"

    log_info "Deduplication: starting..."

    # Determine plain output path (strip .gz suffix when caller wants compression)
    local plain_out
    if [[ "$compress_output" == "true" ]]; then
        plain_out="${output_file%.gz}"
    else
        plain_out="$output_file"
    fi

    local exit_code
    if [[ "$input_file" == *.gz ]]; then
        local temp_in="${plain_out}_input_temp.fastq"
        gzip -dc "$input_file" > "$temp_in"
        _fastq_collapse_core "$temp_in" "$plain_out"
        exit_code=$?
        rm -f "$temp_in"
    else
        _fastq_collapse_core "$input_file" "$plain_out"
        exit_code=$?
    fi

    if [[ $exit_code -eq 0 && -s "$plain_out" ]]; then
        if [[ "$compress_output" == "true" ]]; then
            gzip -c "$plain_out" > "$output_file"
            rm -f "$plain_out"
            log_info "Deduplication complete (gzipped): $output_file"
        else
            log_info "Deduplication complete: $plain_out"
        fi
        return 0
    else
        rm -f "$plain_out" "$output_file"
        log_warning "Deduplication failed or produced empty output."
        return 1
    fi
}

# ═══════════════════════════════════════════════════════════════════════════
# eCLIP UMI Helpers (called by run_fastp in modules.sh)
# ═══════════════════════════════════════════════════════════════════════════

# detect_eclip_umi_length — infer UMI length from first read header
# Input header format: @UMI:REST_OF_ID
# Returns: detected length (echoes integer)
#
# Args: $1 = input.fastq.gz
#       $2 = user_umi_len (optional, for mismatch warning; 0 = not specified)
detect_eclip_umi_length() {
    local input_fastq="$1"
    local user_umi_len="${2:-0}"

    local first_header
    if [[ "$input_fastq" == *.gz ]]; then
        first_header=$(gunzip -c "$input_fastq" | head -1)
    else
        first_header=$(head -1 "$input_fastq")
    fi

    local id_part="${first_header%% *}"   # Remove comment after space
    local id_no_at="${id_part#@}"         # Remove leading @
    local umi="${id_no_at%%:*}"           # Get part before first colon
    local detected_len=${#umi}

    if [[ "$user_umi_len" -gt 0 && "$detected_len" -ne "$user_umi_len" ]]; then
        log_warning "Detected UMI length ($detected_len) differs from specified ($user_umi_len). Using detected."
    fi

    log_info "eCLIP UMI length detected: $detected_len nt"
    echo "$detected_len"
}

# reformat_eclip_umi_to_sequence — move UMI from header into sequence (CTK workflow)
# Required before collapse so UMI participates in duplicate detection.
#
# Input:  @UMI:READ_ID  /  SEQUENCE  /  +  /  QUAL
# Output: @READ_ID      /  UMI+SEQUENCE  /  +  /  UMI_QUAL(Phred40)+QUAL
#
# Args: $1 = input.fastq.gz
#       $2 = output.fastq.gz
#       $3 = umi_len (integer)
reformat_eclip_umi_to_sequence() {
    local input_fastq="$1"
    local output_fastq="$2"
    local umi_len="$3"

    log_info "Moving eCLIP UMI from header to sequence (CTK documentation workflow)..."
    log_info "UMI length: $umi_len nt"

    # Generate quality string for UMI (I = Phred 40, high quality)
    local umi_qual
    umi_qual=$(printf 'I%.0s' $(seq 1 "$umi_len"))

    # Write to temp file first, then gzip (avoids macOS pipe buffering issues)
    local temp_output="${output_fastq%.gz}"

    local read_cmd="cat"
    if [[ "$input_fastq" == *.gz ]]; then read_cmd="gzip -dc"; fi

    $read_cmd "$input_fastq" | awk -v umi_len="$umi_len" -v umi_qual="$umi_qual" '
    BEGIN { OFS="" }
    {
        if (NR % 4 == 1) {
            # Header line: @UMI:READ_ID COMMENT
            n = split($0, parts, " ")
            id_part = parts[1]
            comment = ""
            if (n > 1) { for (i=2; i<=n; i++) comment = comment " " parts[i] }

            id_no_at = substr(id_part, 2)
            colon_pos = index(id_no_at, ":")
            if (colon_pos > 0) {
                umi = substr(id_no_at, 1, colon_pos-1)
                rest = substr(id_no_at, colon_pos+1)
                saved_umi = umi
                print "@" rest comment
            } else {
                saved_umi = ""
                print $0
            }
        } else if (NR % 4 == 2) {
            # Sequence line: prepend UMI
            print saved_umi $0
        } else if (NR % 4 == 0) {
            # Quality line: prepend UMI quality scores
            print umi_qual $0
        } else {
            # + line
            print $0
        }
    }' > "$temp_output"

    if [[ -s "$temp_output" ]]; then
        gzip -f "$temp_output"
    fi

    if [[ -s "$output_fastq" ]]; then
        log_info "UMI moved to sequence: $output_fastq"
    else
        log_error "UMI reformat failed - output is empty"
        return 1
    fi
}

# strip_eclip_barcode — extract UMI from sequence and attach to read ID after collapse
# Uses CTK stripBarcode.pl.  Result: @READ_ID#count#UMI (CTK-compatible format)
#
# Args: $1 = input.fastq.gz
#       $2 = output.fastq.gz
#       $3 = umi_len (integer)
strip_eclip_barcode() {
    local input_fastq="$1"
    local output_fastq="$2"
    local umi_len="$3"

    log_info "Stripping UMI from sequence with stripBarcode.pl (len=$umi_len)..."

    stripBarcode.pl -format fastq -len "$umi_len" "$input_fastq" - 2>> "${LOG_FILE}" | gzip > "$output_fastq"

    if [[ -s "$output_fastq" ]]; then
        log_info "UMI stripped and attached to read ID: $output_fastq"
    else
        log_error "stripBarcode.pl failed - output is empty"
        return 1
    fi
}

# ── Run-scoped tally of large-library dedup warnings ──────────────────────────
# Each sample runs as its own child process in batch mode, so the warnings
# cannot be accumulated in a shell variable. _fastq_collapse_core appends to
# $DEDUP_WARN_FILE instead, and this prints the collected result once at the
# end of the run — to the console and, via console_msg, to the detailed log.
#
# Called from the [COMPLETE] block. Safe to call when nothing warned.
print_dedup_warning_summary () {
    [[ -n "${DEDUP_WARN_FILE:-}" && -s "${DEDUP_WARN_FILE:-/nonexistent}" ]] || return 0
    local n
    n=$(wc -l < "$DEDUP_WARN_FILE" | tr -d ' ')
    console_msg "\n\033[0;33m[DEDUP MEMORY NOTICE]\033[0m"
    console_msg "  ${n} librar$([[ $n -eq 1 ]] && echo y || echo ies) exceeded 30M reads and used the hash-based dedup engine,"
    console_msg "  whose RAM use scales with unique reads and is not capped:"
    while IFS=$'\t' read -r s r; do
        console_msg "    $(printf '%-34s %15s reads' "$s" "$(printf "%'d" "$r" 2>/dev/null || echo "$r")")"
    done < "$DEDUP_WARN_FILE"
    console_msg "  If this run was killed during deduplication, re-run with --low-memory"
    console_msg "  to use the flat-RAM sort engine instead."
    rm -f "$DEDUP_WARN_FILE"
}
