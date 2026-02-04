#!/bin/bash

# lib/modules.sh - Analysis modules for CLIPittyClip

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ═══════════════════════════════════════════════════════════════════════════
# Group File Parsing for CIMS/CITS Aggregation
# ═══════════════════════════════════════════════════════════════════════════

# Check if GNU parallel is installed
has_gnu_parallel() {
    command -v parallel &>/dev/null && parallel --version 2>&1 | head -1 | grep -q "GNU parallel"
}

# Calculate optimal number of parallel jobs based on available RAM and file size
# Arguments: $1 = user_threads, $2 = input_file, $3 = min_ram_per_job_gb (default: 2)
# Output: echoes optimal job count
calculate_optimal_parallel_jobs() {
    local user_threads="$1"
    local input_file="$2"
    local min_ram_per_job="${3:-2}"  # Default: 2GB per job minimum
    
    # Get available RAM in GB (use 'available' column from free)
    local available_ram_gb
    available_ram_gb=$(free -g 2>/dev/null | awk '/^Mem:/ {print $7}')
    
    # Fallback if free -g fails
    if [[ -z "$available_ram_gb" ]] || [[ "$available_ram_gb" -eq 0 ]]; then
        available_ram_gb=8  # Conservative default
    fi
    
    # Get file size and estimate RAM per job
    # Rule of thumb: jobs processing large files need more RAM
    local file_lines=0
    if [[ -f "$input_file" ]]; then
        file_lines=$(wc -l < "$input_file" 2>/dev/null || echo 0)
    fi
    
    # Estimate: for files with >1M lines, use more RAM per job
    local ram_per_job=$min_ram_per_job
    if [[ "$file_lines" -gt 1000000 ]]; then
        ram_per_job=4  # Large files need 4GB per job
    elif [[ "$file_lines" -gt 5000000 ]]; then
        ram_per_job=6  # Very large files need 6GB per job
    fi
    
    # Calculate RAM-based job limit
    local ram_based_jobs=$((available_ram_gb / ram_per_job))
    ram_based_jobs=$((ram_based_jobs > 0 ? ram_based_jobs : 1))  # At least 1
    
    # Final = minimum of user threads and RAM-based limit
    local optimal_jobs=$((user_threads < ram_based_jobs ? user_threads : ram_based_jobs))
    optimal_jobs=$((optimal_jobs > 0 ? optimal_jobs : 1))  # At least 1
    
    echo "$optimal_jobs"
}

# parse_groups_file - Parse groups.txt and output sample→group mapping
# Input: groups_file path, output_map temp file path
# Format: sample_name<TAB>group_name (lines starting with # are comments)
parse_groups_file() {
    local groups_file="$1"
    local output_map="$2"
    
    log_info "Parsing groups file: $groups_file"
    
    > "$output_map"  # Clear output file
    
    while IFS=$'\t' read -r sample group || [[ -n "$sample" ]]; do
        # Skip comments and empty lines
        [[ "$sample" =~ ^#.*$ ]] && continue
        [[ -z "$sample" ]] && continue
        
        # Strip common extensions
        sample="${sample%.fastq.gz}"
        sample="${sample%.fq.gz}"
        sample="${sample%.fastq}"
        sample="${sample%.fq}"
        
        # Write mapping
        echo -e "${sample}\t${group}" >> "$output_map"
    done < "$groups_file"
    
    local group_count=$(cut -f2 "$output_map" | sort -u | wc -l | tr -d ' ')
    local sample_count=$(wc -l < "$output_map" | tr -d ' ')
    log_info "Parsed $sample_count samples into $group_count groups"
}

# ═══════════════════════════════════════════════════════════════════════════

# 0. Demultiplexing
function run_demultiplexing {
    local input_fastq="$1"
    local barcode_file="$2"
    local run_sample_size="$3"
    # Note: dedup_mode parameter ignored - dedup is now handled by caller
    local mismatches="${DEMUX_MISMATCHES:-1}" # Default to 1 if unset/passed global

    log_info "Starting Demultiplexing with cutadapt..."
    log_info "Barcode File: $barcode_file"
    log_info "Allowed Mismatches: $mismatches"
    
    # Input is either the original pooled file or already-deduplicated file
    local work_input="$input_fastq"
    # We call the script relative to the repo root
    local checker_script="${SCRIPT_DIR}/check_barcodes.sh"
    
    if [[ ! -x "$checker_script" ]]; then
        # Fallback if not executable or found (chmod just in case)
        chmod +x "$checker_script" 2>/dev/null
    fi

    "$checker_script" -f "$barcode_file" -m "$mismatches"
    if [ $? -ne 0 ]; then
        log_error "Barcode collisions detected. Aborting pipeline to prevent data mix-up."
        exit 1
    fi
    log_info "Barcode safety check PASSED."

    # 2. Calculate average barcode length to set cutadapt error rate

    # 2. Calculate average barcode length to set cutadapt error rate
    # cutadapt -e is a rate (0.1 = 10%). 
    # rate = mismatches / length
    # We grep the first barcode to estimate length (assuming uniform)
    local first_seq=$(awk '{print $2; exit}' "$barcode_file")
    local bc_len=${#first_seq}
    
    # Calculate rate using awk for floating point
    local error_rate=$(awk "BEGIN {print $mismatches / $bc_len}")
    
    # Cap strictness if 0 errors
    if [[ "$mismatches" -eq 0 ]]; then
        error_rate=0
    fi
    
    log_info "Calculated cutadapt error rate: $error_rate ($mismatches errors in ${bc_len}bp)"

    # Create demux output directory
    mkdir -p demux_fastq
    
    # Convert barcodes for cutadapt
    local fasta_barcodes="barcodes.fasta"
    awk '{print ">"$1"\n"$2}' "$barcode_file" > "$fasta_barcodes"

    local cmd="cutadapt \
        -e $error_rate --no-indels \
        -m 1 \
        -g file:$fasta_barcodes \
        -o \"demux_fastq/{name}.fastq.gz\" \
        $work_input \
        -j ${THREADS:-1}"
    
    log_info "Running demultiplexing..."
    execute_cmd "$cmd"

    # Check outputs
    count=$(ls demux_fastq/*.fastq.gz 2>/dev/null | wc -l)
    if [ "$count" -eq 0 ]; then
        log_error "Demultiplexing failed. No output files created."
        exit 1
    fi
    log_info "Demultiplexing complete. Created $count sample files in 'demux_fastq/'."
    
    # Cleanup Deduplicated Temp File if it exists
    if [[ "$work_input" != "$input_fastq" ]] && [[ -f "$work_input" ]]; then
        log_info "Deleting temporary deduplicated file: $work_input"
        rm -f "$work_input"
    fi
}

# Filter BAM to canonical chromosomes only (chr1-22, X, Y, M)
# Removes contigs like GL000220.1, KI270733.1, etc.
# This reduces data size and improves downstream processing speed
filter_canonical_chromosomes() {
    local input_bam="$1"
    local output_bam="$2"
    
    log_info "Filtering to canonical chromosomes (chr1-22, X, Y, M)..."
    
    # Create list of canonical chromosomes
    local chr_list="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10"
    chr_list+=" chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
    chr_list+=" chr20 chr21 chr22 chrX chrY chrM"
    
    # Count reads before filtering
    local before_count=$(samtools view -c "$input_bam")
    
    # Filter BAM to canonical chromosomes
    samtools view -b "$input_bam" $chr_list > "$output_bam" 2>> "${LOG_FILE}"
    
    if [[ $? -eq 0 && -s "$output_bam" ]]; then
        samtools index "$output_bam"
        local after_count=$(samtools view -c "$output_bam")
        local filtered=$((before_count - after_count))
        log_info "Chromosome filtering complete: $after_count reads kept, $filtered reads on contigs removed"
    else
        log_warning "Chromosome filtering failed. Using unfiltered BAM."
        cp "$input_bam" "$output_bam"
        samtools index "$output_bam"
    fi
}

# Detect eCLIP UMI length from first read of FASTQ file
# Returns the detected UMI length (e.g., 5 or 10)
detect_eclip_umi_length() {
    local input_fastq="$1"
    local user_umi_len="${2:-0}"
    
    # Get first read header and extract UMI (before first colon)
    local first_header
    if [[ "$input_fastq" == *.gz ]]; then
        first_header=$(gunzip -c "$input_fastq" | head -1)
    else
        first_header=$(head -1 "$input_fastq")
    fi
    
    # Parse: @UMI:REST_OF_ID -> extract UMI
    local id_part="${first_header%% *}"  # Remove comment after space
    local id_no_at="${id_part#@}"        # Remove leading @
    local umi="${id_no_at%%:*}"          # Get part before first colon
    local detected_len=${#umi}
    
    if [[ "$user_umi_len" -gt 0 && "$detected_len" -ne "$user_umi_len" ]]; then
        log_warning "Detected UMI length ($detected_len) differs from specified ($user_umi_len). Using detected."
    fi
    
    log_info "eCLIP UMI length detected: $detected_len nt"
    echo "$detected_len"
}

# Reformat eCLIP: Move UMI from header to sequence (CTK documentation workflow)
# This is required for correct fastq2collapse.pl behavior
# Input:  @UMI:READ_ID
#         SEQUENCE
# Output: @READ_ID
#         UMI+SEQUENCE (UMI prepended)
#         +
#         UMI_QUAL+QUAL (high quality prepended)
reformat_eclip_umi_to_sequence() {
    local input_fastq="$1"
    local output_fastq="$2"
    local umi_len="$3"
    
    log_info "Moving eCLIP UMI from header to sequence (CTK documentation workflow)..."
    log_info "UMI length: $umi_len nt"
    
    # Generate quality string for UMI (I = Phred 40, high quality)
    local umi_qual=$(printf 'I%.0s' $(seq 1 $umi_len))
    
    # Write to temp file first, then gzip (avoids macOS pipe buffering issues)
    local temp_output="${output_fastq%.gz}"
    
    gunzip -c "$input_fastq" | awk -v umi_len="$umi_len" -v umi_qual="$umi_qual" '
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
                # Store UMI for sequence line, output clean header
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
    
    # Explicitly gzip the temp file (creates .gz, removes original)
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

# Strip UMI from sequence and attach to read ID after collapse
# Uses CTK stripBarcode.pl to extract UMI from 5' end of sequence
# Result: @READ_ID#count#UMI (CTK-compatible format)
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

# 1. Preprocessing with fastp
run_preprocessing() {
    local input_file="$1"
    local output_prefix="$2"
    local umi_len="$3" 
    local adapter3="$4"
    local threads="$5"
    local sample_size="$6"
    local dedup_mode="$7"
    local eclip_mode="${8:-false}"  # eCLIP mode: UMI in header, use adapter FASTA
    
    update_status_first "Preprocessing"
    log_info "Starting preprocessing..."
    
    # Get path to eCLIP adapter FASTA (same dir as this script)
    local script_dir="$(dirname "${BASH_SOURCE[0]}")"
    local eclip_adapters_fasta="${script_dir}/eclip_adapters.fa"
    
    #==========================================================================
    # eCLIP MODE: CTK documentation workflow
    # 1. UMI header → sequence
    # 2. fastp (adapter trim, no --umi)
    # 3. fastq2collapse.pl (collapse with UMI in sequence)
    # 4. stripBarcode.pl (extract UMI to header after count)
    #==========================================================================
    if [[ "$eclip_mode" == "true" ]]; then
        log_info "eCLIP mode: Preprocessing workflow (collapse before trim)"
        
        # Step 1: Detect and validate UMI length
        local detected_umi_len=$(detect_eclip_umi_length "$input_file" "$umi_len")
        umi_len="$detected_umi_len"
        
        # Step 2: Move UMI from header to sequence
        update_status_first "eCLIP Preprocessing"
        echo -ne "(UMI to Sequence"
        local umi_seq_file="${output_prefix}_umi_in_seq.fastq.gz"
        reformat_eclip_umi_to_sequence "$input_file" "$umi_seq_file" "$umi_len"
        if [[ ! -s "$umi_seq_file" ]]; then
            log_error "Failed to reformat eCLIP UMI to sequence"
            exit 1
        fi
        
        # Step 3: Collapse exact duplicates FIRST (with UMI in sequence)
        # This matches non-eCLIP flow where collapse happens before fastp
        echo -ne " → Collapsing"
        local collapsed_file="${output_prefix}_collapsed.fastq"
        local collapsed_file_gz="${output_prefix}_collapsed.fastq.gz"
        gzip -dc "$umi_seq_file" > "${output_prefix}_umi_temp.fastq"
        fastq2collapse.pl "${output_prefix}_umi_temp.fastq" "$collapsed_file" 2>> "${LOG_FILE}"
        if [[ ! -s "$collapsed_file" ]]; then
            log_error "fastq2collapse.pl failed"
            exit 1
        fi
        gzip -c "$collapsed_file" > "$collapsed_file_gz"
        rm -f "${output_prefix}_umi_temp.fastq" "$collapsed_file" "$umi_seq_file"
        
        # Step 4: Strip UMI from sequence, attach to header after count
        echo -ne " → Extract UMI"
        local stripped_file="${output_prefix}_stripped.fastq.gz"
        strip_eclip_barcode "$collapsed_file_gz" "$stripped_file" "$umi_len"
        rm -f "$collapsed_file_gz"
        if [[ ! -s "$stripped_file" ]]; then
            log_error "stripBarcode.pl failed"
            exit 1
        fi
        
        # Step 5: Adapter trimming with fastp (AFTER collapse for efficiency)
        echo -ne " → Adapter Trim) > "
        local final_file="${output_prefix}_cleaned.fastq.gz"
        local fastp_cmd="fastp -i ${stripped_file} -o ${final_file} \
            --thread ${threads} \
            --adapter_fasta ${eclip_adapters_fasta} \
            --length_required 20 \
            --cut_tail --cut_tail_mean_quality 5 \
            --overlap_len_require 1 \
            --html ${output_prefix}_fastp.html \
            --json ${output_prefix}_fastp.json"
        if [ "$sample_size" -gt 0 ]; then fastp_cmd+=" --reads_to_process $sample_size"; fi
        log_info "Running: $fastp_cmd"
        execute_cmd "$fastp_cmd"
        rm -f "$stripped_file"  # Cleanup intermediate
        
        # Check fastp succeeded
        if [[ ! -s "$final_file" ]]; then
            log_error "fastp failed to create cleaned file"
            exit 1
        fi
        
        log_info "eCLIP preprocessing complete: $final_file"
        log_info "Read ID format: READ#count#UMI (CTK-compatible)"
        return 0
    fi
    
    #==========================================================================
    # NON-eCLIP MODE: Standard workflow
    # 1. Dedup with fastq2collapse.pl (if enabled)
    # 2. fastp (with --umi for UMI extraction)
    #==========================================================================
    log_info "Standard mode: fastp with UMI extraction"
    
    # Define Fastp Command Function
    run_fastp_cmd() {
        local in_f="$1"
        local cmd="fastp -i ${in_f} -o ${output_prefix}_cleaned.fastq.gz \
            --thread ${threads} \
            --length_required 16 \
            --average_qual 30 \
            --html ${output_prefix}_fastp.html \
            --json ${output_prefix}_fastp.json"
        # Standard mode: single adapter and UMI extraction
        if [ -n "$adapter3" ]; then cmd+=" --adapter_sequence ${adapter3}"; fi
        if [ "$umi_len" -gt 0 ]; then cmd+=" --umi --umi_loc=read1 --umi_len=${umi_len} --umi_delim=#"; fi
        if [ "$sample_size" -gt 0 ]; then cmd+=" --reads_to_process $sample_size"; fi
        
        log_info "Running: $cmd"
        execute_cmd "$cmd"
    }

    local final_input="$input_file"
    local temp_dedup=""
    
    # 1. Deduplication (fastq2collapse.pl) - Tracks duplicate counts for CTK
    if [ "$dedup_mode" == "true" ]; then
        update_status "  Deduplicating (CTK)" 
        log_info "Deduplication: Running fastq2collapse.pl (tracks duplicate counts for tag2collapse)..."
        temp_dedup="${output_prefix}_dedup_temp.fastq"
        local temp_dedup_gz="${output_prefix}_dedup_temp.fastq.gz"
        
        if [[ "$input_file" == *.gz ]]; then
            local temp_input="${output_prefix}_input_temp.fastq"
            gzip -dc "$input_file" > "$temp_input"
            fastq2collapse.pl "$temp_input" "$temp_dedup" 2>> "${LOG_FILE}"
            rm -f "$temp_input"
        else
            fastq2collapse.pl "$input_file" "$temp_dedup" 2>> "${LOG_FILE}"
        fi
        
        if [[ $? -eq 0 && -s "$temp_dedup" ]]; then
             gzip -c "$temp_dedup" > "$temp_dedup_gz"
             rm -f "$temp_dedup"
             final_input="$temp_dedup_gz"
             temp_dedup="$temp_dedup_gz"
             log_info "Deduplication complete."
             update_status "  Preprocessing (Fastp)"
        else
             log_warning "fastq2collapse.pl failed. Reverting to original input."
             rm -f "$temp_dedup" "$temp_dedup_gz"
             temp_dedup=""
        fi
    fi

    # 2. Run Fastp
    run_fastp_cmd "$final_input"
    local exit_code=$?
    
    # Cleanup temp dedup file
    if [[ -n "$temp_dedup" && -f "$temp_dedup" ]]; then
        rm "$temp_dedup"
    fi

    if [ $exit_code -eq 0 ]; then
        log_info "Preprocessing complete."
    else
        log_error "fastp failed."
        exit 1
    fi
}

# 1b. ncRNA Pre-filtering with Bowtie2
# Filters out rRNA, tRNA, and other ncRNA reads before genome alignment
# Input: FASTQ from fastp
# Output: Unmapped reads (for genome alignment), Mapped reads (QC)
run_ncrna_filter() {
    local input_fastq="$1"
    local output_unmapped="$2"    # Reads that didn't map to ncRNA (continue to genome)
    local output_dir="$3"         # Directory for ncRNA mapping outputs
    local index_dir="$4"
    local threads="$5"
    local sample_name="$6"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    local ncrna_bam="${output_dir}/${sample_name}_ncrna.bam"
    local ncrna_stats="${output_dir}/${sample_name}_ncrna_stats.txt"
    
    update_status "ncRNA Filter"
    
    # Run Bowtie2 mapping to ncRNA index
    # --un-gz: write unmapped reads to file (these go to genome mapping)
    # Mapped reads are saved to BAM for QC
    local bt2_cmd="bowtie2 -x \"${index_dir}/ncrna\" \
        -U \"$input_fastq\" \
        --un-gz \"$output_unmapped\" \
        -p $threads \
        2> \"$ncrna_stats\" \
        | samtools view -bS - > \"$ncrna_bam\""
    
    log_info "Running ncRNA filter: $bt2_cmd"
    
    if execute_cmd "$bt2_cmd"; then
        # Index the BAM for potential downstream use
        samtools index "$ncrna_bam" 2>/dev/null || true
        
        # Extract alignment rate from stats and display on console
        local align_rate=$(grep "overall alignment rate" "$ncrna_stats" | grep -oE "[0-9]+\.[0-9]+%" || echo "N/A")
        local total_reads=$(grep "reads; of these:" "$ncrna_stats" | grep -oE "^[0-9]+" || echo "N/A")
        local aligned_reads=$(grep "aligned exactly 1 time" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        local multi_aligned=$(grep "aligned >1 times" "$ncrna_stats" | grep -oE "^[[:space:]]*[0-9]+" | tr -d ' ' || echo "0")
        local ncrna_reads=$((aligned_reads + multi_aligned))
        
        log_info "ncRNA alignment rate: $align_rate"
        # Note: Per-sample stats logged to file; summary table shown after batch
        
        return 0
    else
        log_error "ncRNA filtering failed"
        return 1
    fi
}

# 2. Mapping with STAR
run_mapping_star() {
    local input_fastq="$1"
    local output_prefix="$2"
    local genome_dir="$3"
    local threads="$4"
    local mismatch_max="$5"

    update_status "Mapping (STAR)"
    log_info "Starting mapping with STAR..."

    local read_command="gunzip -c"
    if [[ "$input_fastq" != *.gz ]]; then
        read_command="cat"
    fi

    # Create temp directory for STAR (prevents FIFO errors on exFAT/NTFS drives)
    local star_tmp="${TMPDIR:-/tmp}/star_$(basename "$output_prefix")_$$"
    mkdir -p "$star_tmp"
    log_info "STAR temp directory: $star_tmp"

    local cmd="STAR --runThreadN ${threads} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${input_fastq} \
        --readFilesCommand ${read_command} \
        --outFileNamePrefix ${output_prefix}. \
        --outTmpDir ${star_tmp}/STARtmp \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNmax ${mismatch_max} \
        --outSAMattributes NH HI AS nM NM MD \
        --alignEndsType EndToEnd \
        $ADV_ALIGNER_ARGS"

    log_info "Running: $cmd"
    execute_cmd "$cmd"
    local star_exit=$?
    
    # Cleanup temp directory (always, even on failure)
    rm -rf "$star_tmp"
    
    if [ $star_exit -ne 0 ]; then
        log_error "STAR mapping failed. Check the log file for details."
        exit 1
    fi
    
    samtools index "${output_prefix}.Aligned.sortedByCoord.out.bam"
    if [ $? -ne 0 ]; then
        log_error "samtools index failed. The BAM file might be empty or invalid."
        exit 1
    fi

    log_info "Mapping complete. Output: ${output_prefix}.Aligned.sortedByCoord.out.bam"
}


# 2c. Parse Alignment for CIMS/CITS
# parseAlignment.pl requires SAM input (not BAM)
# This function converts BAM→SAM, optionally runs calmd for MD tags, then parses
run_parse_alignment() {
    local bam_file="$1"
    local output_bed="$2"
    local mutation_file="$3"
    local genome_index="$4"
    
    update_status "Processing Alignment"
    log_info "Parsing alignment for CIMS/CITS..."
    log_info "Input BAM: $bam_file"
    
    # Step 1: Convert BAM to SAM (parseAlignment.pl requires SAM input)
    local sam_file="${bam_file%.bam}.sam"
    log_info "Converting BAM to SAM for parseAlignment.pl..."
    samtools view -h "$bam_file" > "$sam_file"
    
    if [[ ! -s "$sam_file" ]]; then
        log_error "BAM to SAM conversion failed or empty output: $sam_file"
        return 1
    fi
    
    local work_sam="$sam_file"
    
    # Step 2: (Optional) Run samtools calmd for MD tag standardization
    local ref_fasta=$(find "$genome_index" -name "*.fa" -o -name "*.fasta" 2>/dev/null | head -n 1)
    
    if [[ -n "$ref_fasta" ]]; then
        log_info "Reference FASTA found: $ref_fasta"
        log_info "Running 'samtools calmd' to standardize MD tags..."
        
        local calmd_sam="${sam_file%.sam}_calmd.sam"
        # calmd can take BAM input and output SAM
        samtools calmd "$bam_file" "$ref_fasta" > "$calmd_sam" 2>/dev/null
        
        if [[ -s "$calmd_sam" ]]; then
            work_sam="$calmd_sam"
            log_info "Using calmd-processed SAM: $calmd_sam"
        else
            log_warning "samtools calmd failed or empty. Using native tags."
            rm -f "$calmd_sam" 2>/dev/null
        fi
    else
        log_info "Reference FASTA not found. Relying on aligner's native MD tags."
    fi

    # Step 3: Run parseAlignment.pl
    # Options:
    #   -v: verbose
    #   --map-qual 1: Require unique mapping (MAPQ >= 1)
    #   --min-len 16: Minimum read length
    #   --mutation-file: Output mutation file for CIMS
    log_info "Running parseAlignment.pl..."
    local cmd="parseAlignment.pl -v --map-qual 1 --min-len 16 --mutation-file '$mutation_file' '$work_sam' '$output_bed'"
    
    execute_cmd "$cmd"
    local parse_exit=$?
    
    # Step 4: Cleanup temp SAM files
    rm -f "$sam_file" 2>/dev/null
    if [[ -f "${sam_file%.sam}_calmd.sam" ]]; then
        rm -f "${sam_file%.sam}_calmd.sam" 2>/dev/null
    fi
    
    if [[ $parse_exit -ne 0 ]]; then
        log_error "parseAlignment.pl failed with exit code $parse_exit"
        return $parse_exit
    fi
    
    # Verify outputs
    if [[ ! -s "$output_bed" ]]; then
        log_warning "parseAlignment.pl produced empty BED file: $output_bed"
    else
        local bed_count=$(wc -l < "$output_bed")
        log_info "parseAlignment.pl output: $bed_count tags in $output_bed"
    fi
    
    if [[ ! -s "$mutation_file" ]]; then
        log_warning "parseAlignment.pl produced empty mutation file: $mutation_file"
    else
        local mut_count=$(wc -l < "$mutation_file")
        log_info "parseAlignment.pl output: $mut_count mutations in $mutation_file"
    fi
    
    log_info "Alignment parsing complete."
}

# 2b. Mapping with Bowtie2
run_mapping_bowtie2() {
    local input_file="$1"
    local output_prefix="$2"
    local genome_index="$3"
    local threads="$4"
    
    update_status "Mapping (Bowtie2)"
    log_info "Starting mapping with Bowtie2..."
    log_info "Input: $input_file"
    log_info "Index: $genome_index"
    
    # Verify index (look for .1.bt2 or .1.bt2l, excluding ncRNA patterns)
    # First try top-level (maxdepth 1), then fall back to deeper search
    local found_idx=""
    
    # Try .1.bt2 at top level first, excluding ncRNA
    found_idx=$(find "$genome_index" -maxdepth 1 -name "*.1.bt2" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    
    # If not found, try .1.bt2l at top level
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -maxdepth 1 -name "*.1.bt2l" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    
    # Fall back to deeper search (excluding ncRNA subfolder)
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -name "*.1.bt2" ! -path "*/ncRNA/*" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    if [[ -z "$found_idx" ]]; then
        found_idx=$(find "$genome_index" -name "*.1.bt2l" ! -path "*/ncRNA/*" ! -name "*ncrna*" ! -name "*.rev.*" 2>/dev/null | head -n 1)
    fi
    
    if [[ -z "$found_idx" ]]; then
        log_error "Bowtie2 index files (*.1.bt2) not found in $genome_index"
        return 1
    fi
    
    # Construct base name for index
    # Standard: /path/hg38.1.bt2 -> Base: /path/hg38
    # We strip the .1.bt2 suffix
    local idx_base="${found_idx%.1.bt2}"
    if [[ "$idx_base" == "$found_idx" ]]; then
         idx_base="${found_idx%.1.bt2l}"
    fi

    local sam_file="${output_prefix}.sam"
    local bam_file="${output_prefix}.Aligned.sortedByCoord.out.bam" # Match STAR naming for compatibility
    
    local cmd="bowtie2 -p $threads --sam-opt-config 'md' -x '$idx_base' -U '$input_file' -S '$sam_file' $ADV_ALIGNER_ARGS"
        
    log_info "Running Bowtie2..."
    execute_cmd "$cmd"
    
    if [ $? -ne 0 ]; then
        log_error "Bowtie2 alignment failed."
        return 1
    fi
    
    # Convert to BAM -> Sort -> Index
    # Note: "Processing Alignment" status is now in run_parse_alignment
    
    local sort_cmd="samtools view -bS '$sam_file' | samtools sort -@ $threads -o '$bam_file' -"
    execute_cmd "$sort_cmd"
    
    if [ -f "$bam_file" ]; then
        samtools index "$bam_file"
        rm -f "$sam_file" # Cleanup SAM
    else
        log_error "BAM conversion failed."
        return 1
    fi
    log_info "Mapping complete. Output: $bam_file"
}

# 3. PCR Duplicate Removal (using -c for chromosome-based processing to prevent OOM)
run_collapse_pcr() {
    local input_bed="$1"
    local output_bed="$2"
    local umi_len="${3:-0}"  # Optional: UMI length, defaults to 0

    update_status "Collapsing"
    log_info "Collapsing PCR duplicates with CTK tag2collapse.pl..."
    
    # Only use --random-barcode when UMI is present in read names
    local barcode_flag=""
    local weight_flags=""
    if [ "${umi_len:-0}" -gt 0 ]; then
        barcode_flag="--random-barcode"
        # CTK weight flags: use pre-collapse counts from fastq2collapse.pl
        # --weight: use tag count as weight
        # --weight-in-name: read count from #count# in read ID
        # -EM 30: EM iterations for barcode collapse
        # --seq-error-model alignment: estimate error from alignment
        weight_flags="--weight --weight-in-name -EM 30 --seq-error-model alignment"
        log_info "UMI mode (length=$umi_len): using random barcode collapse with EM weighting"
    else
        log_info "No UMI: using position-only collapse"
    fi
    
    # Create temp cache directory path for -c option (tag2collapse.pl creates it)
    local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/collapse_cache.XXXXXX")
    
    # Use -big (memory-mapped BIG format) and -c (chromosome-based) for max memory efficiency
    local cmd="$CONDA_PREFIX/bin/perl $(which tag2collapse.pl) -big -c \"${cache_dir}\" --keep-tag-name --keep-max-score ${barcode_flag} ${weight_flags} \
        \"${input_bed}\" \"${output_bed}\""

    execute_cmd "$cmd"
    local exit_code=$?
    
    # Cleanup cache directory
    rm -rf "$cache_dir"

    if [ $exit_code -eq 0 ] && [[ -s "$output_bed" ]]; then
        local output_count=$(wc -l < "$output_bed")
        log_info "Collapsing complete. Output: $output_count tags"
    else
        log_error "PCR duplicate removal failed."
        exit 1
    fi
}

# 4. Peak Calling (HOMER)
run_peak_calling() {
    local input_bed="$1"
    local out_dir="$2"
    local peak_dist="$3"
    local peak_size="$4"
    local frag_len="$5"
    local log_file="${out_dir}_homer.log"

    # Only called in aggregation or single sample (non-batch)
    update_status "Peaks"
    log_info "Calling peaks with HOMER..."
    
    # Capture output to specific log file for extraction later
    echo "Running HOMER makeTagDirectory..." > "$log_file"
    makeTagDirectory "${out_dir}" "${input_bed}" -single -format bed >> "$log_file" 2>&1
    
    echo "Running HOMER findPeaks..." >> "$log_file"
    findPeaks "${out_dir}" -o auto -style factor -L 2 -localSize 10000 -strand separate \
        -minDist "${peak_dist}" -size "${peak_size}" -fragLength "${frag_len}" $ADV_HOMER_ARGS >> "$log_file" 2>&1
        
    log_info "Peak calling complete. Log saved to $log_file"
}

# 4b. Add Enhanced Columns to Peak Coverage Matrix
# Adds: BC (groups), Raw Group Counts, Normalized Counts, BedGraph Stats
# Column Order: BC -> Raw Counts -> Normalized Counts -> BG Stats
add_matrix_columns() {
    local peak_matrix="$1"      # Path to peakCoverage.txt
    local peaks_bed="$2"        # Path to peaks_Sorted.bed
    local bg_dir="$3"           # BedGraph directory
    local scale_file="$4"       # Scale factors TSV
    local groups_file="$5"      # Optional groups file
    
    log_info "Adding enhanced columns to peak matrix..."
    
    # Validate inputs
    if [[ ! -f "$peak_matrix" ]]; then
        log_error "Peak matrix not found: $peak_matrix"
        return 1
    fi
    
    # Extract sample names from scale_factors.tsv (reliable source of actual samples)
    # NOT from matrix header (which may include CTK columns)
    if [[ ! -f "$scale_file" ]]; then
        log_warning "Scale factors file not found: $scale_file. Enhanced columns will be limited."
        local samples=()
    else
        local samples=($(cut -f1 "$scale_file"))
    fi
    
    log_info "  Samples detected: ${samples[*]}"
    
    # Prepare temp files for new columns
    local new_cols_file=$(mktemp)
    local new_header=""
    
    # -------------------------------------------
    # STEP 1: Biological Complexity (BC) - Groups Only
    # -------------------------------------------
    if [[ -n "$groups_file" && -f "$groups_file" ]]; then
        log_info "  Calculating Biological Complexity (BC)..."
        local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
        
        for group in $unique_groups; do
            log_info "    Group: $group"
            # Get sample column indices for this group
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file" | tr '\n' ' ')
            
            # For each peak (row), count samples with count > 0
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                # Build map of sample name -> column index (1-based, offset by 6)
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
            }
            NR==1 { print "BC_'"$group"'"; next }
            {
                bc=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map) {
                        c = col_map[s]
                        if($c + 0 > 0) bc++
                    }
                }
                print bc
            }
            ' "$peak_matrix" > "bc_${group}.col"
            
            # Add to new columns
            if [[ ! -s "$new_cols_file" ]]; then
                cat "bc_${group}.col" > "$new_cols_file"
            else
                paste "$new_cols_file" "bc_${group}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
            fi
            rm -f "bc_${group}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 2: Normalized Read Counts (Per Sample)
    # -------------------------------------------
    if [[ -f "$scale_file" ]]; then
        log_info "  Adding normalized read counts..."
        
        for sample in "${samples[@]}"; do
            # Get scale factor for this sample using awk for reliable tab-delimited matching
            local sf=$(awk -F'\t' -v s="$sample" '$1==s {print $3; exit}' "$scale_file")
            if [[ -z "$sf" ]]; then
                log_warning "    Scale factor not found for $sample, using 1.0"
                sf="1.0"
            fi
            
            # Get column index for this sample (1-based)
            local col_idx=0
            for i in "${!samples[@]}"; do
                if [[ "${samples[$i]}" == "$sample" ]]; then
                    col_idx=$((i + 7))  # Offset by 6 base columns + 1 for 1-indexing
                    break
                fi
            done
            
            # Calculate normalized value: raw_count * scale_factor
            awk -F'\t' -v col="$col_idx" -v sf="$sf" '
            NR==1 { print "NormedTC_'"${sample}"'"; next }
            { printf "%.4f\n", $col * sf }
            ' "$peak_matrix" > "normed_${sample}.col"
            
            # Add to new columns
            if [[ ! -s "$new_cols_file" ]]; then
                cat "normed_${sample}.col" > "$new_cols_file"
            else
                paste "$new_cols_file" "normed_${sample}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
            fi
            rm -f "normed_${sample}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 3: Group Columns (Raw Sum + Normalized Sum) - Groups Only
    # -------------------------------------------
    if [[ -n "$groups_file" && -f "$groups_file" ]]; then
        log_info "  Adding group aggregate columns..."
        local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
        
        for group in $unique_groups; do
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file" | tr '\n' ' ')
            
            # Sum raw counts for group
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
            }
            NR==1 { print "TC_'"$group"'"; next }
            {
                sum=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map) sum += $col_map[s]
                }
                print sum
            }
            ' "$peak_matrix" > "grp_raw_${group}.col"
            
            paste "$new_cols_file" "grp_raw_${group}.col" > "${new_cols_file}.tmp"
            mv "${new_cols_file}.tmp" "$new_cols_file"
            rm -f "grp_raw_${group}.col"
            
            # Sum normalized counts for group
            # This needs to reference the normalized columns we just added
            # For simplicity, recalculate using scale factors
            awk -F'\t' -v samples="$group_samples" -v allsamples="${samples[*]}" -v sf_file="$scale_file" '
            BEGIN {
                split(samples, gs, " ")
                split(allsamples, as, " ")
                for(i=1; i<=length(as); i++) col_map[as[i]] = i + 6
                # Load scale factors
                while((getline line < sf_file) > 0) {
                    split(line, sf_parts, "\t")
                    # Extract sample name from path
                    n = split(sf_parts[1], path_parts, "/")
                    sname = path_parts[n]
                    scale[sname] = sf_parts[3]
                }
            }
            NR==1 { print "NormedTC_'"$group"'"; next }
            {
                sum=0
                for(i in gs) {
                    s = gs[i]
                    if(s in col_map && s in scale) {
                        sum += $col_map[s] * scale[s]
                    }
                }
                printf "%.4f\n", sum
            }
            ' "$peak_matrix" > "grp_normed_${group}.col"
            
            paste "$new_cols_file" "grp_normed_${group}.col" > "${new_cols_file}.tmp"
            mv "${new_cols_file}.tmp" "$new_cols_file"
            rm -f "grp_normed_${group}.col"
        done
    fi
    
    # -------------------------------------------
    # STEP 4: BedGraph Stats (Sum/Avg/Max) Per Sample
    # -------------------------------------------
    if [[ -d "$bg_dir" ]] && [[ ${#samples[@]} -gt 0 ]]; then
        log_info "  Adding BedGraph statistics..."
        
        # Split peaks by strand and sort for bedtools compatibility (same as STEP 5)
        awk -F'\t' '$6=="+"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_pos.tmp.bed
        awk -F'\t' '$6=="-"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_neg.tmp.bed
        
        for sample in "${samples[@]}"; do
            local bg_pos="${bg_dir}/${sample}_pos.bedgraph"
            local bg_neg="${bg_dir}/${sample}_neg.bedgraph"
            
            if [[ ! -f "$bg_pos" || ! -f "$bg_neg" ]]; then
                log_warning "    BedGraph not found for $sample, skipping."
                continue
            fi
            
            # Sort bedgraphs for bedtools compatibility (same as STEP 5)
            sort -k1,1 -k2,2n "$bg_pos" > "${sample}_pos_sorted.bg.tmp"
            sort -k1,1 -k2,2n "$bg_neg" > "${sample}_neg_sorted.bg.tmp"
            
            for stat in sum mean max; do
                # Run bedtools map on sorted files (same as STEP 5)
                bedtools map -a peaks_pos.tmp.bed -b "${sample}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_pos_${stat}.tmp" 2>/dev/null
                bedtools map -a peaks_neg.tmp.bed -b "${sample}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_neg_${stat}.tmp" 2>/dev/null
                
                # Build lookup and match to original peak order (same logic as STEP 5)
                awk -F'\t' '
                BEGIN { 
                    while((getline < "bg_pos_'"$stat"'.tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                    while((getline < "bg_neg_'"$stat"'.tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                }
                NR==1 { print "Cov" toupper(substr("'"$stat"'",1,1)) substr("'"$stat"'",2) "_'"$sample"'"; next }
                {
                    key = $1"\t"$2"\t"$3
                    if($6 == "+") print (key in pos) ? pos[key] : 0
                    else print (key in neg) ? neg[key] : 0
                }
                ' "$peaks_bed" > "bg_${sample}_${stat}.col"
                
                paste "$new_cols_file" "bg_${sample}_${stat}.col" > "${new_cols_file}.tmp"
                mv "${new_cols_file}.tmp" "$new_cols_file"
                rm -f "bg_${sample}_${stat}.col"
            done
            
            rm -f "${sample}_pos_sorted.bg.tmp" "${sample}_neg_sorted.bg.tmp"
        done
        
        rm -f peaks_pos.tmp.bed peaks_neg.tmp.bed "bg_pos_"*.tmp "bg_neg_"*.tmp
        
        # -------------------------------------------
        # STEP 5: Group BedGraph Stats (from combined bedgraph) - Groups Only
        # -------------------------------------------
        if [[ -n "$groups_file" && -f "$groups_file" ]]; then
            log_info "  Adding group BedGraph statistics..."
            local combined_bg_dir="${bg_dir}/COMBINED_BEDGRAPH"
            local unique_groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort -u)
            
            # Re-split and sort peaks for bedtools compatibility
            awk -F'\t' '$6=="+"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_pos.tmp.bed
            awk -F'\t' '$6=="-"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_neg.tmp.bed
            
            for group in $unique_groups; do
                local grp_bg_pos="${combined_bg_dir}/${group}_combined_pos.bedgraph"
                local grp_bg_neg="${combined_bg_dir}/${group}_combined_neg.bedgraph"
                
                if [[ ! -f "$grp_bg_pos" || ! -f "$grp_bg_neg" ]]; then
                    log_warning "    Combined BedGraph not found for $group, skipping."
                    continue
                fi
                
                # Sort group bedgraphs for bedtools compatibility (and ensure TABs)
                tr ' ' '\t' < "$grp_bg_pos" | sort -k1,1 -k2,2n > "${group}_pos_sorted.bg.tmp"
                tr ' ' '\t' < "$grp_bg_neg" | sort -k1,1 -k2,2n > "${group}_neg_sorted.bg.tmp"
                
                for stat in sum mean max; do
                    bedtools map -a peaks_pos.tmp.bed -b "${group}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_pos_${stat}.tmp"
                    bedtools map -a peaks_neg.tmp.bed -b "${group}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 > "bg_neg_${stat}.tmp"
                    
                    awk -F'\t' '
                    BEGIN { 
                        while((getline < "bg_pos_'"$stat"'.tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                        while((getline < "bg_neg_'"$stat"'.tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                    }
                    NR==1 { print "Cov" toupper(substr("'$stat'",1,1)) substr("'$stat'",2) "_'$group'"; next }
                    {
                        key = $1"\t"$2"\t"$3
                        if($6 == "+") print (key in pos) ? pos[key] : 0
                        else print (key in neg) ? neg[key] : 0
                    }
                    ' "$peaks_bed" > "bg_${group}_${stat}.col"
                    
                    paste "$new_cols_file" "bg_${group}_${stat}.col" > "${new_cols_file}.tmp"
                    mv "${new_cols_file}.tmp" "$new_cols_file"
                    rm -f "bg_${group}_${stat}.col"
                done
                
                rm -f "${group}_pos_sorted.bg.tmp" "${group}_neg_sorted.bg.tmp"
            done
            
            rm -f peaks_pos.tmp.bed peaks_neg.tmp.bed "bg_pos_"*.tmp "bg_neg_"*.tmp
        fi
    fi
    
    # -------------------------------------------
    # FINAL: Paste all new columns to matrix
    # -------------------------------------------
    if [[ -s "$new_cols_file" ]]; then
        paste "$peak_matrix" "$new_cols_file" > "${peak_matrix}.enhanced"
        mv "${peak_matrix}.enhanced" "$peak_matrix"
        log_info "Enhanced columns added to $peak_matrix"
    fi
    
    # -------------------------------------------
    # REORDER: Group columns by type (prefix)
    # Order: base -> TC_ -> NormedTC_ -> BC_ -> DEL_ -> SUB_ -> TRUNC_ -> BG*
    # -------------------------------------------
    log_info "Reordering columns by type..."
    awk -F'\t' '
    BEGIN { OFS="\t" }
    NR==1 {
        # Parse header and categorize columns by prefix
        for(i=1; i<=NF; i++) {
            h = $i
            if(h ~ /^(chr|start|end|name|score|strand)$/) { order[i] = 1; base[i] = h }
            else if(h ~ /^TC_/)      { order[i] = 2; tc[i] = h }
            else if(h ~ /^NormedTC_/) { order[i] = 3; ntc[i] = h }
            else if(h ~ /^BC_/)      { order[i] = 4; bc[i] = h }
            else if(h ~ /^DEL_/)     { order[i] = 5; del[i] = h }
            else if(h ~ /^SUB_/)     { order[i] = 6; sub_[i] = h }
            else if(h ~ /^TRUNC_/)   { order[i] = 7; trunc[i] = h }
            else if(h ~ /^BG/)       { order[i] = 8; bg[i] = h }
            else                     { order[i] = 9; other[i] = h }
            headers[i] = h
        }
        # Build column order
        n = 0
        for(o=1; o<=9; o++) {
            for(i=1; i<=NF; i++) {
                if(order[i] == o) { col_order[++n] = i }
            }
        }
        total_cols = n
        # Print reordered header
        for(j=1; j<=total_cols; j++) {
            printf "%s%s", headers[col_order[j]], (j<total_cols ? OFS : ORS)
        }
        next
    }
    {
        # Print reordered data
        for(j=1; j<=total_cols; j++) {
            printf "%s%s", $col_order[j], (j<total_cols ? OFS : ORS)
        }
    }
    ' "$peak_matrix" > "${peak_matrix}.reordered"
    mv "${peak_matrix}.reordered" "$peak_matrix"
    log_info "Columns reordered by type prefix"
    
    rm -f "$new_cols_file"
}

# 4b. Add CTK columns to peak coverage matrix
# Adds _del, _sub, _trunc columns for each sample/group with CTK outputs
add_ctk_columns_to_peak_matrix() {
    local peak_matrix="$1"         # Input/output: peak coverage matrix file
    local peaks_bed="$2"           # Sorted peaks BED file for bedtools
    local ctk_dir="$3"             # Directory containing CTK outputs
    local cims_fdr="${4:-0.05}"    # FDR threshold for CIMS filtering

    local cits_pval="${5:-0.05}"   # P-value threshold for CITS filtering
    local groups_file="$6"         # Optional: Groups file for aggregation
    
    log_info "Adding CTK site counts to peak matrix..."
    log_info "CTK directory: $ctk_dir"
    

    
    # Helper function to add column from CTK file
    add_ctk_column() {
        local ctk_file="$1"
        local column_name="$2"
        local ctk_type="$3"  # "cims" or "cits"
        local threshold="$4"
        
        if [[ ! -s "$ctk_file" ]]; then
            log_warning "CTK file not found or empty: $ctk_file"
            return 1
        fi
        
        local temp_filtered="${ctk_file}.filtered.bed"
        local temp_coverage="${ctk_file}.coverage.txt"
        
        # Filter by threshold and convert to BED
        if [[ "$ctk_type" == "cims" ]]; then
            # CIMS: Column 9 is FDR, skip header line
            grep -v "^#" "$ctk_file" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$temp_filtered"
        else
            # CITS: P-value is embedded in name as [P=value]
            grep -v "^#" "$ctk_file" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$temp_filtered"
        fi
        
        local site_count=$(wc -l < "$temp_filtered")
        log_info "  $column_name: $site_count significant sites"
        
        if [[ "$site_count" -eq 0 ]]; then
            # Add column of zeros
            local num_rows=$(($(wc -l < "$peak_matrix") - 1))  # Subtract header
            echo "$column_name" > temp_col.txt
            yes "0" | head -n "$num_rows" >> temp_col.txt
        else
            # Count sites per peak using bedtools
            bedtools coverage -s -a "$peaks_bed" -b "$temp_filtered" > "$temp_coverage"
            
            # Extract count column (7th) and combine with header
            echo "$column_name" > temp_col.txt
            awk '{print $7}' "$temp_coverage" >> temp_col.txt
        fi
        
        # Paste to matrix
        paste "$peak_matrix" temp_col.txt > temp_matrix.txt
        mv temp_matrix.txt "$peak_matrix"
        
        # Cleanup
        rm -f "$temp_filtered" "$temp_coverage" temp_col.txt
    }
    

    
    # -------------------------------------------------------------------------
    # Scenario A: Group Aggregation Mode (if groups_file provided)
    # -------------------------------------------------------------------------
    if [[ -n "$groups_file" ]]; then
        log_info "Using Groups File for CTK Aggregation..."
        
        # Parse groups
        local groups_map=$(mktemp)
        parse_groups_file "$groups_file" "$groups_map"
        local unique_groups=$(cut -f2 "$groups_map" | sort -u)
        
        for group in $unique_groups; do
            # Get samples for this group
            local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
            
            # --- Aggregate CIMS (Deletions) ---
            local group_del_bed="${ctk_dir}/${group}_aggregated_CIMS_del.txt"
            local precalc_del="${ctk_dir}/${group}/CIMS/${group}_CIMS_del.txt"
            if [[ -s "$precalc_del" ]]; then
                cp "$precalc_del" "$group_del_bed"
            else
                > "$group_del_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    # Look for sample CIMS file in CTK dir or subdirs
                    local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_del.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_del_bed"
                done
            fi
            if [[ -s "$group_del_bed" ]]; then
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            else
                # still add empty column for consistency?
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CIMS (Substitutions) ---
            local group_sub_bed="${ctk_dir}/${group}_aggregated_CIMS_sub.txt"
            local precalc_sub="${ctk_dir}/${group}/CIMS/${group}_CIMS_sub.txt"
             if [[ -s "$precalc_sub" ]]; then
                cp "$precalc_sub" "$group_sub_bed"
            else
                > "$group_sub_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_sub.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_sub_bed"
                done
            fi
            if [[ -s "$group_sub_bed" ]]; then
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            else
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CITS (Truncations) ---
            local group_cits_bed="${ctk_dir}/${group}_aggregated_CITS.txt"
            local precalc_cits="${ctk_dir}/${group}/CITS/${group}_CITS.bed"
            # Try .bed first, then .txt
            if [[ ! -s "$precalc_cits" ]]; then precalc_cits="${ctk_dir}/${group}/CITS/${group}_CITS.txt"; fi
            
             if [[ -s "$precalc_cits" ]]; then
                cp "$precalc_cits" "$group_cits_bed"
            else
                > "$group_cits_bed"
                local samples=$(awk -F'\t' -v g="$group" '$2==g {print $1}' "$groups_map" | tr '\n' ' ')
                for sample in $samples; do
                    local s_file=$(find "$ctk_dir" -name "${sample}_CITS.txt" 2>/dev/null | head -n 1)
                    [[ -s "$s_file" ]] && cat "$s_file" >> "$group_cits_bed"
                done
            fi
            if [[ -s "$group_cits_bed" ]]; then
                add_ctk_column "$group_cits_bed" "TRUNC_${group}" "cits" "$cits_pval"
            else
                add_ctk_column "$group_cits_bed" "TRUNC_${group}" "cits" "$cits_pval"
            fi
            
            # Cleanup aggregated files
            rm -f "$group_del_bed" "$group_sub_bed" "$group_cits_bed"
        done
        
        rm "$groups_map"
        return 0
    fi
    
    # -------------------------------------------------------------------------
    # Scenario B: Default / Legacy Mode (All found in dir)
    # -------------------------------------------------------------------------
    # Find and process CIMS deletion files
    # Check both: 1) Direct structure: ctk_dir/CIMS/  2) Group structure: ctk_dir/*/CIMS/
    local cims_found=false
    
    # First try direct structure
    if [[ -d "${ctk_dir}/CIMS" ]]; then
        cims_found=true
        for cims_del_file in "${ctk_dir}/CIMS/"*_CIMS_del.txt; do
            if [[ -f "$cims_del_file" ]]; then
                local name=$(basename "$cims_del_file" _CIMS_del.txt)
                add_ctk_column "$cims_del_file" "DEL_${name}" "cims" "$cims_fdr"
            fi
        done
        
        for cims_sub_file in "${ctk_dir}/CIMS/"*_CIMS_sub.txt; do
            if [[ -f "$cims_sub_file" ]]; then
                local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                add_ctk_column "$cims_sub_file" "SUB_${name}" "cims" "$cims_fdr"
            fi
        done
    fi
    
    # If no direct CIMS folder, try group subfolders (ctk_dir/*/CIMS/)
    if [[ "$cims_found" == "false" ]]; then
        for group_dir in "${ctk_dir}"/*/; do
            if [[ -d "${group_dir}CIMS" ]]; then
                cims_found=true
                for cims_del_file in "${group_dir}CIMS/"*_CIMS_del.txt; do
                    if [[ -f "$cims_del_file" ]]; then
                        local name=$(basename "$cims_del_file" _CIMS_del.txt)
                        add_ctk_column "$cims_del_file" "DEL_${name}" "cims" "$cims_fdr"
                    fi
                done
                
                for cims_sub_file in "${group_dir}CIMS/"*_CIMS_sub.txt; do
                    if [[ -f "$cims_sub_file" ]]; then
                        local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                        add_ctk_column "$cims_sub_file" "SUB_${name}" "cims" "$cims_fdr"
                    fi
                done
            fi
        done
    fi
    
    # Find and process CITS files
    # Check both: 1) Direct structure: ctk_dir/CITS/  2) Group structure: ctk_dir/*/CITS/
    local cits_found=false
    
    # First try direct structure
    if [[ -d "${ctk_dir}/CITS" ]]; then
        cits_found=true
        for cits_file in "${ctk_dir}/CITS/"*_CITS.txt "${ctk_dir}/CITS/"*_CITS.bed; do
            if [[ -f "$cits_file" ]]; then
                local name=$(basename "$cits_file" | sed 's/_CITS\.\(txt\|bed\)$//')
                add_ctk_column "$cits_file" "TRUNC_${name}" "cits" "$cits_pval"
            fi
        done
    fi
    
    # If no direct CITS folder, try group subfolders
    if [[ "$cits_found" == "false" ]]; then
        for group_dir in "${ctk_dir}"/*/; do
            if [[ -d "${group_dir}CITS" ]]; then
                cits_found=true
                for cits_file in "${group_dir}CITS/"*_CITS.txt "${group_dir}CITS/"*_CITS.bed; do
                    if [[ -f "$cits_file" ]]; then
                        local name=$(basename "$cits_file" | sed 's/_CITS\.\(txt\|bed\)$//')
                        add_ctk_column "$cits_file" "TRUNC_${name}" "cits" "$cits_pval"
                    fi
                done
            fi
        done
    fi
    
    log_info "CTK columns added to peak matrix."
}

# 5. CIMS Analysis (CTK)
# Detects crosslinking-induced mutation sites
# ═══════════════════════════════════════════════════════════════════════════
# CTK CIMS/CITS ANALYSIS FUNCTIONS
# Verified workflow based on CTK documentation (Dec 2024)
# ═══════════════════════════════════════════════════════════════════════════

# CTK Preprocessing: Filter mutations and extract mutation types
# Called after tag2collapse.pl, before CIMS/CITS
run_ctk_preprocessing() {
    local collapsed_bed="$1"
    local raw_mutation_file="$2"
    local output_dir="$3"
    
    # Status update removed - preprocessing is already done in main pipeline
    log_info "Preprocessing mutations for CTK analysis..."
    
    mkdir -p "$output_dir"
    
    # Step 1: Filter mutations to only those in collapsed tags
    # selectRow.pl uses zero-based column indexing: column 3 = read name
    log_info "Filtering mutations to collapsed tags (selectRow.pl -q 3 -f 3)..."
    local matched_file="${output_dir}/mutations_matched.txt"
    
    selectRow.pl -q 3 -f 3 "$raw_mutation_file" "$collapsed_bed" > "$matched_file"
    
    if [[ ! -s "$matched_file" ]]; then
        log_warning "No matching mutations found after filtering."
        return 1
    fi
    
    local matched_count=$(wc -l < "$matched_file")
    log_info "Matched mutations: $matched_count"
    
    # Step 2: Extract mutation types using getMutationType.pl
    log_info "Extracting deletion mutations..."
    local del_file="${output_dir}/deletions.bed"
    getMutationType.pl -t del "$matched_file" "$del_file"
    
    log_info "Extracting substitution mutations..."
    local sub_file="${output_dir}/substitutions.bed"
    getMutationType.pl -t sub "$matched_file" "$sub_file"
    
    # Report counts
    if [[ -s "$del_file" ]]; then
        local del_count=$(wc -l < "$del_file")
        log_info "Deletions extracted: $del_count"
    else
        log_warning "No deletions found."
    fi
    
    if [[ -s "$sub_file" ]]; then
        local sub_count=$(wc -l < "$sub_file")
        log_info "Substitutions extracted: $sub_count"
    else
        log_warning "No substitutions found."
    fi
    
    log_info "CTK preprocessing complete. Output: $output_dir"
}

# 5. CIMS Analysis (CTK) - with parallel chromosome processing
# Detects crosslinking-induced mutation sites
# Input: collapsed BED + mutation file (deletions or substitutions)
# Output: CIMS.txt with significant mutation sites
run_cims() {
    local input_collapsed_bed="$1"
    local mutation_bed="$2"          # Already BED6 from getMutationType.pl
    local output_file="$3"
    local cims_iterations="${4:-10}"  # Default: 10 iterations
    local cims_fdr="${5:-0.001}"      # Default: FDR 0.001
    local threads="${THREADS:-4}"     # Use global THREADS or default
    
    log_info "Running CIMS analysis..."
    log_info "Input tags: $input_collapsed_bed"
    log_info "Input mutations: $mutation_bed"
    log_info "Iterations: $cims_iterations, FDR threshold: $cims_fdr"
    
    # Verify inputs exist
    if [[ ! -s "$input_collapsed_bed" ]]; then
        log_error "CIMS: Collapsed BED file is empty or missing: $input_collapsed_bed"
        return 1
    fi
    if [[ ! -s "$mutation_bed" ]]; then
        log_error "CIMS: Mutation BED file is empty or missing: $mutation_bed"
        return 1
    fi
    
    # Set up CTK environment
    export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
    
    # Check if parallel processing is available and beneficial
    local use_parallel="false"
    local chr_count=$(cut -f1 "$input_collapsed_bed" | sort -u | wc -l)
    
    if has_gnu_parallel && [[ $chr_count -gt 1 ]]; then
        use_parallel="true"
        log_info "Using parallel processing ($chr_count chromosomes, $threads threads)"
    fi
    
    if [[ "$use_parallel" == "true" ]]; then
        # Parallel mode: split BOTH collapsed.bed AND mutation file by chromosome
        local chunk_dir=$(mktemp -d "${TMPDIR:-/tmp}/cims_parallel.XXXXXX")
        
        # Split collapsed BED by chromosome (standard chromosomes only: chr1-22, X, Y, M)
        # This filters out contigs like KI270737.1, GL000220.1, etc.
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$input_collapsed_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".bed")}'
        
        # CRITICAL FIX: Also split mutation file by chromosome to match collapsed.bed chunks
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$mutation_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".mut.bed")}'
        
        # Create processing script (updated to use per-chromosome mutation file)
        local process_script="$chunk_dir/run_cims_chunk.sh"
        cat > "$process_script" << 'CIMS_SCRIPT'
#!/bin/bash
chunk_file="$1"
iterations="$2"
output_dir="$3"
chunk_name=$(basename "$chunk_file" .bed)
# Use matching mutation file for this chromosome
mutation_chunk="${output_dir}/${chunk_name}.mut.bed"
output_chunk="${output_dir}/${chunk_name}.cims.txt"
cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cims_cache.XXXXXX")
export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
# Only run CIMS if both chunk files exist and are non-empty
if [[ -s "$chunk_file" && -s "$mutation_chunk" ]]; then
    CIMS.pl -big -c "$cache_dir" -n "$iterations" "$chunk_file" "$mutation_chunk" "$output_chunk" >/dev/null 2>&1
fi
rm -rf "$cache_dir" 2>/dev/null
CIMS_SCRIPT
        chmod +x "$process_script"
        
        # Calculate optimal parallel jobs based on RAM and file size
        local optimal_jobs=$(calculate_optimal_parallel_jobs "$threads" "$input_collapsed_bed")
        local parallel_jobs=$(( optimal_jobs < chr_count ? optimal_jobs : chr_count ))
        log_info "CIMS: Using $parallel_jobs/$threads threads based on available RAM"
        
        # Pass only chunk_file, iterations, and chunk_dir (mutation file is derived from chunk name)
        ls "$chunk_dir"/chr_*.bed 2>/dev/null | grep -v '\.mut\.bed$' | \
            parallel --memfree 4G -j "$parallel_jobs" \
            "$process_script" {} "$cims_iterations" "$chunk_dir"
        
        # Merge results (skip empty files, keep header from first)
        local first_file="true"
        > "$output_file"
        for result in "$chunk_dir"/*.cims.txt; do
            if [[ -s "$result" ]]; then
                if [[ "$first_file" == "true" ]]; then
                    cat "$result" >> "$output_file"
                    first_file="false"
                else
                    # Skip header line for subsequent files
                    tail -n +2 "$result" >> "$output_file"
                fi
            fi
        done
        
        # Cleanup
        rm -rf "$chunk_dir"
    else
        # Sequential mode (fallback)
        if [[ "$use_parallel" == "false" ]] && has_gnu_parallel; then
            log_info "Single chromosome detected, using sequential processing"
        elif ! has_gnu_parallel; then
            log_info "GNU parallel not available, using sequential processing"
        fi
        # Use -big -c for memory efficiency in sequential mode
        local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cims_cache.XXXXXX")
        local cmd="CIMS.pl -big -c '$cache_dir' -v -n $cims_iterations '$input_collapsed_bed' '$mutation_bed' '$output_file'"
        execute_cmd "$cmd"
        rm -rf "$cache_dir" 2>/dev/null
    fi
    
    # Process results
    if [[ -s "$output_file" ]]; then
        local raw_count=$(grep -v "^#" "$output_file" | wc -l)
        log_info "CIMS complete: $raw_count total sites in $output_file"
        
        # Filter for significance only if FDR < 1
        if (( $(echo "$cims_fdr < 1" | bc -l) )); then
            log_info "Filtering CIMS results by FDR < $cims_fdr..."
            local temp_file="${output_file}.tmp"
            awk -F'\t' -v fdr="$cims_fdr" 'NR==1 || $9 < fdr' "$output_file" | \
                sort -k9,9n -k8,8nr -k7,7n > "$temp_file"
            mv "$temp_file" "$output_file"
            
            local filtered_count=$(grep -v "^#" "$output_file" | wc -l)
            log_info "CIMS filtered: $filtered_count sites (FDR < $cims_fdr)"
            
            # Cleanup unwanted CIMS intermediates
            rm -f mutations_matched.txt deletions.bed substitutions.bed 2>/dev/null
        fi
    else
        echo -e "[WARNING] CIMS produced empty output: $output_file" >> "${LOG_FILE}"
        echo -ne "${YELLOW}[WARNING] Empty Output${NC} > "
        return 1
    fi
    
    return 0
}

# 6. CITS Analysis (CTK) - with parallel chromosome processing
# Detects crosslinking-induced truncation sites
# Input: collapsed BED + deletion file (used to EXCLUDE read-through tags)
# Output: CITS.txt with significant truncation sites (single-nucleotide singletons)
run_cits() {
    local input_collapsed_bed="$1"
    local deletion_bed="$2"           # Used to exclude read-through tags
    local output_file="$3"            # Should be .txt extension
    local cits_pvalue="${4:-0.001}"   # Default: p-value 0.001
    local cits_gap="${5:-25}"         # Default: gap 25 for clustering
    local threads="${THREADS:-4}"     # Use global THREADS or default
    
    log_info "Running CITS analysis..."
    log_info "Input tags: $input_collapsed_bed"
    log_info "Deletion file (to exclude): $deletion_bed"
    log_info "P-value: $cits_pvalue, Gap: $cits_gap"
    
    # Verify inputs exist
    if [[ ! -s "$input_collapsed_bed" ]]; then
        log_error "CITS: Collapsed BED file is empty or missing: $input_collapsed_bed"
        return 1
    fi
    if [[ ! -s "$deletion_bed" ]]; then
        log_warning "CITS: Deletion file is empty. Running without read-through filtering."
    fi
    
    # Set up CTK environment
    export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
    
    # Check if parallel processing is available and beneficial
    local use_parallel="false"
    local chr_count=$(cut -f1 "$input_collapsed_bed" | sort -u | wc -l)
    
    if has_gnu_parallel && [[ $chr_count -gt 1 ]]; then
        use_parallel="true"
        log_info "Using parallel processing ($chr_count chromosomes, $threads threads)"
    fi
    
    local cits_raw="${output_file%.txt}_raw.bed"
    
    if [[ "$use_parallel" == "true" ]]; then
        # Parallel mode: split BOTH collapsed.bed AND deletion file by chromosome
        local chunk_dir=$(mktemp -d "${TMPDIR:-/tmp}/cits_parallel.XXXXXX")
        
        # Split collapsed BED by chromosome (standard chromosomes only: chr1-22, X, Y, M)
        # This filters out contigs like KI270737.1, GL000220.1, etc.
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$input_collapsed_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".bed")}'
        
        # CRITICAL FIX: Also split deletion file by chromosome to match collapsed.bed chunks
        grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' "$deletion_bed" | \
            awk -v dir="$chunk_dir" '{print > (dir"/chr_"$1".del.bed")}'
        
        # Create processing script (updated to use per-chromosome deletion file)
        local process_script="$chunk_dir/run_cits_chunk.sh"
        cat > "$process_script" << 'CITS_SCRIPT'
#!/bin/bash
chunk_file="$1"
pvalue="$2"
gap="$3"
output_dir="$4"
chunk_name=$(basename "$chunk_file" .bed)
# Use matching deletion file for this chromosome
deletion_chunk="${output_dir}/${chunk_name}.del.bed"
output_chunk="${output_dir}/${chunk_name}.cits.bed"
cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cits_cache.XXXXXX")
export PERL5LIB="${CONDA_PREFIX}/lib/czplib:$PERL5LIB"
# Only run CITS if collapsed chunk exists (deletion chunk can be empty)
if [[ -s "$chunk_file" ]]; then
    if [[ -s "$deletion_chunk" ]]; then
        CITS.pl -big -c "$cache_dir" -p "$pvalue" --gap "$gap" "$chunk_file" "$deletion_chunk" "$output_chunk" >/dev/null 2>&1
    else
        # No deletions for this chromosome - run without deletion filtering
        touch "$output_chunk"
    fi
fi
rm -rf "$cache_dir" 2>/dev/null
CITS_SCRIPT
        chmod +x "$process_script"
        
        # Calculate optimal parallel jobs based on RAM and file size
        local optimal_jobs=$(calculate_optimal_parallel_jobs "$threads" "$input_collapsed_bed")
        local parallel_jobs=$(( optimal_jobs < chr_count ? optimal_jobs : chr_count ))
        log_info "CITS: Using $parallel_jobs/$threads threads based on available RAM"
        
        # Pass only chunk_file, pvalue, gap, and chunk_dir (deletion file is derived from chunk name)
        ls "$chunk_dir"/chr_*.bed 2>/dev/null | grep -v '\.del\.bed$' | \
            parallel --memfree 4G -j "$parallel_jobs" \
            "$process_script" {} "$cits_pvalue" "$cits_gap" "$chunk_dir"
        
        # Merge results
        > "$cits_raw"
        for result in "$chunk_dir"/*.cits.bed; do
            if [[ -s "$result" ]]; then
                cat "$result" >> "$cits_raw"
            fi
        done
        
        # Cleanup
        rm -rf "$chunk_dir"
    else
        # Sequential mode (fallback)
        if [[ "$use_parallel" == "false" ]] && has_gnu_parallel; then
            log_info "Single chromosome detected, using sequential processing"
        elif ! has_gnu_parallel; then
            log_info "GNU parallel not available, using sequential processing"
        fi
        # Use -big -c for memory efficiency in sequential mode
        local cache_dir=$(mktemp -u "${TMPDIR:-/tmp}/cits_cache.XXXXXX")
        local cmd="CITS.pl -big -c '$cache_dir' -p $cits_pvalue --gap $cits_gap -v '$input_collapsed_bed' '$deletion_bed' '$cits_raw'"
        execute_cmd "$cmd"
        rm -rf "$cache_dir" 2>/dev/null
    fi
    
    # Process results
    if [[ -s "$cits_raw" ]]; then
        local raw_count=$(wc -l < "$cits_raw")
        log_info "CITS raw output: $raw_count sites"
        
        # Filter to single-nucleotide sites (singleton) - this is the main output
        # Add header for consistency with CIMS output
        echo -e "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand" > "$output_file"
        awk '{if($3-$2==1) {print $0}}' "$cits_raw" >> "$output_file"
        
        local singleton_count=$(($(wc -l < "$output_file") - 1))  # Subtract header
        log_info "CITS complete: $singleton_count singleton sites in $output_file"
        
        # Remove intermediate raw file
        rm -f "$cits_raw"
        rm -f mutations_matched.txt deletions.bed substitutions.bed 2>/dev/null
    else
        echo -e "[WARNING] CITS produced empty output" >> "${LOG_FILE}"
        echo -ne "${YELLOW}[WARNING] Empty Output${NC} > "
        return 1
    fi
    
    return 0
}

# 7. Flanked BED Generation for Motif Analysis
# Generates ±10nt flanked regions around CIMS/CITS sites
# Creates flanked BED file alongside the input file (same directory)
# Users can run their own motif analysis tools on these files
generate_flanked_bed() {
    local input_bed="$1"
    local flank_nt="${2:-10}"         # Default: ±10 nucleotides
    
    if [[ ! -s "$input_bed" ]]; then
        return 1
    fi
    
    # Create flanked file alongside input: sample_CIMS_sub.txt → sample_CIMS_sub_flanked.bed
    local flanked_bed="${input_bed%.txt}_flanked.bed"
    
    awk -v n="$flank_nt" 'BEGIN{OFS="\t"} {
        start = $2 - n
        if (start < 0) start = 0
        print $1, start, $3 + n, $4, $5, $6
    }' "$input_bed" > "$flanked_bed"
    
    local site_count=$(wc -l < "$flanked_bed")
    log_info "Generated flanked BED (±${flank_nt}nt, $site_count sites): $(basename "$flanked_bed")"
}

# 8. Full CTK Analysis Pipeline
# Orchestrates the complete CIMS/CITS workflow based on RUN_CIMS and RUN_CITS flags
run_ctk_full_analysis() {
    local bam_file="$1"
    local output_dir="$2"
    local genome_fasta="$3"
    local cims_iterations="${4:-10}"
    local cims_fdr="${5:-1}"
    local cits_pvalue="${6:-1}"
    local cits_gap="${7:-25}"
    local motif_flank="${8:-10}"
    local run_motif="${9:-yes}"
    local run_cims="${10:-true}"
    local run_cits="${11:-true}"
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK CIMS/CITS FULL ANALYSIS"
    log_info "  CIMS: $run_cims | CITS: $run_cits"
    log_info "═══════════════════════════════════════════════════════════════"
    
    # Create directories based on what's enabled
    mkdir -p "$output_dir/preprocessing"
    [[ "$run_cims" == "true" ]] && mkdir -p "$output_dir/CIMS"
    [[ "$run_cits" == "true" ]] && mkdir -p "$output_dir/CITS"
    [[ "$run_motif" == "yes" ]] && mkdir -p "$output_dir/motif_analysis"
    
    local sample_name=$(basename "${bam_file%.bam}" | sed 's/.Aligned.sortedByCoord.out//')
    
    # Phase 1: Preprocessing
    log_info "Phase 1: Parsing alignment..."
    local tags_bed="${output_dir}/preprocessing/${sample_name}_tags.bed"
    local mutation_file="${output_dir}/preprocessing/${sample_name}_mutations.txt"
    
    # run_parse_alignment handles calmd + parseAlignment.pl
    run_parse_alignment "$bam_file" "$tags_bed" "$mutation_file" "$(dirname "$genome_fasta")"
    
    if [[ ! -s "$tags_bed" ]]; then
        log_error "parseAlignment.pl failed to produce tags. Aborting CTK analysis."
        return 1
    fi
    
    # Phase 1b: Collapse tags
    log_info "Phase 1b: Collapsing PCR duplicates..."
    local collapsed_bed="${output_dir}/preprocessing/${sample_name}_collapsed.bed"
    tag2collapse.pl --keep-tag-name "$tags_bed" "$collapsed_bed"
    
    if [[ ! -s "$collapsed_bed" ]]; then
        log_error "tag2collapse.pl failed. Aborting CTK analysis."
        return 1
    fi
    
    # Phase 1c: CTK Preprocessing (selectRow + getMutationType)
    log_info "Phase 1c: Filtering and extracting mutation types..."
    run_ctk_preprocessing "$collapsed_bed" "$mutation_file" "${output_dir}/preprocessing"
    
    local del_bed="${output_dir}/preprocessing/deletions.bed"
    local sub_bed="${output_dir}/preprocessing/substitutions.bed"
    
    # Phase 2: CIMS Analysis (only if enabled)
    if [[ "$run_cims" == "true" ]]; then
        log_info "Phase 2: CIMS Analysis..."
        
        if [[ -s "$del_bed" ]]; then
            log_info "Running CIMS on deletions..."
            run_cims "$collapsed_bed" "$del_bed" \
                "${output_dir}/CIMS/${sample_name}_CIMS_del.txt" \
                "$cims_iterations" "$cims_fdr"
        fi
        
        if [[ -s "$sub_bed" ]]; then
            log_info "Running CIMS on substitutions..."
            run_cims "$collapsed_bed" "$sub_bed" \
                "${output_dir}/CIMS/${sample_name}_CIMS_sub.txt" \
                "$cims_iterations" "$cims_fdr"
        fi
    else
        log_info "Phase 2: CIMS Analysis... SKIPPED (not enabled)"
    fi
    
    # Phase 3: CITS Analysis (only if enabled)
    if [[ "$run_cits" == "true" ]]; then
        log_info "Phase 3: CITS Analysis..."
        
        if [[ -s "$del_bed" ]]; then
            run_cits "$collapsed_bed" "$del_bed" \
                "${output_dir}/CITS/${sample_name}_CITS.bed" \
                "$cits_pvalue" "$cits_gap"
        else
            log_warning "No deletion file for CITS. Skipping."
        fi
    else
        log_info "Phase 3: CITS Analysis... SKIPPED (not enabled)"
    fi
    
    # Phase 4: Flanked BED Generation (for user's motif analysis)
    if [[ "$run_motif" == "yes" ]]; then
        log_info "Phase 4: Generating flanked BED files..."
        
        if [[ "$run_cims" == "true" ]]; then
            local cims_del_sig="${output_dir}/CIMS/${sample_name}_CIMS_del_significant.bed"
            local cims_sub_sig="${output_dir}/CIMS/${sample_name}_CIMS_sub_significant.bed"
            [[ -s "$cims_del_sig" ]] && generate_flanked_bed "$cims_del_sig" "$motif_flank"
            [[ -s "$cims_sub_sig" ]] && generate_flanked_bed "$cims_sub_sig" "$motif_flank"
        fi
        
        if [[ "$run_cits" == "true" ]]; then
            local cits_singleton="${output_dir}/CITS/${sample_name}_CITS_singleton.bed"
            [[ -s "$cits_singleton" ]] && generate_flanked_bed "$cits_singleton" "$motif_flank"
        fi
    fi
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK ANALYSIS COMPLETE"
    log_info "  Output: $output_dir"
    log_info "═══════════════════════════════════════════════════════════════"
}

# 9. CTK Analysis Pipeline (Streamlined - uses pre-existing collapsed.bed and mutations.txt)
# This function reuses the standard pipeline's outputs to avoid duplicate preprocessing
run_ctk_analysis() {
    local collapsed_bed="$1"          # From standard pipeline
    local mutation_file="$2"          # From standard pipeline
    local output_dir="$3"
    local genome_fasta="$4"
    local sample_name="$5"
    local cims_iterations="${6:-10}"
    local cims_fdr="${7:-0.05}"
    local cits_pvalue="${8:-0.05}"
    local cits_gap="${9:-25}"
    local motif_flank="${10:-10}"
    local run_motif="${11:-yes}"
    local run_cims="${12:-true}"
    local run_cits="${13:-true}"
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK CIMS/CITS ANALYSIS"
    log_info "  Sample: $sample_name"
    log_info "  CIMS: $run_cims | CITS: $run_cits"
    log_info "═══════════════════════════════════════════════════════════════"
    
    # Create directories based on what's enabled
    [[ "$run_cims" == "true" ]] && mkdir -p "$output_dir/CIMS"
    [[ "$run_cits" == "true" ]] && mkdir -p "$output_dir/CITS"
    
    # Verify inputs exist
    if [[ ! -s "$collapsed_bed" ]]; then
        log_error "CTK: Collapsed BED file is empty or missing: $collapsed_bed"
        return 1
    fi
    if [[ ! -s "$mutation_file" ]]; then
        log_warning "CTK: Mutation file is empty or missing: $mutation_file"
        log_warning "CTK: CIMS/CITS requires mutation information. Skipping."
        return 1
    fi
    
    # Step 1: CTK Preprocessing (selectRow + getMutationType)
    log_info "Step 1: Filtering and extracting mutation types..."
    run_ctk_preprocessing "$collapsed_bed" "$mutation_file" "$output_dir"
    
    local del_bed="${output_dir}/deletions.bed"
    local sub_bed="${output_dir}/substitutions.bed"
    
    # Step 2: CIMS Analysis (only if enabled)
    if [[ "$run_cims" == "true" ]]; then
        update_status "CIMS"
        
        local cims_del_file="${output_dir}/CIMS/${sample_name}_CIMS_del.txt"
        local cims_sub_file="${output_dir}/CIMS/${sample_name}_CIMS_sub.txt"
        
        if [[ -s "$del_bed" ]]; then
            run_cims "$collapsed_bed" "$del_bed" "$cims_del_file" \
                "$cims_iterations" "$cims_fdr"
        fi
        
        if [[ -s "$sub_bed" ]]; then
            run_cims "$collapsed_bed" "$sub_bed" "$cims_sub_file" \
                "$cims_iterations" "$cims_fdr"
        fi
        
        # Generate flanked BED for CIMS results (for user's motif analysis)
        if [[ "$run_motif" == "yes" ]]; then
            [[ -s "$cims_del_file" ]] && generate_flanked_bed "$cims_del_file" "$motif_flank"
            [[ -s "$cims_sub_file" ]] && generate_flanked_bed "$cims_sub_file" "$motif_flank"
        fi
    fi
    
    # Step 3: CITS Analysis (only if enabled)
    if [[ "$run_cits" == "true" ]]; then
        update_status "CITS"
        
        local cits_file="${output_dir}/CITS/${sample_name}_CITS.txt"
        
        if [[ -s "$del_bed" ]]; then
            run_cits "$collapsed_bed" "$del_bed" "$cits_file" \
                "$cits_pvalue" "$cits_gap"
                
            # Generate flanked BED for CITS results (for user's motif analysis)
            if [[ "$run_motif" == "yes" ]]; then
                [[ -s "$cits_file" ]] && generate_flanked_bed "$cits_file" "$motif_flank"
            fi
        else
            log_warning "No deletion file for CITS. Skipping."
        fi
    fi
    
    log_info "═══════════════════════════════════════════════════════════════"
    log_info "  CTK ANALYSIS COMPLETE"
    log_info "  Output: $output_dir"
    log_info "═══════════════════════════════════════════════════════════════"
}

# 7. Coverage Analysis (Bedgraph)
# 7. Coverage Analysis (Bedgraph)
run_coverage() {
    local input_bed="$1"      
    local output_prefix="$2"
    local genome_file="$3"     
    local bam_file="$4"       # [NEW] Explicit BAM path

    update_status "Bedgraph"
    log_info "generating normalized, filtered bedGraphs..."
    
    # Validation
    if [[ ! -f "$bam_file" ]]; then
         log_error "Combined/Normalized BedGraph requires BAM input, but could not locate: $bam_file"
         return 1
    fi
    
    # 2. Calculate Scale Factor for Normalization
    # Count mapped reads (primary alignments only? or all mapped?)
    # usually -F 4 (mapped). User script used -F 4.
    local mapped=$(samtools view -c -F 4 "$bam_file")
    
    if [[ "$mapped" -eq 0 ]]; then
        log_warning "No mapped reads found in $bam_file. Skipping BedGraph."
        return 0
    fi
    
    local scale=$(echo "scale=6; 1000000 / $mapped" | bc)
    log_info "  Mapped Reads: $mapped | Scale Factor: $scale"
    
    # Store scale factor for later use (normalized counts)
    # Use same directory as bedgraph output (output_prefix parent) to avoid OUTPUT_ROOT scoping issues
    local bg_dir=$(dirname "$output_prefix")
    local scale_file="${bg_dir}/scale_factors.tsv"
    echo -e "$(basename "$output_prefix")\t${mapped}\t${scale}" >> "$scale_file"
    
    # 3. Generate BedGraph (Filtered & Normalized)
    # Using 'cigar !~ "N"' to remove junction reads
    # Streaming samtools -> bedtools genomecov
    
    # Positive Strand
    samtools view -h -e 'cigar !~ "N"' "$bam_file" | \
    bedtools genomecov -ibam stdin -bg -strand + -scale "$scale" > "${output_prefix}_pos.bedgraph"
    
    # Negative Strand
    samtools view -h -e 'cigar !~ "N"' "$bam_file" | \
    bedtools genomecov -ibam stdin -bg -strand - -scale "$scale" > "${output_prefix}_neg.bedgraph"
    
    
    log_info "Bedgraphs generated: ${output_prefix}_pos.bedgraph, ${output_prefix}_neg.bedgraph"
}

# 8. Combined BedGraph Generation (Group Averaging)
run_combined_bedgraph() {
    local output_dir="$1"
    local groups_file="$2"
    local bedgraph_dir="$3"
    
    log_info "Generating combined average bedgraphs..."
    
    mkdir -p "$bedgraph_dir/COMBINED_BEDGRAPH"
    
    # Identify Groups
    # If groups file provided, use it. Else, maybe regex?
    # For robust implementation, we'll assume GROUPS_FILE structure: SampleName<TAB>GroupName
    
    if [[ -z "$groups_file" || ! -f "$groups_file" ]]; then
        log_warning "No valid groups file provided. Skipping combined bedgraph generation."
        return 0
    fi
    
    # Extract unique groups
    local groups=$(awk '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$groups_file" | sort | uniq)
    
    for group in $groups; do
        log_info "Processing Group: $group"
        
        # Get samples for this group
        local samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file")
        
        # Process Positive and Negative strands separately
        for strand in "pos" "neg"; do
            local bg_files=""
            local count=0
            
            for sample in $samples; do
                # Construct expected filename from run_coverage output
                # {sample}_pos.bedgraph
                local f="$bedgraph_dir/${sample}_${strand}.bedgraph"
                if [[ -f "$f" ]]; then
                    bg_files="$bg_files $f"
                    ((count++))
                fi
            done
            
            if [[ $count -gt 0 ]]; then
                # Union and Average
                # Unionbedg produces: chrom start end val1 val2 ... valN
                # Column 1,2,3 are coords. Columns 4 to 3+N are values.
                # Average = sum(col 4..NF) / (NF-3)
                
                local output_file="$bedgraph_dir/COMBINED_BEDGRAPH/${group}_combined_${strand}.bedgraph"
                
                bedtools unionbedg -i $bg_files | \
                awk -v N="$count" 'BEGIN{OFS="\t"} {sum=0; for(i=4;i<=NF;i++) sum+=$i; print $1,$2,$3,sum/N}' \
                > "$output_file"
                
                log_info "  Generatd: $(basename "$output_file") ($count replicates)"
            else
                log_warning "  No bedgraph files found for group $group ($strand)"
            fi
        done
    done
}

# ═══════════════════════════════════════════════════════════════════════════
# Group-Based CTK Analysis (Bash 3.x Compatible)
# Aggregates samples by group before running CIMS/CITS
# ═══════════════════════════════════════════════════════════════════════════

run_group_ctk_analysis() {
    local groups_file="$1"
    local output_root="$2"
    local genome_index="$3"
    local cims_iterations="$4"
    local cims_fdr="$5"
    local cits_pvalue="$6"
    local cits_gap="$7"
    local motif_flank="$8"
    local run_motif="$9"
    local run_cims="${10}"
    local run_cits="${11}"
    
    console_msg "\n[GROUP CTK ANALYSIS]"
    
    # 1. Parse groups file into temp file (sample<TAB>group format)
    local groups_map=$(mktemp)
    parse_groups_file "$groups_file" "$groups_map"
    
    # 2. Create temp file for group→samples mapping
    local group_samples_file=$(mktemp)
    
    # 3. Find all sample analysis directories and assign to groups
    for sample_dir in *_analysis; do
        [[ ! -d "$sample_dir" ]] && continue
        
        local sample_name="${sample_dir%_analysis}"
        
        # Look up group for this sample using grep (bash 3.x compatible)
        local group=$(grep -w "^${sample_name}" "$groups_map" 2>/dev/null | cut -f2)
        
        if [[ -z "$group" ]]; then
            # Sample not in groups file → individual group (use sample name)
            group="$sample_name"
            log_info "Sample '$sample_name' not in groups file, treating as individual"
        fi
        
        # Append to group_samples_file: group<TAB>sample
        echo -e "${group}\t${sample_name}" >> "$group_samples_file"
    done
    
    # 4. Determine CTK output folder name
    local ctk_folder_name
    if [[ "$run_cims" == "true" ]] && [[ "$run_cits" == "true" ]]; then
        ctk_folder_name="5_CTK_Analysis"
    elif [[ "$run_cims" == "true" ]]; then
        ctk_folder_name="5_CIMS_Analysis"
    elif [[ "$run_cits" == "true" ]]; then
        ctk_folder_name="5_CITS_Analysis"
    fi
    
    # 5. Get unique groups and process alphabetically (skip unknown)
    local unique_groups=$(cut -f1 "$group_samples_file" | sort -u | grep -v "^unknown$")
    
    for group in $unique_groups; do
        # Get all samples for this group
        local samples=$(grep -w "^${group}" "$group_samples_file" | cut -f2 | tr '\n' ' ')
        local sample_count=$(echo $samples | wc -w | tr -d ' ')
        
        # Print group header on same line (no newline) so status updates appear after it
        if [[ "$sample_count" -eq 1 ]]; then
            printf "  > Processing %s: " "$group"
        else
            printf "  > Processing %s (%d samples): " "$group" "$sample_count"
        fi
        
        # 6. Create group CTK directory
        local group_ctk_dir="$output_root/$ctk_folder_name/$group"
        mkdir -p "$group_ctk_dir"
        
        # 7. Aggregate collapsed.bed files
        local group_collapsed="$group_ctk_dir/${group}_collapsed.bed"
        > "$group_collapsed"  # Clear file
        for sample in $samples; do
            # Use find to locate collapsed.bed in sample's analysis directory
            local sample_dir="${sample}_analysis"
            if [[ -d "$sample_dir" ]]; then
                local sample_collapsed=$(find "$sample_dir" -name "*_collapsed.bed" 2>/dev/null | head -n 1)
                if [[ -s "$sample_collapsed" ]]; then
                    cat "$sample_collapsed" >> "$group_collapsed"
                    log_info "Added $sample_collapsed to group $group"
                else
                    log_warning "Collapsed BED not found for sample: $sample"
                fi
            else
                log_warning "Analysis directory not found: $sample_dir"
            fi
        done
        
        # 8. Aggregate mutations.txt files
        local group_mutations="$group_ctk_dir/${group}_mutations.txt"
        > "$group_mutations"  # Clear file
        for sample in $samples; do
            local sample_dir="${sample}_analysis"
            if [[ -d "$sample_dir" ]]; then
                local sample_mutations=$(find "$sample_dir" -name "*_mutations.txt" 2>/dev/null | head -n 1)
                if [[ -s "$sample_mutations" ]]; then
                    cat "$sample_mutations" >> "$group_mutations"
                fi
            fi
        done
        
        # 9. Get genome fasta for motif analysis (prioritize genome/primary, exclude ncrna)
        local genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*genome*.fa" -o -name "*genome*.fasta" \) 2>/dev/null | head -n 1)
        if [[ -z "$genome_fasta" ]]; then
            genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*primary*.fa" -o -name "*primary*.fasta" \) 2>/dev/null | head -n 1)
        fi
        if [[ -z "$genome_fasta" ]]; then
            genome_fasta=$(find "$genome_index" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) ! -name "*ncrna*" 2>/dev/null | head -n 1)
        fi
        
        # 10. Run CTK analysis on aggregated data
        if [[ -s "$group_collapsed" ]]; then
            run_ctk_analysis "$group_collapsed" "$group_mutations" \
                "$group_ctk_dir" "$genome_fasta" "$group" \
                "$cims_iterations" "$cims_fdr" "$cits_pvalue" "$cits_gap" \
                "$motif_flank" "$run_motif" "$run_cims" "$run_cits"
        else
            log_error "CTK: Collapsed BED file is empty or missing: $group_collapsed"
        fi
        
        update_status_done
    done
    
    # Cleanup temp files
    rm -f "$groups_map" "$group_samples_file"
    
    console_msg "  > Group CTK analysis complete"
}
