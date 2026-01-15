#!/bin/bash

# lib/modules.sh - Analysis modules for CLIPittyClip

source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# ═══════════════════════════════════════════════════════════════════════════
# Group File Parsing for CIMS/CITS Aggregation
# ═══════════════════════════════════════════════════════════════════════════

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

# 1. Preprocessing with fastp
run_preprocessing() {
    local input_file="$1"
    local output_prefix="$2"
    local umi_len="$3" 
    local adapter3="$4"
    local threads="$5"
    local sample_size="$6"
    local dedup_mode="$7"
    
    update_status_first "Preprocessing"
    log_info "Starting preprocessing with fastp..."
    
    # Define Fastp Command Function
    run_fastp_cmd() {
        local in_f="$1"
        local cmd="fastp -i ${in_f} -o ${output_prefix}_cleaned.fastq.gz \
            --thread ${threads} \
            --length_required 16 \
            --average_qual 30 \
            --html ${output_prefix}_fastp.html \
            --json ${output_prefix}_fastp.json"
            
        if [ -n "$adapter3" ]; then cmd+=" --adapter_sequence ${adapter3}"; fi
        if [ "$umi_len" -gt 0 ]; then cmd+=" --umi --umi_loc=read1 --umi_len=${umi_len} --umi_delim=#"; fi
        if [ "$sample_size" -gt 0 ]; then cmd+=" --reads_to_process $sample_size"; fi
        
        # Note: DEDUP is handled externally via seqkit before fastp if requested.
        # But if dedup was somehow passed to this function, we ignore fastp's internal dedup 
        # because seqkit is superior for this use case.
        
        log_info "Running: $cmd"
        execute_cmd "$cmd"
    }

    # Preparation: Handle Input File (Sampling / Dedup / Repair)
    # We create a pipeline of commands to prepare the input for fastp
    
    local final_input="$input_file"
    local temp_dedup=""
    
    # 1. Deduplication (SeqKit) - High Priority to reduce load
    if [ "$dedup_mode" == "true" ]; then
        update_status "  Deduplicating (SeqKit)" 
        log_info "Deduplication: Running seqkit rmdup (by sequence, ignoring quality)..."
        temp_dedup="${output_prefix}_dedup_temp.fastq.gz"
        
        local decompress="gzip -dc"
        if command -v pigz &> /dev/null; then decompress="pigz -dc"; fi
        
        # Use seqkit rmdup -s (by seq) -o output.gz
        seqkit rmdup -s -o "$temp_dedup" "$input_file" 2>> "${LOG_FILE}"
        
        if [[ $? -eq 0 && -s "$temp_dedup" ]]; then
             final_input="$temp_dedup"
             log_info "Deduplication complete."
             update_status "  Preprocessing (Fastp)"
        else
             log_warning "SeqKit Deduplication failed or produced empty file. Reverting to original input."
             rm -f "$temp_dedup"
        fi
    fi

    # 2. Run Fastp on prepared input
    run_fastp_cmd "$final_input"
    exit_code=$?

    # 3. Failure Recovery: Sanitize FASTQ (Empty Read / Truncation Fix)
    if [ $exit_code -ne 0 ]; then
        log_warning "fastp failed on $(basename "$final_input"). Attempting to repair malformed/empty reads..."
        
        local repaired_file="${output_prefix}_repaired.fastq.gz"
        local source_for_repair="$final_input" 
        # Note: If dedup failed, final_input is original. If dedup succeeded but fastp failed, final_input is deduped.
        # We repair whatever failed.
        
        local decompress="gzip -dc"
        if command -v pigz &> /dev/null; then decompress="pigz -dc"; fi
        
        # AWK Script: filter empty sequences and ensure 4-line records
        # Logic: Read 4 lines. If seq (line 2) > 0 length, print. 
        # Implicitly drops truncated records at end of file.
        $decompress "$source_for_repair" | awk '{h=$0; getline s; getline p; getline q; if (length(s) > 0 && length(q) > 0) {print h; print s; print p; print q}}' | gzip > "$repaired_file"
        
        if [[ -s "$repaired_file" ]]; then
            log_info "Repair successful. Retrying fastp with repaired file..."
            run_fastp_cmd "$repaired_file"
            exit_code=$?
            
            if [ $exit_code -eq 0 ]; then
                rm "$repaired_file"
            fi
        else
            log_error "Repair failed (empty output). The input file might be completely corrupt."
        fi
    fi
    
    # Cleanup temp dedup file
    if [[ -n "$temp_dedup" && -f "$temp_dedup" ]]; then
        rm "$temp_dedup"
    fi

    if [ $exit_code -eq 0 ]; then
        log_info "Preprocessing complete."
    else
        log_error "fastp failed even after repair attempt."
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

    local read_command="gzip -dc"
    if [[ "$input_fastq" != *.gz ]]; then
        read_command="cat"
    fi

    local cmd="STAR --runThreadN ${threads} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${input_fastq} \
        --readFilesCommand ${read_command} \
        --outFileNamePrefix ${output_prefix}. \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNmax ${mismatch_max} \
        --outSAMattributes NH HI AS nM NM MD \
        --alignEndsType EndToEnd \
        $ADV_ALIGNER_ARGS"

    log_info "Running: $cmd"
    execute_cmd "$cmd"
    if [ $? -ne 0 ]; then
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

# 3. PCR Duplicate Removal
run_collapse_pcr() {
    local input_bed="$1"
    local output_bed="$2"
    local umi_len="${3:-0}"  # Optional: UMI length, defaults to 0

    update_status "Collapsing"
    log_info "Collapsing PCR duplicates with CTK tag2collapse.pl..."
    
    # Only use --random-barcode when UMI is present in read names
    local barcode_flag=""
    if [ "${umi_len:-0}" -gt 0 ]; then
        barcode_flag="--random-barcode"
        log_info "UMI mode (length=$umi_len): using random barcode collapse"
    else
        log_info "No UMI: using position-only collapse"
    fi
    
    local cmd="$CONDA_PREFIX/bin/perl $(which tag2collapse.pl) --keep-tag-name --keep-max-score ${barcode_flag} \
        \"${input_bed}\" \"${output_bed}\""

    execute_cmd "$cmd"

    if [ $? -eq 0 ]; then
        log_info "Collapsing complete."
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
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file")
            
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
            # Get scale factor for this sample
            local sf=$(grep -P "^${sample}\t|/${sample}\t" "$scale_file" | tail -1 | cut -f3)
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
            local group_samples=$(awk -v g="$group" '{gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2)} $2==g {print $1}' "$groups_file")
            
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
        
        # Split peaks by strand and SORT for bedtools compatibility
        awk -F'\t' '$6=="+"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_pos.tmp.bed
        awk -F'\t' '$6=="-"' "$peaks_bed" | sort -k1,1 -k2,2n > peaks_neg.tmp.bed
        
        for sample in "${samples[@]}"; do
            local bg_pos="${bg_dir}/${sample}_pos.bedgraph"
            local bg_neg="${bg_dir}/${sample}_neg.bedgraph"
            
            if [[ ! -f "$bg_pos" || ! -f "$bg_neg" ]]; then
                log_warning "    BedGraph not found for $sample, skipping."
                continue
            fi
            
            # Sort bedgraphs for bedtools compatibility
            sort -k1,1 -k2,2n "$bg_pos" > "${sample}_pos_sorted.bg.tmp"
            sort -k1,1 -k2,2n "$bg_neg" > "${sample}_neg_sorted.bg.tmp"
            
            for stat in sum mean max; do
                # Map positive strand peaks to pos bedgraph
                bedtools map -a peaks_pos.tmp.bed -b "${sample}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 2>/dev/null > "bg_pos_${stat}.tmp"
                # Map negative strand peaks to neg bedgraph
                bedtools map -a peaks_neg.tmp.bed -b "${sample}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 2>/dev/null > "bg_neg_${stat}.tmp"
                
                # Combine and sort to match original peak order
                # Add line numbers for sorting
                awk -F'\t' '{print NR"\t"$0}' "$peaks_bed" > peaks_numbered.tmp
                awk -F'\t' 'NR==FNR {pos[$1"\t"$2"\t"$3]=$NF; next} {key=$1"\t"$2"\t"$3; print (key in pos) ? pos[key] : "0"}' "bg_pos_${stat}.tmp" peaks_pos.tmp.bed > "pos_vals.tmp"
                awk -F'\t' 'NR==FNR {neg[$1"\t"$2"\t"$3]=$NF; next} {key=$1"\t"$2"\t"$3; print (key in neg) ? neg[key] : "0"}' "bg_neg_${stat}.tmp" peaks_neg.tmp.bed > "neg_vals.tmp"
                
                # Merge based on strand in original peaks
                awk -F'\t' 'NR==FNR && $6=="+" {pos[$1"\t"$2"\t"$3]=FNR; next}
                             NR==FNR && $6=="-" {neg[$1"\t"$2"\t"$3]=FNR; next}
                             FNR==NR {posv[FNR]=$0; next}
                             {negv[FNR]=$0}
                             END {
                                 for(i=1; i<=NR; i++) {
                                     # placeholder
                                 }
                             }' "$peaks_bed" "pos_vals.tmp" "neg_vals.tmp"
                
                # Simpler approach: iterate through peaks and pick from correct file
                awk -F'\t' -v stat_val="$stat" -v sample_val="$sample" '
                BEGIN { 
                    while((getline < "bg_pos_"stat_val".tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                    while((getline < "bg_neg_"stat_val".tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                }
                NR==1 { print "BG" toupper(substr(stat_val,1,1)) substr(stat_val,2) "_" sample_val; next }
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
            
            rm -f "bg_pos_"*.tmp "bg_neg_"*.tmp pos_vals.tmp neg_vals.tmp peaks_numbered.tmp
            rm -f "${sample}_pos_sorted.bg.tmp" "${sample}_neg_sorted.bg.tmp"
        done
        
        rm -f peaks_pos.tmp.bed peaks_neg.tmp.bed
        
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
                
                # Sort group bedgraphs for bedtools compatibility
                sort -k1,1 -k2,2n "$grp_bg_pos" > "${group}_pos_sorted.bg.tmp"
                sort -k1,1 -k2,2n "$grp_bg_neg" > "${group}_neg_sorted.bg.tmp"
                
                for stat in sum mean max; do
                    bedtools map -a peaks_pos.tmp.bed -b "${group}_pos_sorted.bg.tmp" -c 4 -o "$stat" -null 0 2>/dev/null > "bg_pos_${stat}.tmp"
                    bedtools map -a peaks_neg.tmp.bed -b "${group}_neg_sorted.bg.tmp" -c 4 -o "$stat" -null 0 2>/dev/null > "bg_neg_${stat}.tmp"
                    
                    awk -F'\t' '
                    BEGIN { 
                        while((getline < "bg_pos_'"$stat"'.tmp") > 0) { pos[$1"\t"$2"\t"$3] = $NF }
                        while((getline < "bg_neg_'"$stat"'.tmp") > 0) { neg[$1"\t"$2"\t"$3] = $NF }
                    }
                    NR==1 { print "BG" toupper(substr("'$stat'",1,1)) substr("'$stat'",2) "_'$group'"; next }
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
            grep -v "^#" "$ctk_file" | awk -F'\t' -v fdr="$threshold" \
                '$9+0 < fdr {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$temp_filtered"
        else
            # CITS: P-value is embedded in name as [P=value]
            grep -v "^#" "$ctk_file" | awk -F'\t' -v pval="$threshold" '{
                if (match($4, /\[P=([0-9.e+-]+)\]/, arr)) {
                    if (arr[1]+0 < pval) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
                }
            }' > "$temp_filtered"
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
            local samples=$(grep -P "\t${group}$" "$groups_map" | cut -f1)
            
            # --- Aggregate CIMS (Deletions) ---
            local group_del_bed="${ctk_dir}/${group}_aggregated_CIMS_del.txt"
            > "$group_del_bed"
            for sample in $samples; do
                # Look for sample CIMS file in CTK dir or subdirs
                local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_del.txt" 2>/dev/null | head -n 1)
                [[ -s "$s_file" ]] && cat "$s_file" >> "$group_del_bed"
            done
            if [[ -s "$group_del_bed" ]]; then
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            else
                # still add empty column for consistency?
                add_ctk_column "$group_del_bed" "DEL_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CIMS (Substitutions) ---
            local group_sub_bed="${ctk_dir}/${group}_aggregated_CIMS_sub.txt"
            > "$group_sub_bed"
            for sample in $samples; do
                local s_file=$(find "$ctk_dir" -name "${sample}_CIMS_sub.txt" 2>/dev/null | head -n 1)
                [[ -s "$s_file" ]] && cat "$s_file" >> "$group_sub_bed"
            done
            if [[ -s "$group_sub_bed" ]]; then
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            else
                add_ctk_column "$group_sub_bed" "SUB_${group}" "cims" "$cims_fdr"
            fi
            
            # --- Aggregate CITS (Truncations) ---
            local group_cits_bed="${ctk_dir}/${group}_aggregated_CITS.txt"
            > "$group_cits_bed"
            for sample in $samples; do
                local s_file=$(find "$ctk_dir" -name "${sample}_CITS.txt" 2>/dev/null | head -n 1)
                [[ -s "$s_file" ]] && cat "$s_file" >> "$group_cits_bed"
            done
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
                add_ctk_column "$cims_del_file" "${name}_del" "cims" "$cims_fdr"
            fi
        done
        
        for cims_sub_file in "${ctk_dir}/CIMS/"*_CIMS_sub.txt; do
            if [[ -f "$cims_sub_file" ]]; then
                local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                add_ctk_column "$cims_sub_file" "${name}_sub" "cims" "$cims_fdr"
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
                        add_ctk_column "$cims_del_file" "${name}_del" "cims" "$cims_fdr"
                    fi
                done
                
                for cims_sub_file in "${group_dir}CIMS/"*_CIMS_sub.txt; do
                    if [[ -f "$cims_sub_file" ]]; then
                        local name=$(basename "$cims_sub_file" _CIMS_sub.txt)
                        add_ctk_column "$cims_sub_file" "${name}_sub" "cims" "$cims_fdr"
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
                add_ctk_column "$cits_file" "${name}_trunc" "cits" "$cits_pval"
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
                        add_ctk_column "$cits_file" "${name}_trunc" "cits" "$cits_pval"
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
# Verified workflow based on Zhang Lab CTK documentation (Dec 2024)
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

# 5. CIMS Analysis (CTK)
# Detects crosslinking-induced mutation sites
# Input: collapsed BED + mutation file (deletions or substitutions)
# Output: CIMS.txt with significant mutation sites
run_cims() {
    local input_collapsed_bed="$1"
    local mutation_bed="$2"          # Already BED6 from getMutationType.pl
    local output_file="$3"
    local cims_iterations="${4:-10}"  # Default: 10 iterations
    local cims_fdr="${5:-0.001}"      # Default: FDR 0.001
    
    # Status update moved to caller (run_ctk_analysis)
    # update_status "CIMS Analysis"
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
    
    # Run CIMS.pl
    local cmd="CIMS.pl -v -n $cims_iterations '$input_collapsed_bed' '$mutation_bed' '$output_file'"
    
    execute_cmd "$cmd"
    local cims_exit=$?
    
    # CIMS.pl may exit with error on sparse data but still produce valid output
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

# 6. CITS Analysis (CTK)
# Detects crosslinking-induced truncation sites
# Input: collapsed BED + deletion file (used to EXCLUDE read-through tags)
# Output: CITS.txt with significant truncation sites (single-nucleotide singletons)
run_cits() {
    local input_collapsed_bed="$1"
    local deletion_bed="$2"           # Used to exclude read-through tags
    local output_file="$3"            # Should be .txt extension
    local cits_pvalue="${4:-0.001}"   # Default: p-value 0.001
    local cits_gap="${5:-25}"         # Default: gap 25 for clustering
    
    # Status update moved to caller (run_ctk_analysis)
    # update_status "CITS Analysis"
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
    
    # Run CITS.pl - outputs intermediate file first
    local cits_raw="${output_file%.txt}_raw.bed"
    local cmd="CITS.pl -p $cits_pvalue --gap $cits_gap -v '$input_collapsed_bed' '$deletion_bed' '$cits_raw'"
    
    execute_cmd "$cmd"
    local cits_exit=$?
    
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

# 7. Motif Enrichment Analysis (HOMER)
# Extracts flanking regions and runs findMotifsGenome.pl
run_motif_analysis() {
    local input_bed="$1"
    local output_dir="$2"
    local genome_fasta="$3"
    local flank_nt="${4:-10}"         # Default: ±10 nucleotides
    local label="${5:-motif}"
    
    # Status update removed - motif analysis is part of CIMS/CITS
    log_info "Running motif enrichment analysis..."
    log_info "Input: $input_bed"
    log_info "Flanking region: ±${flank_nt}nt"
    
    if [[ ! -s "$input_bed" ]]; then
        log_warning "Motif analysis: Input BED is empty, skipping."
        return 1
    fi
    
    mkdir -p "$output_dir"
    
    # Extend regions by ±flank_nt
    local flanked_bed="${output_dir}/${label}_flanked.bed"
    awk -v n="$flank_nt" 'BEGIN{OFS="\t"} {
        start = $2 - n
        if (start < 0) start = 0
        print $1, start, $3 + n, $4, $5, $6
    }' "$input_bed" > "$flanked_bed"
    
    local site_count=$(wc -l < "$flanked_bed")
    log_info "Prepared $site_count sites with ±${flank_nt}nt flanks"
    
    # Run HOMER findMotifsGenome.pl
    local homer_output="${output_dir}/${label}_homer"
    
    if command -v findMotifsGenome.pl &> /dev/null; then
        log_info "Running HOMER findMotifsGenome.pl (RNA mode)..."
        
        # HOMER requires genome to be in its database or a FASTA path
        findMotifsGenome.pl "$flanked_bed" "$genome_fasta" "$homer_output" \
            -size given -rna -p "${THREADS:-1}" 2>> "${LOG_FILE:-/dev/null}"
        
        if [[ -d "$homer_output" ]]; then
            log_info "HOMER complete: $homer_output"
        else
            log_warning "HOMER may have failed. Check logs."
        fi
    else
        log_warning "HOMER (findMotifsGenome.pl) not found. Skipping motif analysis."
        log_info "To install: conda install -c bioconda homer"
    fi
    
    log_info "Motif analysis complete."
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
    
    # Phase 4: Motif Analysis
    if [[ "$run_motif" == "yes" ]]; then
        log_info "Phase 4: Motif Enrichment Analysis..."
        
        if [[ "$run_cims" == "true" ]]; then
            local cims_del_sig="${output_dir}/CIMS/${sample_name}_CIMS_del_significant.bed"
            local cims_sub_sig="${output_dir}/CIMS/${sample_name}_CIMS_sub_significant.bed"
            
            if [[ -s "$cims_del_sig" ]]; then
                run_motif_analysis "$cims_del_sig" "${output_dir}/motif_analysis" \
                    "$genome_fasta" "$motif_flank" "CIMS_del"
            fi
            
            if [[ -s "$cims_sub_sig" ]]; then
                run_motif_analysis "$cims_sub_sig" "${output_dir}/motif_analysis" \
                    "$genome_fasta" "$motif_flank" "CIMS_sub"
            fi
        fi
        
        if [[ "$run_cits" == "true" ]]; then
            local cits_singleton="${output_dir}/CITS/${sample_name}_CITS_singleton.bed"
            
            if [[ -s "$cits_singleton" ]]; then
                run_motif_analysis "$cits_singleton" "${output_dir}/motif_analysis" \
                    "$genome_fasta" "$motif_flank" "CITS"
            fi
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
        
        # Run motif on CIMS results (if enabled) - per-sample
        if [[ "$run_motif" == "yes" ]]; then
            if [[ -s "$cims_del_file" ]]; then
                run_motif_analysis "$cims_del_file" "${output_dir}/CIMS/CIMS_del_motif" \
                    "$genome_fasta" "$motif_flank" "CIMS_del"
            fi
            
            if [[ -s "$cims_sub_file" ]]; then
                run_motif_analysis "$cims_sub_file" "${output_dir}/CIMS/CIMS_sub_motif" \
                    "$genome_fasta" "$motif_flank" "CIMS_sub"
            fi
        fi
    fi
    
    # Step 3: CITS Analysis (only if enabled)
    if [[ "$run_cits" == "true" ]]; then
        update_status "CITS"
        
        local cits_file="${output_dir}/CITS/${sample_name}_CITS.txt"
        
        if [[ -s "$del_bed" ]]; then
            run_cits "$collapsed_bed" "$del_bed" "$cits_file" \
                "$cits_pvalue" "$cits_gap"
                
            # Run motif on CITS results (if enabled) - per-sample
            if [[ "$run_motif" == "yes" ]]; then
                if [[ -s "$cits_file" ]]; then
                    run_motif_analysis "$cits_file" "${output_dir}/CITS/CITS_motif" \
                        "$genome_fasta" "$motif_flank" "CITS"
                fi
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
    
    update_status "Combined BedGraph"
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
                awk -v N="$count" '{sum=0; for(i=4;i<=NF;i++) sum+=$i; print $1,$2,$3,sum/N}' \
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
