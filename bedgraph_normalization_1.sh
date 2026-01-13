#!/bin/bash
## Competition CLIP Analysis: 
## Generate Normalized and Combined BedGraph for Sample Types
## Adapted from CoCLIP pipeline
## Last Edit: 2026-01-12

## Step 1: Generate mapped depth file from BAM files
## Step 2: Filter collapsed BED to peaks, generate strand-specific normalized bedgraphs
## Step 3: Combine replicates by group using unionbedg
## Step 4: Run bedgraph_normalization.R to average replicates

# Set directories
BASE_DIR="/Volumes/1TB_Data/Competition_CLIP/Luna_4460_Pool_output"
BAM_DIR="${BASE_DIR}/1_BAM"
BED_DIR="${BASE_DIR}/2_COLLAPSED_BED"
PEAKS_DIR="${BASE_DIR}/4_PEAKS/COMBINED_PEAKS"
OUTPUT_DIR="${BASE_DIR}/normalized_bedGraph"
GENOME_FAI="/Volumes/1TB_Data/Annotations/GRCh38.primary_assembly.genome.fa.fai"

# Create output directories
mkdir -p "${OUTPUT_DIR}/filtered_bed"
mkdir -p "${OUTPUT_DIR}/normalized_bedgraph"
mkdir -p "${OUTPUT_DIR}/normalized_combined_bedgraph"

# Define samples to process (excluding noRT controls)
SAMPLES=(
    "hnRNPC_inRBM_WT1"
    "hnRNPC_inRBM_WT2"
    "hnRNPC_inRBM_WT3"
    "hnRNPC_inRBM_MUT1"
    "hnRNPC_inRBM_MUT2"
    "hnRNPC_inRBM_MUT3"
    "RBM25_WT1"
    "RBM25_WT2"
    "RBM25_WT3"
    "RBM25_MUT1"
    "RBM25_MUT2"
    "RBM25_MUT3"
)

# Step 1: Generate depth file (normalization factors = 1e6 / mapped_reads)
echo "Generating mapped depth file..."
DEPTH_FILE="${OUTPUT_DIR}/mapped_depth.txt"
> "${DEPTH_FILE}"  # Clear/create file

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
    if [ -f "${BAM_FILE}" ]; then
        MAPPED_READS=$(samtools view -c -F 4 "${BAM_FILE}")
        NORM_FACTOR=$(echo "scale=10; 1000000 / ${MAPPED_READS}" | bc)
        echo -e "${SAMPLE}\t${NORM_FACTOR}" >> "${DEPTH_FILE}"
        echo "  ${SAMPLE}: ${MAPPED_READS} mapped reads, normalization factor = ${NORM_FACTOR}"
    else
        echo "  WARNING: BAM file not found for ${SAMPLE}"
    fi
done

# Step 2: Filter to peaks and generate strand-specific normalized bedgraphs
echo ""
echo "Filtering to peaks and generating normalized bedgraphs..."

while IFS=$'\t' read -r NAME DEPTH; do
    BED_FILE="${BED_DIR}/${NAME}.bed"
    
    if [ -f "${BED_FILE}" ]; then
        FILTERED_BED="${OUTPUT_DIR}/filtered_bed/${NAME}.filtered.bed"
        
        echo "Processing ${NAME}..."
        
        # Filter to peaks
        echo "  Filtering to peaks..."
        bedtools intersect -s -wa \
            -a "${BED_FILE}" \
            -b "${PEAKS_DIR}/peaks_Sorted.bed" \
            > "${FILTERED_BED}"
        
        # Generate strand-specific bedgraphs
        echo "  Generating normalized bedgraphs (scale factor: ${DEPTH})..."
        bedtools genomecov -split \
            -i "${FILTERED_BED}" \
            -g "${GENOME_FAI}" \
            -bg -strand + -scale ${DEPTH} \
            > "${OUTPUT_DIR}/normalized_bedgraph/${NAME}.filtered.pos.bedgraph"
        
        bedtools genomecov -split \
            -i "${FILTERED_BED}" \
            -g "${GENOME_FAI}" \
            -bg -strand - -scale ${DEPTH} \
            > "${OUTPUT_DIR}/normalized_bedgraph/${NAME}.filtered.rev.bedgraph"
    else
        echo "  WARNING: BED file not found for ${NAME}"
    fi
done < "${DEPTH_FILE}"

# Step 3: Combine replicates by group
echo ""
echo "Combining replicates by group..."

cd "${OUTPUT_DIR}/normalized_bedgraph"

# hnRNPC_inRBM25WT
echo "  Combining hnRNPC_inRBM25WT..."
bedtools unionbedg -i \
    hnRNPC_inRBM_WT1.filtered.pos.bedgraph \
    hnRNPC_inRBM_WT2.filtered.pos.bedgraph \
    hnRNPC_inRBM_WT3.filtered.pos.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/hnRNPC_inRBM25WT.combined.filtered.pos.bedgraph"

bedtools unionbedg -i \
    hnRNPC_inRBM_WT1.filtered.rev.bedgraph \
    hnRNPC_inRBM_WT2.filtered.rev.bedgraph \
    hnRNPC_inRBM_WT3.filtered.rev.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/hnRNPC_inRBM25WT.combined.filtered.rev.bedgraph"

# hnRNPC_inRBM25Mut
echo "  Combining hnRNPC_inRBM25Mut..."
bedtools unionbedg -i \
    hnRNPC_inRBM_MUT1.filtered.pos.bedgraph \
    hnRNPC_inRBM_MUT2.filtered.pos.bedgraph \
    hnRNPC_inRBM_MUT3.filtered.pos.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/hnRNPC_inRBM25Mut.combined.filtered.pos.bedgraph"

bedtools unionbedg -i \
    hnRNPC_inRBM_MUT1.filtered.rev.bedgraph \
    hnRNPC_inRBM_MUT2.filtered.rev.bedgraph \
    hnRNPC_inRBM_MUT3.filtered.rev.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/hnRNPC_inRBM25Mut.combined.filtered.rev.bedgraph"

# RBM25WT
echo "  Combining RBM25WT..."
bedtools unionbedg -i \
    RBM25_WT1.filtered.pos.bedgraph \
    RBM25_WT2.filtered.pos.bedgraph \
    RBM25_WT3.filtered.pos.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/RBM25WT.combined.filtered.pos.bedgraph"

bedtools unionbedg -i \
    RBM25_WT1.filtered.rev.bedgraph \
    RBM25_WT2.filtered.rev.bedgraph \
    RBM25_WT3.filtered.rev.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/RBM25WT.combined.filtered.rev.bedgraph"

# RBM25Mut
echo "  Combining RBM25Mut..."
bedtools unionbedg -i \
    RBM25_MUT1.filtered.pos.bedgraph \
    RBM25_MUT2.filtered.pos.bedgraph \
    RBM25_MUT3.filtered.pos.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/RBM25Mut.combined.filtered.pos.bedgraph"

bedtools unionbedg -i \
    RBM25_MUT1.filtered.rev.bedgraph \
    RBM25_MUT2.filtered.rev.bedgraph \
    RBM25_MUT3.filtered.rev.bedgraph \
    -header > "${OUTPUT_DIR}/normalized_combined_bedgraph/RBM25Mut.combined.filtered.rev.bedgraph"

echo ""
echo "Done! Proceed to bedgraph_normalization.R to finish generating the bedgraph files."
