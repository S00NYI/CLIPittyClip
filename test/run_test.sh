#!/bin/bash

# test/run_test.sh - Verification script for CLIPittyClip 3.0
# Usage: bash test/run_test.sh

# Exit on error
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
TEST_DIR="${REPO_DIR}/test_output"

echo "=========================================================="
echo "Running CLIPittyClip Verification"
echo "=========================================================="

# Check for tools
if ! command -v STAR &> /dev/null; then
    echo "Error: STAR not found. Please 'conda activate CLIPittyClip'"
    exit 1
fi

if ! command -v fastp &> /dev/null; then
    echo "Error: fastp not found. Please 'conda activate CLIPittyClip'"
    exit 1
fi

# Setup test directory
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Create tiny genome (2 chromosomes, unique sequences)
cat > genome.fa <<EOF
>chr1
ATGCGTACGTTAGCCTAGCTAGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTAG
>chr2
CGTACGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
EOF

echo "[2/4] Building STAR index..."
mkdir -p STAR_Index
STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir STAR_Index \
     --genomeFastaFiles genome.fa \
     --genomeSAindexNbases 4 \
     --limitGenomeGenerateRAM 50000000 > star_build.log 2>&1

echo "[3/4] Generating synthetic reads..."
# Generate reads matching the genome (unique match)
# Matches chr1 uniquely (first 21 bases of chr1 above)
cat > reads.fastq <<EOF
@READ1
ATGCGTACGTTAGCCTAGCTA
+
IIIIIIIIIIIIIIIIIIIII
@READ2
CGTACGATCGATCGATCGTAG
+
IIIIIIIIIIIIIIIIIIIII
EOF
gzip reads.fastq

echo "[4/4] Running CLIPittyClip pipeline..."
echo "Command: ${REPO_DIR}/CLIPittyClip.sh -i reads.fastq.gz -x STAR_Index -t 2 -u 4 --verbose"

bash "${REPO_DIR}/CLIPittyClip.sh" -i reads.fastq.gz -x STAR_Index -t 2 -u 4 --verbose

echo "=========================================================="
echo "Verification Complete!"
echo "Check output in: $TEST_DIR"
echo "=========================================================="

# ------------------------------------------------------------------
# Test 2: Sampling Mode
# ------------------------------------------------------------------
echo "=========================================================="
echo "Test 2: Sampling Mode (--sample 1)"
echo "=========================================================="
# Clean previous run
rm -rf reads_analysis

bash "${REPO_DIR}/CLIPittyClip.sh" -i reads.fastq.gz -x STAR_Index -t 2 -u 4 --sample 1 --verbose

if grep -q "Sampling first 1 reads" "reads_analysis/reads_analysis.log"; then
    echo "SUCCESS: Sampling mode detected."
else
    echo "FAILURE: Sampling mode NOT detected."
    exit 1
fi
