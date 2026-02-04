<p align="center">
  <img src="logo.png" alt="CLIPittyClip Logo" width="400">
</p>

# CLIPittyClip: Modern CLIP-seq Analysis Pipeline
**Version 3.0.0**

A comprehensive, single-command CLIP-seq data analysis pipeline from FASTQ to peaks and crosslink sites.

## Overview

CLIPittyClip v3.0 provides a complete, modernized workflow for CLIP-seq analysis:

- **Dual Aligner Support**
  - STAR (default) or Bowtie2
- **Advanced Coverage Analysis**
  - Normalized BedGraph and group averaging
- **Advanced Peak Calling**
  - Aggregated (meta-peak) or individual peak calling support
- **Group/Condition-Based Analysis**
  - Group/Condition-wise peak aggregation and CTK analysis
- **Interactive Wizard**
  - `--wizard` mode for guided configuration

## Pipeline Flow

<p align="center">
  <img src="flowchart.png" alt="CLIPittyClip Pipeline Flow" width="800">
</p>

## Installation

> ⚠️ **Development Version**: This is the `v3-development` branch. For the stable release, use the `main` branch.

### Quick Install (Recommended)

CLIPittyClip provides self-contained installation scripts that automatically install all dependencies including CTK, HOMER, and required Perl modules.

#### 1. Clone Repository
```bash
git clone -b v3-development https://github.com/S00NYI/CLIPittyClip.git
cd CLIPittyClip
```

#### 2. Run Installation Script

**macOS (Intel or Apple Silicon):**
```bash
./install_macos.sh --env clipittyclip --tools-dir ~/Tools
```

**Linux:**
```bash
./install_linux.sh --env clipittyclip --tools-dir ~/Tools
```

> **Note:** The scripts will:
> - Create a conda environment with all dependencies
> - Install CTK and HOMER from source
> - Install required Perl modules (Math::CDF, Bio::SeqIO) via CPAN
> - Configure your shell PATH automatically

#### 3. Activate and Verify
```bash
# Restart terminal or source your shell config
source ~/.zshrc    # macOS
source ~/.bashrc   # Linux

# Activate the environment
conda activate clipittyclip

# Verify installation
which CLIPittyClip.sh
which parseAlignment.pl
which findPeaks
perl -MBio::SeqIO -e 'print "Bio::SeqIO OK\n"'
perl -MMath::CDF -e 'print "Math::CDF OK\n"'
```

### Installation Options

| Option | Default | Description |
|--------|---------|-------------|
| `--env <name>` | `clipittyclip` | Conda environment name |
| `--tools-dir <path>` | `~/Tools` | Directory for CTK/HOMER |
| `--help` | — | Show help message |

### Manual Installation (Advanced)

If the automated scripts don't work for your system, see [docs/manual_installation.md](docs/manual_installation.md) for step-by-step instructions.


## Quick Start

```bash
# Basic analysis with STAR (single file)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 -u 7

# With demultiplexing
CLIPittyClip.sh -i pool.fastq.gz -b barcodes.txt -x /path/to/star_index -t 8

# Pre-demultiplexed samples in a folder
CLIPittyClip.sh -d /path/to/samples_folder/ -x /path/to/star_index -t 8

# Using Bowtie2 instead
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/bt2_index -t 8 -m bowtie2

# With CIMS analysis only
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 --run-cims

# With CITS analysis only
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 --run-cits

# With both CIMS and CITS analysis
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 --run-ctk
# OR equivalently:
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 --run-cims --run-cits

# ENCODE eCLIP analysis (pre-processed files with UMI in header)
CLIPittyClip.sh --eclip -d /path/to/samples/ -x /path/to/star_index -t 8 --run-cims --run-cits
```

## ENCODE eCLIP Mode

For pre-processed ENCODE eCLIP data, use `--eclip` mode:

```bash
# eCLIP preprocessing only
CLIPittyClip.sh --eclip -d /path/to/eclip_fastqs/ -x /path/to/star_index -t 8

# eCLIP with CIMS/CITS analysis (optional, add flags as needed)
CLIPittyClip.sh --eclip -d /path/to/eclip_fastqs/ -x /path/to/star_index -t 8 --run-cims --run-cits
```

**What `--eclip` does (preprocessing only):**
- **Skips UMI extraction** - UMI is already in read header (ENCODE format: `@NCCTGAATGA:...`)
- **Uses 9 standard eCLIP adapters** - Automatically trims all adapter variants
- **Reformats headers for CTK** - Converts to CTK-compatible format for tag2collapse

> **Note:** `--eclip` only affects preprocessing. CIMS/CITS analysis requires separate `--run-cims` and/or `--run-cits` flags.

> **Note:** Dynamic thread scaling for CIMS/CITS (based on available RAM) applies to all modes, not just eCLIP.

**When to use:**
- ENCODE eCLIP data downloaded from `encodeproject.org`
- Files with UMI already moved to read ID (format: `@UMI:REST_OF_ID`)
- Pre-demultiplexed eCLIP samples

**Ignored options in eCLIP mode:**
- `-u` (UMI length) - Detected from header
- `-a` (adapter) - Uses all 9 eCLIP adapters

## Input Modes

| Mode | Flag | Use Case |
|------|------|----------|
| Single file | `-i sample.fastq.gz` | One FASTQ, direct analysis |
| Pooled + barcodes | `-i pool.fastq.gz -b barcodes.txt` | Demultiplex then analyze |
| Pre-demuxed folder | `-d /path/to/folder/` | Batch analyze existing FASTQs |

## Command-Line Options

Run `CLIPittyClip.sh --help` for full usage.

### Input/Output Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-i` | `--input-file` | path | — | Input FASTQ file (required unless using `-d`) |
| `-d` | `--input-dir` | path | — | Input directory with pre-demultiplexed FASTQs |
| `-x` | `--index` | path | — | Genome index directory (required) |
| `-o` | `--output` | string | next to input | Output folder name or full path (see below) |
| `-k` | `--keep` | bool | no | Keep intermediate files |

> **Output Location (`-o`):**
> - **No `-o`**: Output created next to input files (e.g., `/data/reads.fq.gz` → `/data/reads_output/`)
> - **Name only** (`-o HepG2`): Creates folder next to input (e.g., `/data/HepG2/`)
> - **Full path** (`-o /results/HepG2`): Uses exact path specified

### Processing Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-t` | `--threads` | int | 1 | Number of threads |
| `-m` | `--mapper` | string | star | Mapper: `star` or `bowtie2` |
| `-v` | `--verbose` | bool | false | Enable verbose logging |
| `-h` | `--help` | — | — | Show help message |

### Preprocessing Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-u` | `--umi-length` | int | — | UMI length (auto-detected for eCLIP) |
| `-a` | `--adapter` | string | L32 | 3' adapter sequence |
| — | `--eclip` | bool | false | ENCODE eCLIP mode |
| — | `--no-dedup` | bool | — | Disable FASTQ deduplication (default: ON) |

### Demultiplexing Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-b` | `--barcodes` | path | — | Barcode file (enables demultiplexing) |
| — | `--demux-mismatches` | int | 1 | Max barcode mismatches |
| — | `--align-mismatches` | int | 2 | Max alignment mismatches (STAR only) |
| — | `--skip-ncrna` | bool | false | Disable ncRNA pre-filtering |

### CTK Analysis Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| — | `--run-ctk` | bool | — | Enable full CTK CIMS+CITS analysis |
| — | `--run-cims` | bool | — | Enable CIMS analysis (mutation sites) |
| — | `--run-cits` | bool | — | Enable CITS analysis (truncation sites) |
| — | `--cims-iter` | int | 5 | CIMS permutation iterations |
| — | `--cims-fdr` | float | 0.05 | CIMS FDR threshold |
| — | `--cits-pval` | float | 0.05 | CITS p-value threshold |
| — | `--cits-gap` | int | 25 | CITS clustering gap (-1 disables) |
| `-f` | `--flank` | int | 10 | Flanked BED nucleotides |
| — | `--no-motif` | bool | — | Skip flanked BED generation |

### Grouping Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-g` | `--groups` | path | — | Groups file for bedgraph/peak grouping |
| — | `--ctk-group` | bool | false | Enable group CTK analysis |

### Other Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| — | `--notification` | bool | false | Enable system notifications |
| `-w` | `--wizard` | bool | — | Launch interactive wizard |
| `-s` | `--sample` | int | — | Test mode: process only first N reads |

## Group-Based CIMS/CITS Analysis

Use `--ctk-group` with `-g groups.txt` to aggregate replicates/samples before running CIMS/CITS:

```bash
CLIPittyClip.sh -i pool.fq.gz -b barcodes.txt -x index --run-cims --run-cits -g groups.txt --ctk-group
```

**groups.txt format** (tab-separated: `SampleName<TAB>GroupName`):
```
Condition_A_Rep1    Condition_A
Condition_A_Rep2    Condition_A
Condition_B_Rep1    Condition_B
Condition_B_Rep2    Condition_B
```

> **Note:** Samples not listed in the groups file are treated as individual groups (analyzed separately).

> [!WARNING]
> **Memory Requirements:** Group-based CTK analysis aggregates all samples in each group before running CIMS/CITS, which can require significant memory (>64GB RAM recommended for large datasets). If you encounter memory issues or system crashes, consider:
> - Running CTK on individual samples instead of groups
> - Using a machine with more RAM

## Output Structure

```
{INPUT}_output/
├── 0_DEMUX_FASTQ/           # Demultiplexed reads
├── 1_BAM/                   # Aligned BAM files
├── 2_COLLAPSED_BED/         # PCR-deduplicated BED
├── 3_BEDGRAPH/              # Coverage tracks (Normalized Coverage)
│   ├── COMBINED_BEDGRAPH/   # (Group mode) Averaged replicates
├── 4_PEAKS/                 # HOMER peak results
│   ├── COMBINED_PEAKS/      # Peaks across all samples (or aggregated)
│   └── SAMPLE_PEAKS/        # Peaks per sample
│
├── 5_CTK_Analysis/          # When --run-ctk
│   ├── CIMS/
│   ├── CITS/
│   └── motif_analysis/
│ OR
├── 5_CIMS_Analysis/         # When --cims only
│ OR
├── 5_CITS_Analysis/         # When --cits only
│
├── 6_OTHERS/                # Support files
│   ├── STAR_OUTPUT/         # Splice junctions (STAR only)
│   └── ncRNA_Mapping/       # ncRNA mapping results
│
└── REPORTS/                 # Logs and QC reports
    ├── FASTP_REPORT/        # HTML/JSON QC reports
    ├── ALIGNER_LOGS/        # Aligner summaries (STAR/Bowtie2)
    ├── SAMPLES/             # Detailed per-sample logs
    ├── PEAK/                # Peak calling logs
    └── BEDGRAPH/            # BedGraph generation logs (optional)

{INPUT_BASENAME}.log         # Complete console log of the run

```

## Peak Matrix Metrics

The `COMBINED_peakMatrix.txt` file contains up to 54+ advanced metrics summarizing peak intensity, reproducibility, and coverage shape.

| Metric Type | Column Prefix | Scoping | Calculation Method | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Biological Complexity** | `BC_` | **Group** | Count of samples in the group with Raw Count (TC) > 0. | Tells you how reproducible the peak is across replicates. |
| **Total Count (Raw)** | `TC_` | **Sample** | Raw number of read starts (tags) overlapping the peak. | Absolute number of RBP-bound RNA (PCR-duplicates removed). |
| **Total Count (Raw Aggregate)** | `TC_` | **Group** | **SUM** of raw counts for all samples in the group. | Total RBP-bound RNA for the peak within that biological condition/group. |
| **Normalized Count** | `NormedTC_` | **Sample** | `Raw TC * (1,000,000 / Total Mapped Reads)`. | Read abundance relative to the total depth of the sample (RPM). |
| **Normalized Count Aggregate** | `NormedTC_` | **Group** | **SUM** of normalized counts for all samples in the group. | Total intensity of the peak for that group, standardized for sequencing depth. |
| **Coverage Sum** | `CovSum_` | **Sample** | Sum of per-base signal values from a **Normalized BedGraph** across the peak. | Read coverage across the peak per sample. |
| **Coverage Sum (Group Average)** | `CovSum_` | **Group** | Sum of per-base signal from an **Averaged Group BedGraph** across the peak. | Averaged read coverage across the peak for the group. |
| **Coverage Mean** | `CovMean_` | **Both** | `CovSum / Peak Length`. | Average signal intensity per base within the peak. |
| **Coverage Max** | `CovMax_` | **Both** | The **highest** signal value found at any single base within the peak. | Peak height. |

> [!TIP]
> **NormedTC** (additive) measures the total intensity of the group, while **Cov** columns (mean-based) measure the density/shape of the average replicate in that group.

## Console Output

```
[DEDUPLICATING]
  > Deduplicating Pooled Reads (SeqKit)
  > Deduplicating Complete

[DEMULTIPLEXING]
  > Barcodes: barcodes.txt
  > Mismatches Allowed: 2
  > Checking barcodes...
  > All barcodes are unique with 2 mismatches.
  > Demultiplexing Complete

[BATCH ANALYSIS]
   1/3 Sample1 : Preprocessing > Mapping (STAR) > Processing Alignment > Collapsing > Bedgraph > Peaks > CIMS > CITS > Done!
```

## Standalone Tools

### MAPittyMap.sh
Standalone mapping module for aligning FASTQ files to a reference genome.

**Required inputs:**
- `-i <path>`: Input FASTQ file (gzipped)
- `-x <path>`: Path to genome index directory

**Key options:**
- `-t <int>`: Number of threads (default: 1)
- `-m, --mapper <star|bowtie2>`: Alignment tool (default: star)
- `-o <path>`: Output directory
- `-w, --wizard`: Launch interactive configuration wizard

```bash
# Using STAR
MAPittyMap.sh -i reads.fastq.gz -x /path/to/star_index -t 8 -m star

# Using Bowtie2
MAPittyMap.sh -i reads.fastq.gz -x /path/to/bt2_index -t 8 -m bowtie2

# Interactive wizard mode for custom aligner settings
MAPittyMap.sh -i reads.fastq.gz -x /path/to/star_index -t 8 -w
```

---

### PEAKittyPeak.sh
Standalone peak calling using HOMER. Requires a directory containing collapsed BED files.

**Required inputs:**
- Run from a directory containing a `BED/` folder with `.bed` files
- OR use `-i <directory>` to specify input BED directory explicitly

**Key options:**
- `-p <int>`: Min distance between peaks (default: 50)
- `-z <int>`: Peak size (default: 20)
- `-f <int>`: Fragment length (default: 25)
- `-n <string>`: Output name prefix
- `-a <string>`: Additional HOMER findPeaks arguments
- `--aggregate`: Combine all input BED files into a single meta-sample for peak calling
- `--no-aggregate`: Process each BED file individually (Default)
- `--ctk-dir <path>`: Add CIMS/CITS site counts from CTK analysis
- `--ctk-group <file>`: Groups file for CTK aggregation
- `--wizard`: Launch interactive HOMER configuration wizard

```bash
# Individual Peak Calling (Default)
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Analysis --no-aggregate

# Aggregated Peak Calling (Meta-sample)
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Combined --aggregate

# With custom HOMER arguments
PEAKittyPeak.sh -n Combined -a '-style factor -L 2'

# With CTK site counts (adds _del, _sub, _trunc columns)
PEAKittyPeak.sh -n Combined --ctk-dir ./CTK_Analysis/

# Interactive wizard mode for HOMER settings
PEAKittyPeak.sh --wizard
```

## Generating Genome Indices

### STAR Index
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 100
```

### Bowtie2 Index
```bash
bowtie2-build genome.fa /path/to/bt2_index/GRCh38
```

### ncRNA Pre-filtering Index (Recommended)

CLIPittyClip can automatically filter ncRNA reads (rRNA, tRNA, snRNA, snoRNA) before genome alignment to improve peak calling accuracy. This is **enabled by default**.

**Setup:** Place a Bowtie2 index with prefix `ncrna` in either:
- **Recommended:** `<annotation_dir>/ncRNA/` subfolder
- **Legacy:** Directly in `<annotation_dir>/` (same location as `-x`)

**Building the ncRNA index:**

1. Download ncRNA sequences from Ensembl:
```bash
# Human (GRCh38)
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

# Mouse (GRCm39)
wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz
gunzip Mus_musculus.GRCm39.ncrna.fa.gz
```

2. Build Bowtie2 index (recommended: place in ncRNA subfolder):
```bash
mkdir -p /path/to/annotation/ncRNA
bowtie2-build Homo_sapiens.GRCh38.ncrna.fa /path/to/annotation/ncRNA/ncrna
```

**Behavior:**
- **Prioritization**: When searching for genome indices, Bowtie2/STAR wrappers exclude indices matching `*ncrna*` to prevent accidental alignment to the ncRNA subset (fixing 0% alignment issues).
- **Filtering**: Before main alignment, reads are mapped against `<annotation_dir>/ncRNA/ncrna.1.bt2`.
- **Output**: Unfiltered (non-ncRNA) reads continue to genome alignment. ncRNA stats are saved to `REPORTS/ncRNA_Mapping/`.

**To disable:** Use `--skip-ncrna` flag

### Annotation Directory Structure

For CTK CIMS/CITS motif analysis, the annotation directory (`-x`) should contain:

```
/path/to/annotation/
├── GRCh38.primary_assembly.genome.fa    # Reference FASTA (required for motif analysis)
├── Genome                               # STAR index files...
├── SA
├── SAindex
├── genomeParameters.txt
├── chrom.sizes                          # Optional: for bedgraph generation
└── ncRNA/                               # Recommended subfolder
    ├── ncrna.1.bt2                      # Bowtie2 ncRNA index
    ├── ncrna.2.bt2
    ├── ncrna.3.bt2
    ├── ncrna.4.bt2
    ├── ncrna.rev.1.bt2
    └── ncrna.rev.2.bt2
```

> [!NOTE]
> The reference FASTA for motif analysis is auto-detected with priority: `*genome*.fa` > `*primary*.fa` > any `.fa` excluding `*ncrna*`


## License

GPL-3.0 License - See [LICENSE](LICENSE) for details.
