<p align="center">
  <img src="logo.png" alt="CLIPittyClip Logo" width="400">
</p>

# CLIPittyClip: Modern CLIP-seq Analysis Pipeline
**Version 3.3.0**

A comprehensive, single-command CLIP-seq data analysis pipeline from FASTQ to peaks and crosslink sites.

## Update Log (v3.3.0)
This release adds the **Clink** pipeline: a Python-native crosslink site caller using `umi_tools` + `pysam`, contributed by the LunaRNALab fork. Clink runs in parallel with CTK and is independent of Perl infrastructure.

**New Features (v3.3.0):**
- **Clink Pipeline** (`--run-clink`): Python-native crosslink-induced truncation site (CITS) and deletion site (CIMS) caller. Uses `umi_tools directional` BAM deduplication, a pileup engine (`lib/clink/pileup.py`), and separate CITS (`lib/clink/cits.py`) / CIMS (`lib/clink/cims.py`) callers with Benjamini-Hochberg FDR control.
- **Clink-Only Fast Path**: When `--run-clink` is used without CTK flags, `parseAlignment.pl`/`tag2collapse.pl` are skipped entirely. `umi_tools dedup` → `bedtools bamtobed` produces the collapsed BED directly — faster and Perl-free for crosslink calling.
- **Parallel Mode**: `--run-cims-cits --run-clink` runs both CTK and Clink on the same BAM. Outputs land in numbered folders (`5_CTK_Analysis/`, `6_Clink/`, `7_OTHERS/`) for side-by-side comparison.
- **New `lib/clink/` module**: 5 Python scripts (`collapse.py`, `pileup.py`, `cits.py`, `cims.py`, `stats.py`) sourced from `luna/v3.3`, integrated without disturbing any existing CTK hardening.
- **Install script updates**: `pysam` and `umi_tools` added to conda environments; Python pinned to `<3.13` (umi_tools incompatibility).

## Update Log (v3.2.0)
This release introduces a new preprocessing-only standalone script, a GEO deposit mode, significant runtime improvements via lazy decompression, and cleaner output directory management.

**New Features (v3.2.0):**
- **Hash-Based Deduplication Engine**: Replaced `fastq2collapse.pl` (CTK) with a new Python-based collapse engine (`lib/fastq_collapse_hash.py`). Uses an O(n) hash (Counter) approach with no disk spill, delivering significant speedups on large pooled libraries. Output format is identical to `fastq2collapse.pl` (`@id#COUNT`), ensuring full compatibility with downstream CTK tools. A sort-based fallback (`awk | sort | uniq`) is used automatically if Python 3 is unavailable.
- **Modular Deduplication** (`lib/dedup.sh`): Deduplication logic is now a standalone library module sourced by all scripts. Eliminates redundant gzip round-trips — the dedup step accepts and emits plain `.fastq` by default.
- **PREPittyPrep.sh**: New standalone preprocessing script. Runs the full CLIP-seq preprocessing stack (dedup → demux → fastp) and stops before alignment, producing ready-to-map `*_prepped.fastq.gz` files. Supports single-file, pooled+demux, and directory batch modes.
- **GEO Mode** (`PREPittyPrep.sh --geo`): Raw barcode-based demultiplexing only — reads are split by barcode and written exactly as received (no dedup, no adapter trimming). Intended for GEO/SRA deposit of pooled library FASTQs.
- **Lazy Gzip**: All internal pipeline steps now operate on plain (uncompressed) `.fastq` files. Gzip compression is only applied to final output products, eliminating redundant compress/decompress cycles between steps for significant runtime improvement on large libraries.
- **Conditional `0_DEMUX_FASTQ`**: The demux FASTQ output folder is now only created and retained when the `-k` (`--keep`) flag is passed. Without `-k`, split FASTQs are discarded after child processes complete, keeping the output directory clean.
- **Consistent Help Behavior**: `PEAKittyPeak.sh` now shows help and exits when called with no arguments, matching the behavior of `MAPittyMap.sh` and `CLIPittyClip.sh`.

## Update Log (v3.1.0 / v3.0.2)
This release focuses on heavily expanding CTK feature support, providing cleaner outputs, and patching edge cases in directory and single-file processing.

**New Features (v3.1.0):**
- **Peak Caller Selection**: Added `--peak-caller` flag supporting `ctk` (`tag2peak.pl`) and `homer` (default).
- **Custom Peak Caller Arguments**: Added `--peak-caller-args` pass-through for custom `tag2peak.pl` or `findPeaks` arguments.
- **Track Headers**: Added strict sorting and UCSC track headers to all output BedGraph files for immediate visualization compatibility.
- **Uncompressed FASTQ Support**: Added support for plain `.fastq` and `.fq` files in both single-file (`-i`) and directory modes (`-d`). They are auto-gzipped into temporary file, which will be removed after the pipeline finishes.
- **Organized Outputs**: Upgraded the generic `_analysis` directory in the single-file (`-i`) execution into the same organized, numbered folders used in Demultiplex and Batch modes.
- **Flag Renaming**: Renamed `--run-ctk` to `--run-cims-cits` for explicit clarity. (The deprecated `--run-ctk` alias will still function with a warning).
- **Interactive Wizard Mode Overhaul**: The interactive `--wizard` mode fully supports the selection of peak callers (HOMER vs CTK) and features adaptive prompts based on the chosen caller.

**Bug Fixes (v3.0.2):**
- **Directory Mode Deduplication**: Fixed an issue where `-d` batch mode unintentionally skipped PCR deduplication by hardcoding `--no-dedup`. Children processes now inherit the parent `DEDUP_MODE` properly.
- **Weight Assignments**: Avoided applying duplicate weights formatting if deduplication wasn't run.
- **Peak Aggregation**: CTK peak tag files (flat BED) are now fully supported, bypassing the bug where aggregation expected a HOMER `_peaks` directory layout.
- **Completion Output**: Missing CLI completion messaging and `Done!` status lines in single-file modes have been restored.
- **Standalone Safety**: Added a `/dev/null` fallback for `LOG_FILE` variables to resolve crashes when running `PEAKittyPeak.sh` standalone.



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

> [!WARNING]
> **macOS Compatibility Note:**
> STAR version `2.7.11b` (current latest in Bioconda) is **broken on macOS Tahoe (26.2)** when running via Rosetta (fails to spawn sub-processes).
> **Solution:** You MUST install/pin STAR to version **2.7.10b**: `mamba install bioconda::star=2.7.10b` (last tested: 2026-04-08)

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
python3 -c "import pysam; print('pysam OK')"
umi_tools --version
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

# With CIMS analysis only (--genome-fasta strongly recommended for optimal deletion detection)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-cims

# With CITS analysis only
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-cits

# With both CIMS and CITS analysis
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-cims --run-cits
# OR equivalently:
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-cims-cits

# Clink only (Python-native, no Perl CTK required)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-clink

# CTK + Clink in parallel (side-by-side comparison)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index --genome-fasta /path/to/genome.fa -t 8 --run-cims-cits --run-clink

# eCLIP paired-end analysis (post-eclipdemux R2 files, UMI in header)
CLIPittyClip.sh --eclip pe -d /path/to/samples/ -x /path/to/star_index -t 8 --run-cims --run-cits

# eCLIP single-end analysis (raw seCLIP R1 files, UMI in sequence)
CLIPittyClip.sh --eclip se -i reads.fastq.gz -x /path/to/star_index -t 8
```

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
| — | `--align-mismatches` | int | 2 | Absolute mismatch backstop (STAR only; primary filter is fractional 10% of read length) |
| — | `--genome-fasta` | path | — | Reference FASTA; enables `samtools calmd` for accurate MD tags (strongly recommended for CIMS/CITS) |
| `-v` | `--verbose` | bool | false | Enable verbose logging |
| `-h` | `--help` | — | — | Show help message |

### Preprocessing Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-u` | `--umi-length` | int | — | UMI length (auto-detected for eCLIP) |
| `-a` | `--adapter` | string | L32 | 3' adapter sequence |
| — | `--eclip` | string | — | eCLIP mode: `pe` (paired-end) or `se` (single-end) |
| — | `--no-dedup` | bool | — | Disable FASTQ deduplication (default: ON) |
| — | `--filter-ncrna` | bool | false | Enable ncRNA pre-filtering (opt-in) |

## Deduplication

CLIPittyClip uses a custom hash-based deduplication engine (`lib/fastq_collapse_hash.py`) as a drop-in replacement for CTK's `fastq2collapse.pl`.

**Engine selection (automatic):**

| Priority | Engine | Algorithm | Notes |
|----------|--------|-----------|-------|
| 1st | `fastq_collapse_hash.py` | O(n) hash (Python `Counter`) | No disk spill; requires Python 3 |
| 2nd | `awk \| sort \| uniq` | O(n log n) sort-based | Automatic fallback; mirrors `fastq2collapse.pl` exactly |

**Output format:** `@read_id#COUNT` — identical to `fastq2collapse.pl`, fully compatible with all downstream CTK tools (`tag2collapse.pl`, `parseAlignment.pl`, etc.).

**To disable:** Pass `--no-dedup` to skip deduplication entirely (e.g., for pre-deduplicated inputs or eCLIP data already processed upstream).

### Demultiplexing Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-b` | `--barcodes` | path | — | Barcode file (enables demultiplexing) |
| — | `--demux-mismatches` | int | 1 | Max barcode mismatches |

### Peak Calling Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| — | `--peak-caller` | str | homer | Peak caller: `homer` or `ctk` (tag2peak.pl) |
| — | `--peak-caller-args` | str | — | Additional peak caller arguments (quoted string) |

### CTK Analysis Options

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| — | `--run-cims-cits` | bool | — | Enable full CTK CIMS+CITS analysis |
| — | `--run-cims` | bool | — | Enable CIMS analysis (mutation sites) |
| — | `--run-cits` | bool | — | Enable CITS analysis (truncation sites) |
| — | `--cims-iter` | int | 5 | CIMS permutation iterations |
| — | `--cims-fdr` | float | 0.05 | CIMS FDR threshold |
| — | `--cits-pval` | float | 0.05 | CITS p-value threshold |
| — | `--cits-gap` | int | 25 | CITS clustering gap (-1 disables) |
| `-f` | `--flank` | int | 10 | Flanked BED nucleotides |
| — | `--no-motif` | bool | — | Skip flanked BED generation |

> **Note:** `--run-ctk` is a deprecated alias for `--run-cims-cits` (still works, emits a warning).

### Clink Analysis Options

Clink is a Python-native crosslink site caller (`umi_tools` + `pysam`). It runs independently of CTK and can be used alone or in parallel.

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| — | `--run-clink` | bool | — | Enable Clink pipeline (umi_tools dedup + Python pileup → CITS/CIMS) |
| — | `--clink-umi-len` | int | -1 | UMI length for umi_tools (-1 = auto-detect) |
| — | `--clink-fdr` | float | 0.05 | Benjamini-Hochberg FDR threshold |
| — | `--clink-min-cov` | int | 5 | Minimum coverage to test a position |

**Requires:** `pysam`, `numpy`, `scipy`, `umi_tools` (all installed by `install_macos.sh` / `install_linux.sh`).

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

## CIMS/CITS Mapping Parameters

CLIPittyClip tunes STAR alignment specifically for CIMS (crosslink-induced mutation sites) and CITS (crosslink-induced truncation sites) analysis, where reads carrying genuine crosslink deletions must survive alignment.

### Why `--genome-fasta` matters

STAR index directories contain only precomputed binary files — the source FASTA is consumed at index-build time and not stored there. Without the reference FASTA, `samtools calmd` cannot run, and STAR's native MD tags at deletion boundaries (especially in homopolymer runs) may be inconsistent. `parseAlignment.pl` uses both the CIGAR string and the MD tag to classify mutations; inconsistent MD tags → missed or misclassified crosslink deletions.

Provide `--genome-fasta /path/to/genome.fa` to enable `samtools calmd` and recalculate MD tags authoritatively before parsing.

### STAR gap and mismatch parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--outFilterMismatchNoverReadLmax` | `0.1` | Fractional filter: allows ~3 mismatches in a 30 bp read, ~2 in a 20 bp read. Scales with post-trim read length — equivalent to BWA's `-n 0.06` philosophy. Replaces the old absolute `--outFilterMismatchNmax 2`, which discarded reads carrying a genuine crosslink deletion + 2 sequencing errors (NM = 3). |
| `--outFilterMismatchNmax` | `5` | Hard backstop only; blocks extreme cases but is not the primary filter. |
| `--scoreDelOpen` / `--scoreDelBase` | `-1` | Reduces deletion open/extension penalty to match the substitution cost. STAR's default (`-2/-2`) makes a 1-nt deletion score worse than a substitution on short (20–30 bp) iCLIP reads, causing crosslink-induced deletions to be realigned as mismatches. |
| `--scoreInsOpen` / `--scoreInsBase` | `-1` | Symmetrically reduces insertion penalty (CTK CIMS processes all indels). |

> **Bowtie2 note:** Bowtie2 is not splice-aware and its gap penalties are not tuned for CIMS/CITS. It is supported for speed but STAR is strongly recommended when running CIMS or CITS analysis. A warning is emitted automatically.

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

### CLIPittyClip.sh
```
{INPUT}_output/
├── 0_DEMUX_FASTQ/           # Demultiplexed reads  ← only kept with -k
├── 1_BAM/                   # Aligned BAM files
├── 2_COLLAPSED_BED/         # PCR-deduplicated BED
├── 3_BEDGRAPH/              # Coverage tracks (Normalized Coverage)
│   ├── COMBINED_BEDGRAPH/   # (Group mode) Averaged replicates
├── 4_PEAKS/                 # HOMER peak results
│   ├── COMBINED_PEAKS/      # Peaks across all samples (or aggregated)
│   └── SAMPLE_PEAKS/        # Peaks per sample
│
├── 5_CTK_Analysis/          # When --run-cims-cits  (or --run-cims / --run-cits)
│   ├── CIMS/
│   ├── CITS/
│   └── motif_analysis/
│ OR
├── 5_CIMS_Analysis/         # When --run-cims only
│ OR
├── 5_CITS_Analysis/         # When --run-cits only
│ OR
├── 5_Clink/                 # When --run-clink only (no CTK flags)
│
├── 6_Clink/                 # When both CTK + Clink (--run-cims-cits --run-clink)
│
├── 5_OTHERS/ or 6_OTHERS/ or 7_OTHERS/   # Support files (number auto-adjusts)
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

> **Output folder numbering:**
>
> | Flags | Folder layout |
> |-------|---------------|
> | *(none)* | `5_OTHERS/` |
> | `--run-cims-cits` | `5_CTK_Analysis/`, `6_OTHERS/` |
> | `--run-clink` | `5_Clink/`, `6_OTHERS/` |
> | `--run-cims-cits --run-clink` | `5_CTK_Analysis/`, `6_Clink/`, `7_OTHERS/` |


## Peak Matrix Metrics

The `COMBINED_PEAK_MATRIX.txt` file contains up to 54+ advanced metrics summarizing peak intensity, reproducibility, and coverage shape.

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

## eCLIP Analysis Modes

CLIPittyClip supports two eCLIP protocols, selected via `--eclip <pe|se>`.

---

### eCLIP Paired-End Mode (`--eclip pe`)

For ENCODE eCLIP paired-end data after inline-barcode demultiplexing by `eclipdemux`.

```bash
# PE eCLIP — batch (typical for ENCODE downloads)
CLIPittyClip.sh --eclip pe -d /path/to/eclip_fastqs/ -x /path/to/star_index -t 8

# PE eCLIP with CIMS/CITS analysis
CLIPittyClip.sh --eclip pe -d /path/to/eclip_fastqs/ -x /path/to/star_index -t 8 --run-cims --run-cits
```

**Expected input:** Read 2, post-eclipdemux format — UMI already moved to read header by `eclipdemux`.
Read header format: `@NTACGTTGAT:NB501168:530:HJ3WMBGXF:...`

**Preprocessing workflow:**
1. **Validate input** — confirms R2 with UMI in header (colon-prefix format); exits with a clear error if wrong file is supplied
2. **UMI to sequence** — prepends the UMI back onto the read sequence for deduplication
3. **Deduplicating** — exact-duplicate collapse using hash-based engine (UMI-aware)
4. **Extract UMI** — strips UMI from sequence into CTK-compatible header format (`@READ#count#UMI`)
5. **Adapter Trim** — fastp with full eCLIP adapter set (`eclip_adapters.fa`: inline barcodes + TruSeq R2)

**Ignored options:** `-u` (UMI length is auto-detected from header), `-a` (adapter is hardcoded to `eclip_adapters.fa`)

> **Note:** `--eclip pe` only affects preprocessing. CIMS/CITS analysis requires separate `--run-cims`/`--run-cits` flags.

---

### eCLIP Single-End Mode (`--eclip se`)

For single-end eCLIP (seCLIP) data following the Blue et al. 2022 protocol.

```bash
# SE eCLIP — single file
CLIPittyClip.sh --eclip se -i sample_R1.fastq.gz -x /path/to/star_index -t 8

# SE eCLIP — batch directory
CLIPittyClip.sh --eclip se -d /path/to/seclip_fastqs/ -x /path/to/star_index -t 8
```

**Expected input:** Raw Read 1, unprocessed — UMI is the first 10nt of the read sequence.
This is the fastq directly from ENCODE before any barcode/UMI processing.

**Hardcoded parameters (Blue et al. 2022 — not user-configurable):**
- UMI length: 10nt
- Adapter: TruSeq Read 1 (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`)

**Preprocessing workflow:**
1. **Validate input** — confirms R1 with UMI in sequence; exits with a clear error if a pre-processed or R2 file is supplied
2. **Deduplicating** — exact-duplicate collapse on full read (UMI + cDNA sequence together)
3. **Extract UMI** — strips first 10nt into CTK-compatible header format (`@READ#count#UMI`)
4. **Adapter Trim** — fastp with TruSeq R1 adapter

**Ignored options:** `-u` (UMI length hardcoded to 10nt), `-a` (TruSeq R1 is hardcoded)

> **Note:** `--eclip se` uses the same downstream alignment, PCR collapse, and peak-calling steps as all other modes. Only preprocessing differs.

## Standalone Tools

### PREPittyPrep.sh
Standalone preprocessing module. Runs dedup → demux (optional) → fastp and stops before alignment, producing clean, ready-to-map `*_prepped.fastq.gz` files. No genome index required.

**Output Structure:**
```text
{INPUT}_prepped/
├── PREPPED_FASTQ/           # Ready-to-map gzipped FASTQs (*_prepped.fastq.gz)
├── 0_DEMUX_FASTQ/           # Split FASTQs  ← only kept with -k (demux mode)
└── REPORTS/
    ├── FASTP_REPORT/        # HTML/JSON QC reports per sample
    └── detailed_output.log
```

**GEO Deposit Output Structure:**
```text
{INPUT}_GEO/
├── sample1.fastq.gz         # Reads split by barcode, unmodified
├── sample2.fastq.gz
├── md5sums.txt              # Standard MD5 checksums for all FASTQs (for GEO upload)
├── ...
└── unknown.fastq.gz         # Unmatched reads
```

**Input modes:**

| Mode | Flags | Behavior |
|------|-------|----------|
| Single file | `-i reads.fastq.gz` | dedup → fastp → `_prepped.fastq.gz` |
| Pooled + demux | `-i pool.fastq.gz -b barcodes.txt` | dedup → cutadapt split → per-sample fastp |
| Directory batch | `-d /path/to/samples/` | per-file: dedup → fastp |
| GEO deposit | `--geo -i pool.fastq.gz -b barcodes.txt` | raw cutadapt split only, no modification. md5 checksum values output as a text file |

**Key options:**
- `-u <int>`: UMI length for extraction (default: 0)
- `-b <path>`: Barcode file — enables demux mode
- `-a <str>`: 3' adapter sequence
- `--bc-len <int>`: Barcode length (auto-detected from `-b`)
- `--spacer-len <int>`: Spacer length after barcode (default: 0)
- `--no-dedup`: Skip deduplication
- `--demux-mismatches <N>`: Barcode mismatch tolerance (default: 1)
- `--filter-ncrna`: Enable ncRNA pre-filtering (requires `-x`)
- `-k`: Keep intermediate files including `0_DEMUX_FASTQ/`
- `--geo`: GEO mode — raw demux only, reads not modified

```bash
# Preprocess a single FASTQ
PREPittyPrep.sh -i reads.fastq.gz -u 7 -t 8

# Demux + preprocess a pooled library
PREPittyPrep.sh -i pool.fastq.gz -b barcodes.txt -u 7 -t 8

# Batch preprocess all FASTQs in a directory
PREPittyPrep.sh -d /path/to/samples/ -u 7 -t 8

# GEO deposit: raw barcode split, no modification
PREPittyPrep.sh -i pool.fastq.gz -b barcodes.txt -u 7 --geo -o my_GEO
```

> **Note:** `PREPittyPrep.sh` output is fully compatible with `MAPittyMap.sh` and `CLIPittyClip.sh -d` for downstream alignment.

---

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
Standalone peak calling module. Supports HOMER (default) and CTK `tag2peak.pl`. Requires a directory containing collapsed BED files.

**Required inputs:**
- Run from a directory containing a `BED/` folder with `.bed` files
- OR use `-i <directory>` to specify input BED directory explicitly

**Key options:**
- `--peak-caller <str>`: Peak caller: `homer` (default) or `ctk`
- `-p <int>`: Min distance between peaks (default: 50)
- `-z <int>`: Peak size (default: 20, HOMER only)
- `-f <int>`: Fragment length (default: 25, HOMER only)
- `-n <string>`: Output name prefix
- `-a <string>`: Additional HOMER findPeaks arguments (HOMER only)
- `--aggregate`: Combine all input BED files into a single meta-sample for peak calling
- `--no-aggregate`: Process each BED file individually (Default)
- `--ctk-dir <path>`: Add CIMS/CITS site counts from CTK analysis
- `--ctk-group <file>`: Groups file for CTK aggregation
- `--wizard`: Launch interactive HOMER configuration wizard

```bash
# Individual Peak Calling (Default, HOMER)
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Analysis --no-aggregate

# Aggregated Peak Calling using CTK tag2peak.pl
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Combined --aggregate --peak-caller ctk

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

**To enable ncRNA filtering:** Pass `--filter-ncrna` to opt in. Omitting the flag leaves all reads unfiltered (the default).

> **Note:** ncRNA filtering in `PREPittyPrep.sh` also requires `-x` (genome index path) alongside `--filter-ncrna`.

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
