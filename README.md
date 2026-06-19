<p align="center">
  <img src="logo.png" alt="CLIPittyClip Logo" width="400">
</p>

# CLIPittyClip: Modern CLIP-seq Analysis Pipeline
**Version 3.5.0**

A comprehensive, single-command CLIP-seq analysis pipeline from raw FASTQ to peaks and crosslink sites. Supports iCLIP, irCLIP, eCLIP, PAR-CLIP, and related variant protocols.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Input Modes](#input-modes)
5. [Output Structure](#output-structure)
6. [Command-Line Reference](#command-line-reference)
7. [Crosslink Site Analysis: CTK and Clink](#crosslink-site-analysis-ctk-and-clink)
8. [Protocol-Specific Modes](#protocol-specific-modes)
9. [Standalone Tools](#standalone-tools)
10. [Genome Index Setup](#genome-index-setup)
11. [Peak Matrix Metrics](#peak-matrix-metrics)
12. [Changelog](#changelog)

---

## Overview

<p align="center">
  <img src="flowchart.png" alt="CLIPittyClip Pipeline Flow" width="800">
</p>

CLIPittyClip runs the complete CLIP-seq stack in a single command. All three stages run automatically:

1. **Preprocessing** — global sequence dedup, adapter trimming (fastp), optional repeat element pre-filter (rRNA/tRNA/TEs), genome alignment (STAR / Bowtie2)
2. **Crosslink sites** — two parallel tracks from the same BAM:
   - **Clink** (`--run-clink`, v3.3): `umi_tools` dedup → `pileup.py` → strand-aware CITS + CIMS (Python, BH FDR)
   - **CTK** (`--run-cims-cits`): `tag2collapse.pl` → `parseAlignment.pl` → `CITS.pl` + `CIMS.pl` (Perl, permutation FDR)
3. **Peaks** — RPM bedgraphs, HOMER or CTK peak calling, consolidated peak matrix, and grouped averaging for downstream analysis and visualization

---

## Installation

> [!WARNING]
> **macOS:** STAR `2.7.11b` is broken on macOS Tahoe via Rosetta. Pin to `2.7.10b`:
> `mamba install bioconda::star=2.7.10b`

### 1. Clone

```bash
git clone https://github.com/LunaRNALab/CLIPittyClip.git
cd CLIPittyClip
```

### 2. Install

**macOS (Intel or Apple Silicon):**
```bash
./install_macos.sh --env clipittyclip --tools-dir ~/Tools
```

**Linux:**
```bash
./install_linux.sh --env clipittyclip --tools-dir ~/Tools
```

The scripts create a `clipittyclip` conda environment with all dependencies: STAR, Bowtie2, samtools, fastp, bedtools, ucsc-bedgraphtobigwig, CTK, HOMER, Perl modules (Math::CDF, Bio::SeqIO), and Python packages (pysam, numpy, scipy, umi_tools).

| Option | Default | Description |
|--------|---------|-------------|
| `--env <name>` | `clipittyclip` | Conda environment name |
| `--tools-dir <path>` | `~/Tools` | Directory for CTK and HOMER |

### 3. Activate and verify

```bash
source ~/.zshrc          # or ~/.bashrc on Linux
conda activate clipittyclip

# Verify
which CLIPittyClip.sh
which parseAlignment.pl
which findPeaks
perl -MMath::CDF -e 'print "Math::CDF OK\n"'
python3 -c "import pysam; print('pysam OK')"
umi_tools --version
```

---

## Quick Start

```bash
# Single sample — STAR alignment
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8

# With UMI length
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index -t 8 -u 7

# Pooled library — demultiplex first, then analyze each sample
CLIPittyClip.sh -i pool.fastq.gz -b barcodes.txt -x /path/to/star_index -t 8

# Pre-demultiplexed folder — batch mode
CLIPittyClip.sh -d /path/to/samples/ -x /path/to/star_index -t 8

# Crosslink sites with Clink (recommended, v3.5)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index \
    --genome-fasta /path/to/genome.fa -t 8 --run-clink

# Head-to-head: CTK vs Clink on the same data
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index \
    --genome-fasta /path/to/genome.fa -t 8 --run-cims-cits --run-clink

# eCLIP paired-end
CLIPittyClip.sh --eclip pe -d /path/to/eclip_r2s/ -x /path/to/star_index -t 8 --run-clink

# PAR-CLIP (T>C crosslink sites)
CLIPittyClip.sh --parclip -i reads.fastq.gz -u 4 -x /path/to/star_index -t 8 --run-clink

# Interactive wizard (guided setup; emits the equivalent CLI for next time)
CLIPittyClip.sh -w
```

---

## Input Modes

| Mode | Flags | Use Case |
|------|-------|----------|
| Single file | `-i sample.fastq.gz` | One FASTQ, direct analysis |
| Pooled + barcodes | `-i pool.fastq.gz -b barcodes.txt` | Demultiplex then analyze each sample |
| Pre-demuxed folder | `-d /path/to/folder/` | Batch-analyze a set of FASTQs |

CLIPittyClip accepts both `.fastq.gz` and plain `.fastq` / `.fq` inputs in all modes.

---

## Output Structure

All results land in a single numbered-folder hierarchy next to your input (or at `-o`).

```
{INPUT}_output/
├── 00_REPORTS/
│   ├── FASTP_REPORT/        ← HTML/JSON QC
│   ├── ALIGNER_LOGS/        ← STAR/Bowtie2 summaries
│   ├── PEAK/                ← peak calling logs
│   └── SAMPLES/             ← per-sample detailed logs
├── 0_DEMUX_FASTQ/           ← demultiplexed reads (only with -k)
├── 01_BAM/                  ← sorted, indexed BAM files
├── 02_COLLAPSED_BED/        ← PCR-deduplicated read BED
├── 03_COVERAGE/             ← RPM-normalized ± strand bedgraphs
│   ├── COMBINED_BEDGRAPH/  ← group-averaged tracks (with -g)
│   └── BigWig/              ← per-sample crosslink-site bigWigs (with --xl-bigwig)
│       ├── {sample}_xl_pos.bw  ← + strand crosslink counts
│       └── {sample}_xl_neg.bw  ← - strand crosslink counts
├── 04_PEAKS/
│   ├── SAMPLE_PEAKS/        ← per-sample peak calls
│   └── COMBINED_PEAKS/      ← aggregated peaks + COMBINED_PEAK_MATRIX.txt
│
├── 05_CTK_Analysis/         ← CTK crosslink sites (--run-cims-cits)
│   ├── {sample}/CIMS/ CITS/
│   └── GROUP_{group}/CIMS/ CITS/  ← group-pooled results (--ctk-group)
│
├── 05_Clink/ or 06_Clink/  ← Clink crosslink sites (--run-clink)
│   ├── {sample}/
│   │   ├── {sample}_dedup.bam
│   │   ├── {sample}_pileup.npz
│   │   ├── {sample}_all_crosslinks.bed  ← all positions with ≥1 truncation (unfiltered)
│   │   ├── {sample}_truncations.bed     ← CITS: truncation sites passing BH FDR cutoff (--clink-fdr)
│   │   ├── {sample}_deletions.bed
│   │   └── {sample}_TtoC.bed  (+ all 12 substitution types)
│   └── GROUP_{group}/               ← group-pooled results (--group-xlsite)
│       ├── GROUP_{group}_pileup.npz
│       ├── GROUP_{group}_all_crosslinks.bed
│       ├── GROUP_{group}_truncations.bed
│       └── GROUP_{group}_deletions.bed  (+ substitution types)
│
└── 05_OTHERS/ → 06_OTHERS/ → 07_OTHERS/  ← number adjusts automatically; always last
    ├── STAR_OUTPUT/
    └── Repeat_Mapping/      ← repeat element BAMs, stats, and quantification TSVs (if --filter-repeat)
```

> **Folder numbering:** `05_OTHERS` with no crosslink analysis; `06_OTHERS` with CTK or Clink; `07_OTHERS` with both. CTK is always `05`, Clink is `05` (no CTK) or `06` (with CTK).

**Output location (`-o`):**
- No `-o`: created next to input (`/data/reads.fq.gz` → `/data/reads_output/`)
- Name only (`-o HepG2`): created next to input as `/data/HepG2/`
- Full path (`-o /results/HepG2`): exact path used

---

## Command-Line Reference

Run `CLIPittyClip.sh --help` for full usage.

### Input / Output

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-i` | `--input-file` | — | Input FASTQ (required unless using `-d`) |
| `-d` | `--input-dir` | — | Directory of pre-demultiplexed FASTQs |
| `-x` | `--index` | — | Genome index directory (required) |
| `-o` | `--output` | next to input | Output folder name or full path |
| `-k` | `--keep` | off | Keep intermediate files |

### Alignment

| Long | Default | Description |
|------|---------|-------------|
| `-m` / `--mapper` | `star` | Aligner: `star` (default) or `bowtie2` |
| `-t` / `--threads` | `1` | Number of threads |
| `--genome-fasta` | — | Reference FASTA — enables `samtools calmd` for accurate MD tags; strongly recommended for crosslink site analysis |
| `--align-mismatches` | `2` | Absolute mismatch backstop (STAR; primary filter is fractional 10% of read length) |

### Preprocessing

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-u` | `--umi-length` | auto | UMI length in bases |
| `-a` | `--adapter` | L32 | 3' adapter sequence |
| `-b` | `--barcodes` | — | Barcode file (enables demultiplexing) |
| — | `--demux-mismatches` | `1` | Max barcode mismatches |
| — | `--eclip` | — | eCLIP mode: `pe` (paired-end) or `se` (single-end) |
| — | `--parclip` | — | PAR-CLIP mode: specialized preprocessing for 4SU CLIP (requires `-u`) |
| — | `--parclip-adapters` | bundled | Custom PAR-CLIP adapter FASTA (default: `lib/parclip_adapters.fa`) |
| — | `--no-dedup` | — | Skip FASTQ deduplication |
| — | `--filter-repeat` | off | Pre-filter repeat element reads: rRNA, tRNA, and transposable elements (opt-in) |
| — | `--bc-len` | — | Barcode length (auto-detected from `-b`) |
| — | `--spacer-len` | `0` | Spacer bases after barcode |

### Peak Calling

| Long | Default | Description |
|------|---------|-------------|
| `--peak-caller` | `homer` | Peak caller: `homer` or `ctk` |
| `--peak-caller-args` | — | Extra arguments passed to peak caller (quoted string) |
| `-f` / `--flank` | `10` | Flanking nucleotides for motif BED |
| `--no-motif` | — | Skip flanked BED generation |

### Grouping

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-g` | `--groups` | — | Groups file for bedgraph/peak aggregation (`SampleName\tGroupName`) |
| — | `--ctk-group` | off | Pool samples by group before running CTK crosslink analysis |
| — | `--group-xlsite` | off | Pool samples by group for Clink CITS/CIMS (and CTK if enabled). Per-sample dedup BAMs are produced first, then merged by group for pileup → CITS/CIMS. Requires `-g`. |

### Crosslink Site Analysis

#### CTK (standard)

| Long | Default | Description |
|------|---------|-------------|
| `--run-cims-cits` | off | Enable full CTK CIMS + CITS |
| `--run-cims` | off | CIMS only (mutation/deletion sites) |
| `--run-cits` | off | CITS only (truncation sites) |
| `--cims-iter` | `5` | CIMS permutation iterations |
| `--cims-fdr` | `0.05` | CIMS FDR threshold |
| `--cits-pval` | `0.05` | CITS p-value threshold |
| `--cits-gap` | `25` | CITS clustering gap (`-1` disables) |

#### Clink (v3.5)

| Long | Default | Description |
|------|---------|-------------|
| `--run-clink` | off | Enable Clink crosslink site analysis |
| `--clink-umi-len` | auto | UMI length for umi_tools (auto-detected if omitted) |
| `--clink-fdr` | `0.05` | Benjamini-Hochberg FDR threshold |
| `--clink-min-cov` | `5` | Minimum coverage to test a position |
| `--clink-multi-map` | off | Rescue multi-mapped reads (NH:i:>1) via single-pass positional assignment before deduplication. Requires `pysam`. See note below. |
| `--xl-bigwig` | off | Generate per-sample strand-specific crosslink-site bigWig files in `03_COVERAGE/BigWig/`. Each file records per-nucleotide truncation event counts (read 5′ ends), **not** RPM read coverage — suitable for BindingSiteFinder and similar tools. Requires `bedGraphToBigWig` in PATH (`ucsc-bedgraphtobigwig` conda package). |

> **Grouped Clink analysis:** combine `--run-clink --group-xlsite -g groups.txt` to produce per-sample dedup BAMs first, then merge by group for pileup → CITS/CIMS. Group results land in `GROUP_<name>/` inside the Clink output directory.

> **Multi-mapper rescue (`--clink-multi-map`):** unique reads (NH:i:1) are deduped first; multi-mapped reads (NH:i:>1) are hard-assigned to the candidate locus with the greatest unique-read pileup support within a ±50 bp window, then passed through normal UMI deduplication, and merged with the unique dedup BAM. Compatible with `--group-xlsite`. Off by default — all analyses presented here use unique reads only.

### Other

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-s` | `--sample` | — | Test mode: process only first N reads |
| `-w` | `--wizard` | — | Launch interactive configuration wizard |
| — | `--notification` | off | System notification on completion |
| `-v` | `--verbose` | off | Verbose logging |
| `-h` | `--help` | — | Show help |

---

## Crosslink Site Analysis: CTK and Clink

CLIPittyClip supports two crosslink site callers that can run in parallel on the same BAM.

**CTK** (`--run-cims-cits`) uses Perl tools from the Chaolin Zhang lab: `tag2collapse.pl` deduplication, `parseAlignment.pl` signal extraction, and permutation-based FDR.

**Clink** (`--run-clink`, v3.5) is CLIPittyClip's Python-native caller: `umi_tools` deduplication, a single-pass `pileup.py` scan shared by both CITS and CIMS, and Benjamini-Hochberg FDR.

```bash
# CTK only
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index \
    --genome-fasta /path/to/genome.fa -t 8 --run-cims-cits

# Clink only (recommended)
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index \
    --genome-fasta /path/to/genome.fa -t 8 --run-clink

# Both in one run
CLIPittyClip.sh -i reads.fastq.gz -x /path/to/star_index \
    --genome-fasta /path/to/genome.fa -t 8 --run-cims-cits --run-clink
```

**Clink output** (`5_Clink/{sample}/`): `_dedup.bam`, `_pileup.npz`, `_all_crosslinks.bed` (every position with ≥1 truncation, unfiltered — suitable for PEKA and BindingSiteFinder), `_truncations.bed` (FDR-filtered CITS), `_deletions.bed`, `_TtoC.bed` (+ all 12 substitution types).

> [!NOTE]
> `--genome-fasta` is strongly recommended for crosslink site analysis. STAR index directories don't store the source FASTA; without it, `samtools calmd` can't recalculate MD tags and crosslink deletions near homopolymers may be missed.

### Group-based crosslink analysis

Pool replicates before running CITS/CIMS to increase depth:

```bash
# CTK
CLIPittyClip.sh -i pool.fq.gz -b barcodes.txt -x index \
    --run-cims-cits -g groups.txt --ctk-group

# Clink
CLIPittyClip.sh -d /path/to/samples/ -x index \
    --run-clink -g groups.txt --group-xlsite
```

**groups.txt** (tab-separated): `SampleName\tGroupName`. Samples not listed run individually.

> [!WARNING]
> Group-based CTK aggregates all samples before CIMS/CITS. This can require >64 GB RAM for large datasets.

---

## Protocol-Specific Modes

### eCLIP Paired-end (`--eclip pe`)

For ENCODE eCLIP data after `eclipdemux`. Supply **Read 2** with UMI already in the header (`@NTACGTTGAT:NB501168:...`).

```bash
CLIPittyClip.sh --eclip pe -d /path/to/eclip_r2s/ -x /path/to/star_index -t 8 --run-clink
```

Preprocessing: validate R2 format → UMI to sequence → hash dedup → extract UMI → fastp (full eCLIP adapter set). UMI length auto-detected from header.

### eCLIP Single-end (`--eclip se`)

For seCLIP (Blue et al. 2022). Supply raw **Read 1** — UMI is the first 10 nt of the sequence. UMI length and adapter hardcoded (10 nt; TruSeq R1).

```bash
CLIPittyClip.sh --eclip se -i sample_R1.fastq.gz -x /path/to/star_index -t 8
```

### PAR-CLIP (`--parclip`)

For 4-thiouridine CLIP data. Read layout: `[UMI][READ][2nt spacer][6mer barcode][adapter]`. Requires `-u` (UMI length). Primary crosslink signal is T→C substitution — use `--run-clink` for CIMS.

```bash
CLIPittyClip.sh --parclip -i reads.fastq.gz -u 4 \
    -x /path/to/star_index -t 8 --run-clink
```

Preprocessing: hash dedup → fastp pass 1 (UMI extraction + barcode+adapter trim using `lib/parclip_adapters.fa`) → fastp pass 2 (2nt spacer trim). A custom adapter FASTA can be supplied with `--parclip-adapters`.

The bundled `lib/parclip_adapters.fa` contains 43 MDX-O barcodes (6mer + `TGGAATTCTCGGGTGCCAAGG`). The leading 2nt random spacer is not included in the FASTA — it is removed in pass 2 after the barcode+adapter is trimmed.

---

## Standalone Tools

### PREPittyPrep.sh

Preprocessing only — dedup → demux (optional) → fastp → ready-to-map FASTQs. No genome index required.

```bash
# Single FASTQ
PREPittyPrep.sh -i reads.fastq.gz -u 7 -t 8

# Pooled library + demux
PREPittyPrep.sh -i pool.fastq.gz -b barcodes.txt -u 7 -t 8

# Batch directory
PREPittyPrep.sh -d /path/to/samples/ -u 7 -t 8

# GEO deposit: raw barcode split, no modification, MD5 checksums
PREPittyPrep.sh -i pool.fastq.gz -b barcodes.txt --geo -o my_GEO

# Interactive wizard
PREPittyPrep.sh -w
```

Output: `{INPUT}_prepped/PREPPED_FASTQ/*_prepped.fastq.gz` + `REPORTS/`

GEO output: `{INPUT}_GEO/{sample}.fastq.gz` + `md5sums.txt`

Key options: `-u` (UMI length), `-b` (barcodes), `-a` (adapter), `--bc-len`, `--spacer-len`, `--no-dedup`, `--geo`, `--filter-repeat`, `-k`, `-w`

---

### MAPittyMap.sh

Standalone alignment + crosslink-site analysis — FASTQ in; BAM + collapsed BED + bedgraph by default, plus optional CTK CIMS/CITS and Clink tracks.

Output tree (`{BASENAME}_mapping/`):

```
1_BAM/                # sorted, indexed BAM
2_BED/                # PCR-collapsed BED (feeds PEAKittyPeak directly)
3_COVERAGE/           # strand-specific normalized bedgraphs
4_CTK_Analysis/       # only with --run-cims / --run-cits
5_Clink_Analysis/     # only with --run-clink
REPORTS/              # alignment logs
```

```bash
# Default: BAM + collapsed BED + bedgraph
MAPittyMap.sh -i reads.fastq.gz -x /path/to/star_index -t 8

# Legacy BAM-only behavior
MAPittyMap.sh -i reads.fastq.gz -x /path/to/star_index --bam-only

# Add CTK CIMS + CITS crosslink-site analysis (genome FASTA recommended)
MAPittyMap.sh -i reads.fastq.gz -x star_index/ --run-cims-cits --genome-fasta ref.fa -t 8

# Add Clink track (Python pileup → CITS/CIMS)
MAPittyMap.sh -i reads.fastq.gz -x star_index/ --run-clink -u 7 -t 8

# Interactive wizard
MAPittyMap.sh -i reads.fastq.gz -x star_index/ -w
```

Key options: `-i`, `-x` (required), `-t`, `-m`, `-o`, `-u`, `--aligner`, `--genome-fasta`, `--filter-repeat`, `--no-chr-filter`, `--bam-only`, `--run-cims`, `--run-cits`, `--run-cims-cits`, `--run-clink`, `--no-motif`, `-w`

The collapsed BED in `2_BED/` is the canonical input to PEAKittyPeak; the CTK output dir is the right argument to PEAKittyPeak's `--ctk-dir` for annotating peaks with CIMS/CITS counts.

---

### PEAKittyPeak.sh

Standalone peak calling from a directory of collapsed BED files.

```bash
# Aggregated peaks (HOMER)
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Combined --aggregate

# Aggregated peaks (CTK tag2peak.pl)
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Combined --aggregate --peak-caller ctk

# With crosslink site counts added to matrix
PEAKittyPeak.sh -i ./2_COLLAPSED_BED -n Combined --aggregate \
    --ctk-dir ./5_CTK_Analysis/

# Interactive wizard
PEAKittyPeak.sh --wizard
```

Key options: `-i` (BED dir), `-n` (output prefix), `--aggregate` / `--no-aggregate`, `--peak-caller`, `--peak-caller-args`, `--ctk-dir`, `--ctk-group`, `-p` (min peak distance), `-z` (peak size), `-f` (fragment length)

---

## Genome Index Setup

### STAR index
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 100
```

### Bowtie2 index
```bash
bowtie2-build genome.fa /path/to/bt2_index/GRCh38
```

### Repeat element pre-filtering index (optional)

Enable with `--filter-repeat`. Reads mapping to the repeat index are diverted to `OTHERS/Repeat_Mapping/` and excluded from genome alignment; unmapped reads continue to STAR normally.

**Index composition (human GRCh38):**

| Source | Accession / URL | Contents |
|--------|----------------|----------|
| NCBI Nucleotide | U13369.1, NR_023363.1 | 45S pre-rRNA (18S + 5.8S + 28S + ITS) and standalone 5S rDNA |
| GtRNAdb | `gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/` | All human tRNA gene sequences (hg38) |
| Dfam | `dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz` | TE consensus sequences, filtered to mammal-relevant clades |

**Building the index (human GRCh38):**

```bash
# 1. rRNA: 45S rDNA unit (18S + 5.8S + 28S) and standalone 5S
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U13369.1&rettype=fasta&retmode=text" \
  > human_45S_rDNA.fa
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_023363.1&rettype=fasta&retmode=text" \
  >> human_45S_rDNA.fa

# 2. tRNA from GtRNAdb
curl -O "https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa"

# 3. Dfam RepeatMasker library
curl -O "https://dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz"
gunzip Dfam-RepeatMasker.lib.gz

# Inspect headers — format is >Name#Class/Family @Clade [S:scores]
grep "^>" Dfam-RepeatMasker.lib | head -10
grep -c "^>" Dfam-RepeatMasker.lib

# 4. Extract only human-relevant clades from Dfam
awk '/^>/{
    keep=0
    if(/(@Vertebrata_vertebrates|@Amniota|@Mammalia|@Theria_mammals|@Eutheria|@Boreoeutheria|@Afrotheria|@Primates|@Haplorrhini|@Catarrhini|@Platyrrhini|@Hominidae|@Homo_sapiens)/) keep=1
} keep{print}' Dfam-RepeatMasker.lib > Dfam-human.fa

# 5. Combine and build the Bowtie2 index
mkdir -p Repeat
cat human_45S_rDNA.fa hg38-tRNAs.fa Dfam-human.fa > Repeat/repeat_sequences.fa
grep -c "^>" Repeat/repeat_sequences.fa   # sanity check total entries
bowtie2-build --threads 8 Repeat/repeat_sequences.fa Repeat/repeat
```

**Building the index (mouse GRCm39):**

Same steps, substituting:
- Mouse 45S rDNA: fetch `BK000964.3` and `NR_003279.1` (5S) in place of the human accessions
- GtRNAdb mouse tRNAs: `mm10-tRNAs.fa` from the `Mmusc10` genome page
- Additional Dfam clade filter terms: `@Rodentia|@Muridae|@Murinae|@Mus_genus|@Mus_musculus`

**Outputs** (written to `OTHERS/Repeat_Mapping/` per sample):

| File | Description |
|------|-------------|
| `{sample}_repeat.bam` | Reads that mapped to the repeat index (sorted, indexed) |
| `{sample}_repeat_stats.txt` | Bowtie2 alignment summary |
| `{sample}_repeat_elements.tsv` | Per-element read counts and RPM (element / class / family / raw_reads / rpm) |
| `{sample}_repeat_families.tsv` | Family-level aggregated counts (class / family / raw_reads / rpm), sorted by read count |

RPM is calculated as `(reads mapped to element / total input reads) × 10⁶`.

### Recommended annotation directory layout

```
/path/to/star_index/
├── Genome                    ← STAR index files
├── SA, SAindex, genomeParameters.txt
├── genome.fa                 ← pass to --genome-fasta
├── chrom.sizes               ← optional, for bedgraph
└── Repeat/
    ├── repeat.1.bt2
    ├── repeat.2.bt2
    ├── repeat.3.bt2
    ├── repeat.4.bt2
    ├── repeat.rev.1.bt2
    └── repeat.rev.2.bt2
```

---

## Peak Matrix Metrics

`COMBINED_PEAK_MATRIX.txt` contains up to 54+ metrics per peak:

| Metric | Prefix | Scope | Description |
|--------|--------|-------|-------------|
| Biological Complexity | `BC_` | Group | Samples in group with raw count > 0 — measures reproducibility |
| Total Count (raw) | `TC_` | Sample | Raw read starts overlapping the peak |
| Total Count (aggregate) | `TC_` | Group | Sum of raw counts across all group samples |
| Normalized Count | `NormedTC_` | Sample | RPM: `raw TC × (1,000,000 / mapped reads)` |
| Normalized Count (aggregate) | `NormedTC_` | Group | Sum of RPM-normalized counts across group |
| Coverage Sum | `CovSum_` | Sample | Sum of per-base bedgraph signal across peak |
| Coverage Sum (group avg) | `CovSum_` | Group | Per-base signal from group-averaged bedgraph |
| Coverage Mean | `CovMean_` | Both | `CovSum / peak length` |
| Coverage Max | `CovMax_` | Both | Highest single-base signal in peak |

> **NormedTC** (additive) measures total group intensity. **Cov** columns (mean-based) measure signal density and shape of the average replicate.

---

## Changelog

### v3.5.0
- **Wizard overhaul**: per-tool interactive wizards (`-w` flag) for CLIPittyClip, PREPittyPrep, MAPittyMap, and PEAKittyPeak. New top-level launcher `CLIPittyClip_wizard.sh` provides a tool-triage entry point. Wizards collect inputs, ask for analysis tracks (HOMER / CTK CIMS-CITS / Clink), let you spot-edit defaults by category, show a diff-vs-defaults view + the equivalent CLI command at the end, and launch the chosen tool automatically.
- **MAPittyMap scope expansion**: alignment module now produces collapsed BED and bedgraphs by default, and optionally runs CTK CIMS/CITS (`--run-cims`, `--run-cits`) and Clink (`--run-clink`) crosslink-site analysis. Output tree mirrors the full pipeline (`1_BAM/`, `2_BED/`, `3_COVERAGE/`, `4_CTK_Analysis/`, `5_Clink_Analysis/`). Legacy BAM-only behavior preserved via `--bam-only`. PEAKittyPeak's role is now purely peak calling; point its `--ctk-dir` at MAPittyMap's `4_CTK_Analysis/` to annotate peaks with crosslink-site counts.
- **Repeat element pre-filter** (`--filter-repeat`): replaces the former ncRNA filter with a dedicated repeat filter targeting truly repetitive sequences that STAR cannot place uniquely
  - Index composition: NCBI 45S rDNA (U13369.1), GtRNAdb tRNAs, Dfam/RepeatMasker TE consensus sequences filtered to mammal-relevant clades
  - miRNA, Y RNA, lncRNA, snRNA, and snoRNA pass through to STAR for normal genomic annotation
  - Off by default (opt-in with `--filter-repeat`); index auto-detected from `{STAR_INDEX}/Repeat/repeat.*.bt2`
- **`run_repeat_quantify()`**: per-sample repeat quantification immediately after filtering
  - Per-element TSV (`_repeat_elements.tsv`): element / class / family / raw_reads / RPM
  - Per-family TSV (`_repeat_families.tsv`): family-level aggregated counts sorted by read depth
  - Handles three reference name formats: Dfam (`Name#Class/Family`), GtRNAdb tRNA headers, NCBI rDNA
- Output directory renamed from `ncRNA_Mapping/` to `Repeat_Mapping/`
- **`03_BEDGRAPH` renamed to `03_COVERAGE`**: directory name now reflects that it holds all coverage-format outputs, including bedgraphs and bigWigs
- **`--xl-bigwig` flag**: generates per-sample strand-specific crosslink-site bigWig files in `03_COVERAGE/BigWig/`
  - Extracts read 5′ ends from dedup BAMs by strand (truncation-event counts, not RPM read coverage)
  - `{sample}_xl_pos.bw` (+ strand) and `{sample}_xl_neg.bw` (− strand) per sample
  - Requires `bedGraphToBigWig` in PATH (`ucsc-bedgraphtobigwig` added to install scripts); skips gracefully if absent
  - Suitable as direct input for BindingSiteFinder and similar nucleotide-resolution tools
- **`all_crosslinks.bed` output from Clink**: every position with ≥1 truncation event written to `{sample}_all_crosslinks.bed` alongside the FDR-filtered `_truncations.bed`; no significance threshold applied — suitable for PEKA and BindingSiteFinder
- Fixed CTK output directory using hardcoded path `5_CTK_Analysis` instead of the `DIR_CTK` variable (broke numbering when other crosslink modules shifted folder numbers)

### v3.4.0
- **PAR-CLIP mode** (`--parclip`): end-to-end support for 4-thiouridine CLIP data
  - Read layout: `[UMI][READ][2nt spacer][6mer barcode][adapter]`; requires `-u` (UMI length)
  - Two-pass fastp: pass 1 extracts UMI and trims barcode+adapter (bundled `lib/parclip_adapters.fa`, 43 MDX-O barcodes); pass 2 removes the 2nt random spacer
  - `--parclip-adapters`: supply a custom adapter FASTA to override the bundled set
- **Strand-aware substitution calling** in `pileup.py`: all 12 substitution types now reported on the correct genomic strand (T→C for PAR-CLIP, G→T for eCLIP oxidative damage, etc.)
- **`--group-peaks`** flag in `PEAKittyPeak.sh` and `CLIPittyClip.sh`: pool per-sample BEDs before peak calling to increase depth for low-replicate experiments
- **Ensembl/UCSC chromosome naming auto-detection**: `filter_canonical_chromosomes()` detects `chr1` vs `1` naming from the BAM header and applies the correct filter set; portable `awk`-based implementation removes the `grep -P` dependency
- Suppressed `tag2collapse.pl` verbose EM output in non-debug runs
- `--clink-multi-map` reverted to single-pass positional assignment (CLAM-style EM removed — too slow for large datasets)
- Fixed `--ctk-group` flag not recognized by `PEAKittyPeak.sh`

### v3.3.1
- **`--clink-multi-map`**: single-pass positional rescue of multi-mapped reads (NH:i:>1) in the Clink collapse step
  - Unique reads (NH:i:1) are deduped first; multi-mappers are hard-assigned to the candidate locus with the greatest unique-read pileup support (±50 bp window), then deduped and merged
  - Off by default — all analyses presented use unique reads only
  - Compatible with `--group-xlsite`
- **`--group-xlsite`** flag propagated correctly through all batch/group code paths
- Removed `benchmark_pipeline.sh` and `benchmark_pipeline_clink.sh` (superseded by test suite)

### v3.3.0
- **Clink pipeline** (`--run-clink`): Python-native crosslink site caller
  - `umi_tools directional` deduplication — more conservative than `tag2collapse.pl` EM algorithm
  - `pileup.py`: single BAM scan → `.npz` shared between CITS and CIMS (no duplicate scanning)
  - `cits.py`: truncation site calling with binomial test + Benjamini-Hochberg FDR
  - `cims.py`: deletions + all 12 substitution types from one pileup
  - `--run-cims-cits --run-clink` runs both pipelines for direct comparison
  - Output folder numbering adjusts automatically (5_Clink or 6_Clink)
  - Hard dependency check at startup with clear install instructions
- Install scripts: `pysam` and `umi_tools` added; Python pinned to `<3.13`

### v3.2.0
- **`--genome-fasta` flag**: explicit reference FASTA for `samtools calmd` (STAR never stores FASTA in its index directory)
- **STAR mismatch tuning**: fractional `--outFilterMismatchNoverReadLmax 0.1` replaces hard `--outFilterMismatchNmax 2`; symmetric `--scoreInsOpen/Base -1` added
- **`scale_factors.tsv` reset** before batch aggregation to prevent peak matrix corruption on reruns
- macOS AppleDouble fix (`._*`) in BAM/BED `find` commands
- Warning emitted when Bowtie2 + crosslink analysis are combined
- Hash-based FASTQ dedup engine (`lib/fastq_collapse_hash.py`) — O(n), no disk spill
- `PREPittyPrep.sh`: standalone preprocessing + GEO deposit mode
- Lazy gzip: internal steps use plain `.fastq`; compression only at final output
- Conditional `0_DEMUX_FASTQ`: only retained with `-k`

### v3.1.0 / v3.0.2
- `--peak-caller` flag: `homer` (default) or `ctk`
- UCSC track headers added to all bedgraph outputs
- Plain `.fastq` / `.fq` support in `-i` and `-d` modes
- `--run-cims-cits` replaces `--run-ctk` (deprecated alias preserved)
- Interactive wizard supports peak caller selection
- Organized numbered output folders in single-file mode
- Bug fixes: batch mode dedup inheritance, peak aggregation, log file paths

---

## License

GPL-3.0 — See [LICENSE](LICENSE) for details.
