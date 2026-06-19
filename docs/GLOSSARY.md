# CLIPittyClip Domain Glossary

> Authoritative vocabulary for CLIPittyClip v3.5.
> Written following Domain-Driven Design (DDD) principles.
> All terms are defined as they are used in this codebase. 
> When a term diverges from how an external tool uses it, that divergence is noted.

---

## Table of Contents

1. [The Core Domain: What CLIPittyClip Does](#1-the-core-domain)
2. [Biological Concepts](#2-biological-concepts)
3. [Read-Level Concepts](#3-read-level-concepts)
4. [Pipeline Stages](#4-pipeline-stages)
5. [Signal Types](#5-signal-types)
6. [Crosslink Site Callers](#6-crosslink-site-callers)
7. [Statistical Concepts](#7-statistical-concepts)
8. [File Types and Data Structures](#8-file-types-and-data-structures)
9. [Experimental Protocols](#9-experimental-protocols)
10. [External Tools Used](#10-external-tools-used)
11. [Configuration and Parameters](#11-configuration-and-parameters)
12. [Bounded Contexts](#12-bounded-contexts)

---

## 1. The Core Domain

CLIPittyClip takes a raw CLIP-seq FASTQ and produces a **peak matrix**: a quantified table of RNA-binding protein (RBP) binding regions across all samples. Along the way the pipeline also produces crosslink site files (CITS, CIMS), strand-separated bedgraphs, and per-sample peak BEDs, all ready for downstream analysis and visualization.

```
Raw FASTQ
   │
   ▼ Preprocessing (dedup, trim, align)
Aligned BAM
   │
   ├──▶ Crosslink sites  →  CITS BED, CIMS BED         (single-nucleotide resolution)
   ├──▶ Bedgraphs        →  RPM-normalized ± tracks     (genome browser visualization)
   └──▶ Peaks            →  Peak BEDs + Peak Matrix     (quantification across samples)
```

The three bounded contexts (major subsystems):

| Context | Input | Primary Output |
|---|---|---|
| **Preprocessing** | Raw FASTQ | Aligned, deduplicated BAM |
| **Crosslink site calling** | BAM | CITS + CIMS BED files per sample |
| **Peak calling** | BAM / BED | Bedgraphs, peak BEDs, peak matrix |

---

## 2. Biological Concepts

### CLIP-seq
**Cross-Linking and ImmunoPrecipitation sequencing.** The experimental technique that generates CLIPittyClip's input data. An RNA-binding protein (RBP) is UV-crosslinked to RNA in living cells, then immunoprecipitated, and the bound RNA fragments are sequenced. The short UV-crosslinked RNA fragments become the sequencing reads.

Variants supported: iCLIP, irCLIP, eCLIP (PE and SE), CoCLIP, BrdU-CLIP.

### RNA-Binding Protein (RBP)
The protein whose binding sites are being mapped. CLIPittyClip is agnostic to which RBP is studied as it operates on the reads that result from the CLIP experiment.

### Crosslink Site
A single genomic position where the RBP was UV-crosslinked to RNA. Crosslink sites are evidence of direct protein-RNA contact.

Distinguished from **peaks**: a peak is a broad enrichment region; a crosslink site is a single-nucleotide position identified by a specific biochemical signature (truncation or mutation) in the sequencing read.

### RT Truncation (Truncation Site / CITS Signal)
When reverse transcriptase (RT) copies the crosslinked RNA into cDNA, it stalls at the crosslink nucleotide and falls off — producing a read whose 5′ end is one nucleotide downstream of the crosslink. The **crosslink position is therefore read_start − 1** for forward-strand reads, or **read_end** for reverse-strand reads.

The accumulation of 5′ read ends at a given position is the **truncation signal**. Sites with statistically enriched truncations are called CITS (Crosslink-Induced Truncation Sites).

See also: `clean_truncations`.

### RT Deletion (Deletion Site / CIMS Signal)
RT sometimes reads through the crosslink rather than truncating, introducing a characteristic 1-nucleotide deletion in the cDNA at the crosslink position. This deletion appears as a CIGAR `D` operation in the aligned read.

Sites with statistically enriched deletions are called CIMS (Crosslink-Induced Mutation Sites). CIMS also covers substitutions (e.g. T→C in PAR-CLIP).

### UV Crosslink
The covalent bond between the RBP and RNA formed by UV irradiation (254 nm for standard CLIP; 365 nm for PAR-CLIP with 4-thiouridine incorporation). The crosslink is the direct cause of both RT truncation and RT deletion signals.

---

## 3. Read-Level Concepts

### Read
A short DNA sequence produced by the sequencer representing one RNA fragment bound by the RBP. Reads enter as FASTQ, are cleaned through preprocessing, aligned to the genome, and emerge as rows in a BAM file.

### UMI (Unique Molecular Identifier)
A short random nucleotide sequence (5–10 nt) embedded in the read before PCR amplification. Allows deduplication after amplification: reads with the same UMI mapping to the same position are PCR duplicates of one molecule.

Read layouts supported:
- Standard: `[UMI][BC][spacer][READ]`
- BC-first (`--bc-first`): `[BC][UMI][spacer][READ]`
- eCLIP PE: UMI pre-extracted to read header by upstream `eclipdemux`

### Barcode (BC)
A sample-specific nucleotide sequence (4–6 nt) in the read used during demultiplexing to route each read to its sample. The barcode identifies the sample; the UMI identifies the molecule.

### PCR Duplicate
An amplification copy of an already-counted molecule. Removed in two passes:
1. **FASTQ-level** — exact sequence collapse before alignment (`_fastq_collapse_core`)
2. **BAM-level** — UMI + genomic position after alignment (`umi_tools dedup` for Clink; `tag2collapse.pl` for CTK)

### Tag
CTK's term for a uniquely mapped, PCR-deduplicated read as a BED interval. The collapsed BED (`*_collapsed.bed`) is a tag file. CLIPittyClip uses "tag" only in the CTK context.

### Collapsed Read / Collapsed Tag
A read that has survived PCR deduplication. Its original molecule count is stored in the read name as `READ#COUNT#UMI` (CTK format). The collapsed BED is the input to all CTK site-calling tools.

### Cluster
A set of reads (tags) whose genomic positions are close enough to be considered a single pile of evidence for one binding event. In CTK, clusters are built by `tag2cluster.pl` using a maximum gap (`-maxgap`). In Clink, clusters are implicit: positions within `cluster_gap` bp (default 25) are grouped for local background estimation.

Clusters are the unit of local context for CITS calling — truncation enrichment is tested relative to the cluster's own read density, not the genome-wide rate.

### Peak
A broad genomic region of significantly enriched read coverage, representing a candidate RBP binding region. Peaks are called by HOMER or CTK's `tag2peak.pl` and are wider than crosslink sites (typically 20–200 bp vs 1 bp).

Peaks are the unit of the **peak matrix**: each row is a peak, each column is a sample, values are read counts. The peak matrix is the primary output for differential binding and downstream analysis.

Crosslink sites (CITS/CIMS) sit inside peaks and provide single-nucleotide precision about where within a peak the RBP contacts the RNA.

### Clean Truncation
A truncation from a read with **no CIGAR D operation**. Used as the signal for Clink CITS. Deletion-carrying reads contribute to CIMS; excluding them from CITS avoids double-counting and mirrors CTK's `removeRow.pl` step.

Stored as the 6th element of the pileup strand tuple: `(positions, coverage, truncations, deletions, clean_truncations, subs)`.

---

## 4. Pipeline Stages

### Preprocessing
Converts raw FASTQ into an aligned, deduplicated BAM. Steps in order:

1. **FASTQ dedup** — exact sequence collapse before alignment
2. **Demultiplexing** — split pooled FASTQ by barcode (cutadapt)
3. **Adapter trimming** — remove 3′ adapter, extract UMI to header (fastp)
4. **ncRNA pre-filter** *(optional)* — discard rRNA/tRNA reads (Bowtie2)
5. **Genome alignment** — STAR or Bowtie2
6. **Chromosome filter** — keep chr1–22, X, Y, M only
7. **BAM-level dedup** — Clink: `umi_tools dedup`; CTK: `tag2collapse.pl`

### Demultiplexing
Splitting a pooled FASTQ into per-sample FASTQs by barcode sequence (cutadapt `--action=none` — reads are routed, not modified). Barcode collision detection (`check_barcodes.sh`) runs first. Two modes: **standard demux** (reads modified downstream) and **GEO demux** (`run_geo_demux`, writes reads as-is for public deposition).

### Adapter Trimming
Removes the 3′ adapter, applies quality filtering (Q30 average, 16 nt minimum), and moves the UMI from the read sequence to the header (delimited `#`). Performed by fastp.

### Alignment / Mapping
Maps trimmed reads to the reference genome. **STAR** (default) is splice-aware; **Bowtie2** uses BWA-aln–equivalent parameters tuned for CIMS. Key CLIP-specific settings: `EndToEnd` alignment mode; reduced gap penalties to preserve deletion-carrying reads.

### calmd (MD Tag Standardization)
`samtools calmd` recalculates MD tags from the reference FASTA, ensuring `parseAlignment.pl` correctly classifies mutations at deletion boundaries. Triggered by `--genome-fasta`. Strongly recommended for CIMS.

### Chromosome Filtering
Keeps only canonical chromosomes (chr1–22, X, Y, M). Removes unplaced contigs that inflate background rates and slow downstream tools.

---

## 5. Signal Types

### Coverage
The number of reads spanning a given genomic position. Used as the denominator (total trials) in binomial significance testing. A position must have `coverage >= min_cov` to be tested.

### Truncation Count
The number of reads whose 5′ end (RT stop) maps to a given position. The raw truncation signal before statistical testing.

- Forward strand: truncation at `reference_start - 1`
- Reverse strand: truncation at `reference_end`

### Clean Truncation Count
Truncation count restricted to reads with no CIGAR `D` operation. This is the signal used by Clink CITS. See *Clean Truncation*.

### Deletion Count
The number of reads carrying a 1-bp CIGAR `D` operation at a given position. The CIMS signal. Does not include intron-spanning CIGAR `N` operations.

### Substitution Count
The number of reads carrying a mismatch of a specific type (e.g. T→C) at a given position. Reported per substitution type (12 possible: AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG). Used for PAR-CLIP (T→C signature) and can be informative for iCLIP.

### Truncation Fraction
`truncations / coverage` at a position. Must be ≥ `min_frac` (default 0.05) to be tested. Pre-filters noisy positions before statistical testing.

### Signal Fraction
Generic term for the per-type fraction at a position: `signal / coverage`. Applies to truncations, deletions, or substitutions.

---

## 6. Crosslink Site Callers

### CTK (CLIP Tool Kit)
The Chaolin Zhang lab's Perl-based toolkit for CLIP-seq crosslink site analysis. CLIPittyClip wraps CTK tools for the `--run-cims-cits` track. Key tools used:

| CTK Tool | Role |
|---|---|
| `parseAlignment.pl` | Converts SAM → BED with mutation annotations |
| `tag2collapse.pl` | PCR dedup using UMI + position (CTK BAM-level dedup) |
| `CIMS.pl` | Calls deletion/mutation sites; permutation-based FDR |
| `CITS.pl` (manual steps) | Calls truncation sites; wrapped manually due to coordinate bug |
| `tag2cluster.pl` | Groups reads into clusters by overlap |
| `tag2peak.pl` | Calls peaks within clusters |
| `bedExt.pl` | Extends/shifts BED intervals |
| `removeRow.pl` | Filters BED rows by a condition (used to remove deletion-carrying reads) |

### Clink
CLIPittyClip's own Python crosslink site caller (`lib/clink/`). Runs parallel to CTK on the same BAM. Key design goals: strand-aware, BH FDR, local cluster background, compatible with `umi_tools` dedup.

Modules:

| Module | Role |
|---|---|
| `pileup.py` | BAM → stranded numpy arrays (coverage, truncations, deletions, subs) |
| `cits.py` | CITS: local cluster background + binomial test + BH FDR |
| `cims.py` | CIMS: global background + binomial test + BH FDR |
| `stats.py` | Shared statistical functions: `estimate_background`, `test_signal`, `bh_fdr`, `write_bed` |

### CITS (Crosslink-Induced Truncation Sites)
Statistical method that identifies positions with significantly elevated RT truncation rates. The enrichment is tested relative to a local background (Clink) or local cluster (CTK), not genome-wide, because truncation rate varies by local read density.

CTK algorithm:
1. Remove deletion-carrying reads (`removeRow.pl -q 3 -f 3`)
2. Shift 5′ read ends by -1 (`bedExt.pl -n up -l -1 -r -1`)
3. Cluster reads by overlap (`tag2cluster.pl -maxgap -1`)
4. Filter low-depth clusters (`awk '$5>2'`)
5. **1bp cluster boundary fix** (added in v3.3): extend each cluster 1bp on its 5′ end
6. Test truncation enrichment within clusters (`tag2peak.pl -gap 25 -p 0.001`)

Clink algorithm:
1. Extract `clean_truncations` from pileup
2. Compute per-position local background rate within read clusters (≤25bp gap)
3. Binomial test vs local rate per strand
4. BH FDR correction

### CIMS (Crosslink-Induced Mutation Sites)
Statistical method that identifies positions with significantly elevated deletion (or substitution) rates. Tested relative to genome-wide background rate.

CTK CIMS: permutation-based FDR (shuffle reads N times, measure false discovery rate empirically).
Clink CIMS: Benjamini-Hochberg FDR on binomial p-values.

### Local Cluster Background
In Clink CITS, instead of a single genome-wide background truncation rate (λ_global), each position is tested against a **per-cluster rate** computed from all reads within ≤25bp. This mirrors CTK's philosophy: truncation enrichment is assessed relative to local read density, not the genome average.

Positions in clusters with zero truncations get `local_rate = 0` and are skipped entirely (no test performed).

Computed by `build_local_rates()` in `cits.py`.

### Global Background Rate (λ)
A genome-wide estimate of the baseline signal fraction. Computed by `estimate_background()`:
- Only positions with `coverage >= min_cov` are included
- Top 1% of positions by signal fraction are excluded (to prevent genuine crosslink sites from inflating the background)
- If the 99th percentile is zero (very sparse signal), falls back to raw rate over all covered positions

Used by Clink CIMS (global λ for deletions/subs). Also computed for CITS as a reference but not used directly for testing.

---

## 7. Statistical Concepts

### Binomial Test
The significance test used by Clink. At each position: `p = P(X ≥ signal | Binomial(coverage, λ))`. Tests whether the observed signal count could arise by chance given the background rate λ and total coverage.

### Benjamini-Hochberg (BH) FDR
Multiple testing correction used by Clink. Controls the expected false positive rate across all tested positions. Produces q-values; a site is significant if `q ≤ fdr_threshold`. Default: 0.05; Chaolin Zhang tutorial default: 0.001.

### Permutation FDR
The FDR method used by CTK's `CIMS.pl`. Reads are shuffled N times (tutorial: 10; publications: 200+) and the false discovery rate is estimated empirically. More conservative than BH FDR, particularly on sparse data.

### Pre-filter
Positions that fail any threshold are excluded before statistical testing — reducing the test burden and preserving FDR budget for testable positions:
- `min_cov` (default 5) — minimum coverage
- `min_frac` (default 0.05) — minimum signal/coverage fraction
- `min_signal` (default 1) — minimum raw signal count

### q-value
The BH-adjusted p-value. A site is called significant when `q ≤ fdr_threshold`. Controls the expected rate of false positives across the full set of calls, not for any individual site.

### Score (BED score column)
Clink output BEDs encode significance as `score = min(-log10(q_value) × 100, 1000)` — the standard 0–1000 BED range. Lower q-value → higher score → more prominent in genome browsers.

---

## 8. File Types and Data Structures

### FASTQ
Raw sequencing reads file. Each read is 4 lines: header (@...), sequence, +, quality scores. CLIPittyClip accepts `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`.

### BAM (Binary Alignment Map)
Compressed binary alignment file. Contains all aligned reads with their genomic positions, CIGAR strings, and optional tags (NH, NM, MD). The central data file between preprocessing and crosslink site calling. Must be sorted by coordinate and indexed (`.bai`).

### BED (Browser Extensible Data)
Tab-separated genomic interval file. CLIPittyClip uses BED6 (chrom, start, end, name, score, strand) for read tags and BED6+ for crosslink site output (adds signal, coverage, fraction, pvalue, qvalue columns).

Coordinate convention: **0-based half-open** [start, end). A single-nucleotide site at position P is `[P, P+1)`.

### Collapsed BED
The BED file of PCR-deduplicated reads produced by `tag2collapse.pl`. Read names are in CTK format: `READ#COUNT#UMI`. This is the primary input to CTK's CIMS and CITS analysis.

### Pileup NPZ
A numpy `.npz` file storing per-position signal arrays for all chromosomes and both strands. Produced by `pileup.py`. Avoids re-scanning the BAM when running multiple callers on the same data.

Array keys follow the pattern: `{chrom}__{strand}__{signal_type}`. Example: `chr1__fwd__coverage`, `chr1__rev__truncations`.

Strand tuple stored per chromosome: `(positions, coverage, truncations, deletions, clean_truncations, subs)` where `subs` is a dict of `(ref, alt) -> array`.

### BedGraph
Per-position coverage track in UCSC BedGraph format (chrom, start, end, value). CLIPittyClip produces RPM-normalized strand-separated bedgraphs (one `+`, one `−` per sample) for genome browser visualization.

### Peak Matrix
A tab-separated table (`COMBINED_PEAK_MATRIX.txt`) with one row per peak and one column per sample. Values are read counts per peak. Enhanced with biological complexity (BC), normalized counts, and group-averaged columns.

### CIGAR String
Alignment description encoded in the SAM/BAM format. Each operation has a code and length:
- `M` — alignment match (can be match or mismatch)
- `D` — deletion from reference (1-nt D = CIMS signal)
- `N` — skipped region (intron, not counted as deletion)
- `I` — insertion to reference
- `S` — soft clip (not aligned)

CLIPittyClip specifically distinguishes `D` (deletion) from `N` (intron) — only `D` is counted as a CIMS mutation signal.

### MD Tag
A SAM/BAM auxiliary tag describing mismatches relative to the reference. Used by `parseAlignment.pl` to classify substitution types. Can be inconsistent when produced by STAR at deletion boundaries — `samtools calmd` standardizes these from the reference FASTA.

---

## 9. Experimental Protocols

### iCLIP
Individual-nucleotide resolution CLIP. Uses UV crosslinking at 254 nm. Signal is primarily RT truncation (CITS) with some deletions (CIMS). Reads are single-end; UMI is at the 5′ end. The protocol that CTK was originally designed for.

### irCLIP / irCLIP2
Improved and irradiation CLIP. Variant of iCLIP with higher efficiency. irCLIP2 uses BC-first layout (`--bc-first`).

### eCLIP (Enhanced CLIP)
Enhanced CLIP developed by the ENCODE project (Van Nostrand lab). Paired-end sequencing; the protein-RNA complex is size-selected, producing well-defined RNA fragments.

Two eCLIP modes in CLIPittyClip:
- **PE mode** (`--eclip pe`): Input is post-`eclipdemux` Read 2. UMI has already been moved to the read header by `eclipdemux`. CLIPittyClip validates this format, re-embeds UMI in sequence for collapse, then strips and trims.
- **SE mode** (`--eclip se`): Input is raw Read 1 from seCLIP (Blue et al. 2022). UMI is the first 10 nt of the sequence.

### PAR-CLIP
Photoactivatable Ribonucleoside-Enhanced CLIP. Uses 4-thiouridine (4SU) incorporation and 365 nm UV. The primary signal is T→C substitution at the crosslink site (not truncation). CLIPittyClip's CIMS module reports all substitution types including T→C.

### CoCLIP
Colocalization CLIP protocol using a 7 nt UMI. The default UMI length in CLIPittyClip examples.

---

## 10. External Tools Used

| Tool | Role in CLIPittyClip |
|---|---|
| **fastp** | Adapter trimming, UMI extraction, quality filtering |
| **STAR** | Splice-aware genome alignment (default) |
| **Bowtie2** | End-to-end genome alignment (CIMS-tuned alternative) |
| **samtools** | BAM manipulation: view, sort, index, calmd |
| **cutadapt** | Barcode-based demultiplexing |
| **umi_tools** | BAM-level UMI deduplication (Clink track) |
| **CTK** | Chaolin Zhang lab's Perl toolkit: parseAlignment, tag2collapse, CIMS, CITS |
| **HOMER** | Peak calling (makeTagDirectory + findPeaks) |
| **pysam** | Python BAM reading in Clink pileup.py |
| **numpy / scipy** | Array math and binomial statistics in Clink |

---

## 11. Configuration and Parameters

### min_cov (--clink-min-cov, --min-cov)
Minimum read coverage at a position for it to be tested. Positions with fewer reads are skipped. Default: 5. Prevents calling sites from noise at low-coverage positions.

### min_frac (--min-frac)
Minimum signal fraction (signal / coverage) at a position. Default: 0.05 (5%). Pre-filter before binomial testing.

### fdr (--clink-fdr, --fdr, --cims-fdr)
Benjamini-Hochberg FDR threshold. Default: 0.05 (5%). Zhang Lab CTK tutorial default for CIMS post-hoc filter: 0.001.

### cluster_gap (--cluster-gap, --cits-gap)
Maximum bp gap between consecutive positions before they are considered separate clusters for local background estimation. Default: 25, matching CTK's `--gap 25`. A smaller gap creates tighter clusters (more conservative local background).

### cims_iter (--cims-iter)
Number of permutation iterations for CTK CIMS FDR. Default in CLIPittyClip: 5 (fast). Zhang Lab CTK tutorial: 10. Publications: 200+.

### Zhang Lab CTK Tutorial Defaults
The parameter set used by the CTK documentation at columbia.edu/~chaolin:
- CIMS: 10 iterations, post-hoc FDR < 0.001
- CITS: p = 0.001, gap = 25, no Bonferroni correction
- Clink (matched): FDR = 0.001, cluster_gap = 25, min_cov = 5, min_frac = 0.05

Encoded in `benchmark_pipeline_defaults.sh`.

---

## 12. Bounded Contexts

CLIPittyClip has three bounded contexts, each with it's own language model:

### Preprocessing Context
**Language**: reads, FASTQ, UMI, barcode, adapter, alignment, deduplication.
**Invariant**: a read must survive exact-sequence dedup → adapter trim → alignment → position-based dedup before entering the signal context.
**Key files**: `CLIPittyClip.sh`, `lib/dedup.sh`, `lib/modules.sh` (fastp, STAR, Bowtie2, calmd, umi_tools functions).

### Crosslink Context (Clink)
**Language**: position, coverage, truncation, deletion, substitution, strand, cluster, local rate, pileup.
**Invariant**: all arrays are strand-separated; a position only exists in the sparse array if at least one event was observed; clean_truncations is always a subset of truncations.
**Key files**: `lib/clink/pileup.py`, `lib/clink/cits.py`, `lib/clink/cims.py`, `lib/clink/stats.py`.

### Crosslink Context (CTK)
**Language**: tag, collapsed tag, mutation file, cluster, peak.
**Invariant**: the collapsed BED carries read counts in the read name; all CTK tools consume and produce BED6 files.
**Key files**: CTK Perl scripts wrapped in `lib/modules.sh` (run_parse_alignment, run_collapse_pcr, run_ctk_cims, run_ctk_cits).

### Peak Context
**Language**: peak, bedgraph, RPM, peak matrix, biological complexity (BC), group.
**Invariant**: peaks are genomic regions, not single-nucleotide; they are quantified across samples in the peak matrix.
**Key files**: `lib/modules.sh` (run_peak_calling_homer, run_peak_calling_ctk, add_matrix_columns), `PEAKittyPeak.sh`.
