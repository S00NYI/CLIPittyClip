# CLIPittyClip

A single-command CLIP-seq analysis pipeline that takes raw sequencing reads and produces a peak matrix — a quantified table of RNA-binding protein (RBP) binding regions across all samples — along with crosslink site files (CITS, CIMS) and bedgraphs for downstream analysis and visualization.

---

## Language

### Core entities

**FASTQ**:
The raw input — a file of sequencing reads, each with a nucleotide sequence and quality scores.
_Avoid_: reads file, input file

**Read**:
One sequencing read: a short DNA sequence representing one RNA fragment captured in the CLIP experiment.
_Avoid_: sequence, fragment (fragment is ambiguous — it means the RNA piece, not the sequenced representation)

**BAM**:
The aligned, deduplicated read file — the central handoff between preprocessing and all downstream analysis.
_Avoid_: alignment file, mapped reads file

**Crosslink Site**:
A single genomic position (1 bp) where the RBP was UV-crosslinked to RNA, identified by a statistically significant biochemical signature in the reads.
_Avoid_: binding site (too vague — peaks are also binding sites), mutation site (too narrow — truncations are not mutations)

**Peak**:
A broad genomic region (typically 20–200 bp) of significantly enriched read coverage representing a candidate RBP binding region. The unit of the peak matrix.
_Avoid_: binding region, enriched region

**Peak Matrix**:
The primary output — a table with one row per peak and one column per sample, values are read counts. Used for differential binding and downstream quantitative analysis.
_Avoid_: count matrix, coverage matrix

**Cluster**:
A set of reads whose positions are close enough to be treated as one pile of evidence for a single binding event. The unit of local context for CITS calling.
_Avoid_: read pile, read stack

**Tag**:
CTK-specific term for a PCR-deduplicated read represented as a BED interval. Only used in the CTK bounded context.
_Avoid_: using "tag" to mean a generic read outside the CTK context

**Collapsed Read**:
A read that has survived PCR deduplication. Its original molecule count is embedded in the read name as `READ#COUNT#UMI`.
_Avoid_: deduplicated read (ambiguous — could mean FASTQ-level or BAM-level dedup)

### Read anatomy

**UMI (Unique Molecular Identifier)**:
A short random nucleotide sequence (5–10 nt) in the read that identifies the original molecule, enabling PCR deduplication.
_Avoid_: barcode (barcode means the sample-demultiplexing sequence, not the UMI)

**Barcode (BC)**:
A sample-specific nucleotide sequence (4–6 nt) in the read used to route each read to its sample during demultiplexing.
_Avoid_: UMI, index (index is a sequencing instrument concept)

**Clean Truncation**:
A truncation event from a read with no CIGAR `D` operation. The signal used by Clink CITS to avoid double-counting reads that carry both a truncation and a deletion.
_Avoid_: pure truncation, non-deletion truncation

### Signal types

**Truncation**:
A 5′ read end caused by RT stalling at a crosslink. The CITS signal. Crosslink position = `read_start − 1` (forward strand) or `read_end` (reverse strand).
_Avoid_: RT stop (ambiguous — any RT stop, not necessarily at a crosslink)

**Deletion**:
A 1-bp CIGAR `D` operation in an aligned read caused by RT reading through a crosslink. The primary CIMS signal.
_Avoid_: mutation (too broad — substitutions are also mutations)

**Substitution**:
A mismatch of a specific nucleotide type (e.g. T→C in PAR-CLIP) caused by modified nucleotides at the crosslink. A secondary CIMS signal.
_Avoid_: SNP, variant

**Coverage**:
The number of reads spanning a genomic position. The denominator in all significance tests.
_Avoid_: depth, read depth (acceptable informally, but "coverage" is the term used in code)

**Signal Fraction**:
`signal / coverage` at a position. Used as a pre-filter (`min_frac`) and as an output column.
_Avoid_: rate (rate is reserved for background rate λ)

**Background Rate (λ)**:
The genome-wide (or cluster-level) baseline signal fraction, estimated from non-crosslink positions. Used as the null hypothesis in binomial testing.
_Avoid_: noise rate, baseline

### Statistical terms

**Local Cluster Background**:
A per-cluster background rate computed from all reads within ≤`cluster_gap` bp, used by Clink CITS instead of a genome-wide λ.
_Avoid_: local rate (acceptable shorthand internally), local lambda

**BH FDR**:
Benjamini-Hochberg false discovery rate correction. The multiple testing method used by Clink. Produces q-values.
_Avoid_: FDR alone when CTK's permutation FDR is also in scope — be specific

**Permutation FDR**:
CTK's empirical FDR method — reads are shuffled N times to estimate the null distribution. Used by `CIMS.pl`.
_Avoid_: FDR alone when Clink's BH FDR is also in scope

**q-value**:
A BH-adjusted p-value. A site is significant when `q ≤ fdr_threshold`.
_Avoid_: adjusted p-value (technically correct but not the term used in code or output columns)

---

## Relationships

- A **FASTQ** is preprocessed into a **BAM**
- A **BAM** is the shared input for both crosslink site calling and peak calling
- A **Read** in a BAM carries a **Truncation** (at its 5′ end), zero or one **Deletion** (CIGAR D), and zero or more **Substitutions** (mismatches)
- A **Clean Truncation** is a **Truncation** from a **Read** with no **Deletion**
- A **Cluster** is a group of **Reads** within `cluster_gap` bp of each other
- A **Crosslink Site** is a single genomic position inside a **Cluster** that passes significance testing
- A **Peak** spans one or more **Clusters** and contains zero or more **Crosslink Sites**
- A **Peak Matrix** has one row per **Peak** and one column per sample; values are **Coverage** counts
- A **Tag** is a **Collapsed Read** in the CTK bounded context; outside CTK, use **Collapsed Read**

---

## Example dialogue

> **Dev:** "Should I use 'reads' or 'tags' when naming this variable that holds post-dedup BAM alignments?"
> **Domain expert:** "Neither. Use 'collapsed reads' if you're in the CTK track (they carry the `#COUNT#UMI` format), or just 'reads' if you're in Clink. 'Tag' is only CTK's internal BED-format term."

> **Dev:** "I want to test whether there's signal at a position — should I use the global λ or local cluster background?"
> **Domain expert:** "Depends on the caller. CITS always uses local cluster background — truncation rate varies too much genome-wide. CIMS uses global λ — deletions are rare everywhere, so the genome-wide rate is a stable null."

> **Dev:** "The user said 'crosslink sites' but is looking at the peak matrix. Are those the same thing?"
> **Domain expert:** "No. Crosslink sites are 1-bp positions from CITS/CIMS. Peaks are broad regions from HOMER or tag2peak.pl. The peak matrix quantifies peaks, not crosslink sites."

---

## Flagged ambiguities

- **"FDR"** — used by both Clink (BH FDR) and CTK (permutation FDR). Always qualify: **BH FDR** or **permutation FDR**.
- **"dedup"** — happens twice: FASTQ-level (exact sequence collapse) and BAM-level (UMI + position). Always specify which stage.
- **"truncation"** vs **"clean truncation"** — `truncations` in the pileup includes all 5′ ends; `clean_truncations` excludes reads with CIGAR D. Clink CITS uses `clean_truncations`. Don't use "truncation" when you mean the clean subset.
- **"peak"** — overloaded in HOMER (internal data structure) vs CLIPittyClip (the genomic output region). In CLIPittyClip, a peak is always the output BED region, never HOMER's internal format.
- **"signal"** — generic term. In code, always prefer the specific type: truncation, deletion, substitution, or coverage.
