# CLIPittyClip.sh options

Source: `CLIPittyClip.sh --help` and the knob registry in `lib/wizard.sh` (v3.5.3).
If the installed version differs, run `CLIPittyClip.sh --help` and use that output.

"Ask when" tells you when the option applies. Do not ask about an option outside its
condition.

## Required

| Option | Meaning | Ask when |
|---|---|---|
| `-i <fastq>` | One FASTQ file (gz or plain) | Input is one file or a pooled file |
| `-d <dir>` | Folder of FASTQ files (batch mode) | Input is a folder |
| `-x <dir>` | STAR or Bowtie2 genome index folder | Always |

Use `-i` or `-d`, never both.

## 1. Preprocessing

| Option | Default | Meaning | Ask when |
|---|---|---|---|
| `-t <int>` | 1 | Threads | Always |
| `-u <int>` | none | UMI length | Not eCLIP |
| `-a <seq or name>` | L32 `GTGTCAGTCACTTCCAGCGG` | 3' adapter sequence. Newer versions also accept `L32`, `L3`, `L3MOD`, `L7`/`L19`, `SRA3`. | Not eCLIP or `--parclip` |
| `-b <file>` | none | Barcode file, `Name<TAB>Sequence`. Turns on demultiplexing. | Pooled input |
| `--bc-len <int>` | from `-b` | Barcode length | Folder input with barcodes still in reads |
| `--spacer-len <int>` | 0 | Bases trimmed after the barcode and UMI | Layout has a spacer |
| `--bc-first` | off | Layout `[BC][UMI]` instead of `[UMI][BC]` | Layout has the barcode first |
| `--trim3 <int>` | 0 | Bases trimmed from the 3' end after adapter removal | Layout has bases between read and adapter (iCLIP3) |
| `--demux-mismatches <int>` | 1 | Barcode mismatches allowed | Pooled input |
| `--min-qual <int>` | 30 | fastp average quality threshold | Settings walk-through |
| `--no-dedup` | off | Turn off PCR deduplication | Data has no UMI |
| `--low-memory` | off | Sort-based dedup with flat RAM | More than about 30 M reads per sample, or low RAM |
| `-m <star\|bowtie2>` | star | Aligner | Index is Bowtie2 |
| `--align-mismatches <int>` | 2 | Max alignment mismatches (STAR only) | Settings walk-through |
| `--eclip <pe\|se>` | off | eCLIP mode | Protocol is eCLIP |
| `--parclip` | off | fPAR-CLIP mode, needs `-u` | Protocol is fPAR-CLIP |
| `--parclip-adapters <fa>` | `lib/parclip_adapters.fa` | Custom PAR-CLIP adapter FASTA | fPAR-CLIP with other linkers |

## 2. Filters and reports

| Option | Default | Meaning | Ask when |
|---|---|---|---|
| `--filter-repeat` | off | Remove reads that map to rRNA, tRNA, and repeat elements before genome alignment. Needs a repeat index in `<index>/` or `<index>/Repeat/`. | User asks about rRNA or repeats |
| `--no-chr-filter` | off | Keep all contigs. Required for scaffold-level genomes. | Genome has no canonical chromosome names |
| `--report-abundant-rna [list]` | off | Report the fraction of reads on abundant ncRNA classes. Removes nothing. Needs `--gtf`. | User asks about abundant RNA |
| `--report-abundant-rna-custom name:path[,...]` | off | Report on user BED or GTF regions | User has a custom region list |
| `--gtf <path>` | none | Ensembl GTF for the report | With `--report-abundant-rna` |
| `--rmsk <path>` | none | RepeatMasker table for the tRNA category | tRNA in the report and no `--trna-bed` |
| `--trna-bed <path>` | none | Genomic tRNA BED (GtRNAdb) | tRNA in the report |

## 3. Peak calling

| Option | Default | Meaning | Ask when |
|---|---|---|---|
| `--peak-caller <homer\|ctk>` | homer | Peak caller | Always |
| `--peak-caller-args "<args>"` | none | Extra arguments for the peak caller | Settings walk-through |

HOMER defaults: `-style factor -L 2 -localSize 1000 -minDist 50 -size 20 -fragLength 25`.
CTK `tag2peak.pl` defaults: `-big -ss --valley-seeking -gap 50`.
To change `-minDist`, `-size`, `-fragLength`, or `-gap`, put the new value in
`--peak-caller-args`. These values have no separate flags.

## 4. CTK CIMS and CITS

Ask only when the user picked a CTK track.

| Option | Default | Meaning |
|---|---|---|
| `--run-cims-cits` | off | CIMS and CITS |
| `--run-cims` | off | CIMS only |
| `--run-cits` | off | CITS only |
| `--cims-iter <int>` | 5 | CIMS permutation iterations |
| `--cims-fdr <float>` | 0.05 | CIMS FDR threshold |
| `--cits-pval <float>` | 0.05 | CITS p-value threshold |
| `--cits-gap <int>` | 25 | CITS clustering gap, -1 turns off clustering |

## 5. Clink

Ask only when the user picked Clink.

| Option | Default | Meaning |
|---|---|---|
| `--run-clink` | off | Clink: umi_tools dedup, pileup, CITS and CIMS with BH FDR |
| `--genome-fasta <fa>` | none | Reference FASTA for `samtools calmd`. Also used by CTK. |
| `--clink-umi-len <int>` | auto | UMI length for umi_tools |
| `--clink-fdr <float>` | 0.05 | BH FDR threshold |
| `--clink-min-cov <int>` | 5 | Minimum coverage to test a position |
| `--clink-multi-map` | off | Assign multi-mapped reads with EM before dedup |
| `--clink-threads <int>` | min(-t, 8) | Pileup workers. RAM grows with workers. |
| `--xl-bigwig` | off | Per-sample crosslink-site bigWigs in `03_COVERAGE/BigWig/` |

## 6. Motif and flanked BED

| Option | Default | Meaning |
|---|---|---|
| `-f <int>` | 10 | Flank size in nt for flanked BED files |
| `--no-motif` | off | Do not make flanked BED files |

## 7. Grouping

Ask only for pooled or folder input.

| Option | Default | Meaning |
|---|---|---|
| `-g <file>` | none | Groups file, `SampleName<TAB>GroupName`. Makes group-averaged tracks and peak-matrix metrics. |
| `--group-xlsite` | off | Pool samples by group for CTK and Clink crosslink-site calls. Needs `-g`. Group CTK can need more than 64 GB RAM. |

## 8. Output and runtime

| Option | Default | Meaning |
|---|---|---|
| `-o <name or path>` | `<input>_output` next to the input | A name puts the folder next to the input. A path is used as given. |
| `-k` | off | Keep intermediate files |
| `-s <int>` | all reads | Test mode: process only the first N reads |
| `--notification` | off | System notification when the run ends |
| `-v` | off | Verbose log |
