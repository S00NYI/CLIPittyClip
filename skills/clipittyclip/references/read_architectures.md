# CLIP read architectures and CLIPittyClip arguments

Source: "An Index of CLIP Read Architectures", J. Luna, v2 (`docs/CLIP_Read_Arch_index_v2.pdf`).
Check that file for the primer sequences and the cloning steps.

Key: `N` = UMI, `X` = barcode, sp = spacer, READ = the insert.

If the user gives a barcode file with `-b`, the pipeline computes the barcode length. Then
`--bc-len` is not necessary.

## Variants by family

### iCLIP family

| Variant | Read layout (5' to 3') | Arguments |
|---|---|---|
| iCLIP (Konig 2010) | `XX NNN READ L3` | `-u 3 --bc-len 2 --bc-first --spacer-len 0 -a TGAGATCGGAAGAGCGGTTCAG` |
| iCLIP2 (Buchbender 2020) | `NNNNN XXXXXX NNNN READ L3` | `-u 5 --bc-len 6 --spacer-len 4 -a TGAGATCGGAAGAGCGGTTCAG` |
| iCLIP3 (Despic 2026) | `NNNNNNNNN READ NNN L7` | `-u 9 --trim3 3 -a AGATCGGAAGAGCACACGTC` |
| irCLIP (Zarnegar 2016) | `NNNNN XXXXXX NNNWWW READ modL3` | `-u 5 --bc-len 6 --spacer-len 6 -a AGATCGGAAGAGCGGTTCAG` |
| irCLIP2 (Ducoli 2025) | `XXXXXX NNNNNNNN T(12) READ L7` | `-u 8 --bc-len 6 --bc-first --spacer-len 12 -a AGATCGGAAGAGCACACGTC` |

Notes:

- iCLIP2: the pipeline trims the second UMI (UMI2) as a spacer.
- iCLIP3: the sequencing core demultiplexes with i5/i7 indexes, so the input is one file
  per sample. `--trim3 3` removes the 3 nt UMI2 between the read and the linker.

### coCLIP and BrdU-CLIP

| Variant | Read layout (5' to 3') | Arguments |
|---|---|---|
| coCLIP (Yi 2024, Luna lab) | `NNNNNNN XXXXXX DDD READ L32` | `-u 7 --bc-len 6 --spacer-len 3` (L32 is the default adapter, omit `-a`) |
| BrdU-CLIP (Weyn-Vanhentenryck 2014) | `XXXXC NNNNNNNN C READ SRA3` | `-u 8 --bc-len 6 --bc-first --spacer-len 1 -a TCGTATGCCGTCTTCTGCTTG` |

Note: for BrdU-CLIP, the first anchor C counts as part of the barcode. The second anchor C
counts as the spacer.

### eCLIP

| Variant | Read layout (5' to 3') | Arguments |
|---|---|---|
| seCLIP (single-end, raw R1) | `NNNNNNNNNN READ InvRiL19` | `--eclip se` |
| peCLIP (paired-end, post-eclipdemux R2) | `READ InvRiL19`, UMI in the read header | `--eclip pe` |

Notes:

- `--eclip se` hardcodes a 10 nt UMI and the TruSeq R1 adapter. The pipeline ignores `-u`
  and `-a`.
- `--eclip pe` reads the UMI from the read header and uses `lib/eclip_adapters.fa`. The
  pipeline ignores `-u` and `-a`.
- The sequencing core demultiplexes eCLIP libraries. The input is a file or a folder,
  not a pooled file.

### PAR-CLIP and HITS-CLIP

| Variant | Read layout (5' to 3') | Arguments |
|---|---|---|
| fPAR-CLIP (Anastasakis 2017) | `NNNN READ NN XXXXXX adapter` | `--parclip -u 4` |
| Legacy PAR-CLIP or HITS-CLIP (no UMI) | `READ SRA3` | `--no-dedup -a TCGTATGCCGTCTTCTGCTTG` |

Notes:

- `--parclip` trims the 3' barcode and the adapter with `lib/parclip_adapters.fa`
  (43 barcoded 3' linkers). Use `--parclip-adapters <fa>` for other linkers.
- `--parclip` requires `-u`.
- Legacy data has no UMI, so PCR deduplication by UMI is not possible. Use `--no-dedup`.

## Adapter sequences

| Name | Sequence | Used by |
|---|---|---|
| L32 | `GTGTCAGTCACTTCCAGCGG` | coCLIP (pipeline default) |
| L3 | `TGAGATCGGAAGAGCGGTTCAG` | iCLIP, iCLIP2 |
| modified L3 (`L3MOD`) | `AGATCGGAAGAGCGGTTCAG` | irCLIP |
| L7 / InvRiL19 | `AGATCGGAAGAGCACACGTC` | iCLIP3, irCLIP2, eCLIP |
| SRA3 | `TCGTATGCCGTCTTCTGCTTG` | BrdU-CLIP, legacy PAR-CLIP |
| fPAR-CLIP 3' adapter | `TGGAATTCTCGGGTGCCAAGG` | fPAR-CLIP (in `parclip_adapters.fa`) |

## Custom layouts

If the protocol is not in this file, ask the user for the read layout from 5' to 3'.
Map each part to an argument:

| Part | Argument |
|---|---|
| UMI length | `-u` |
| Barcode length | `--bc-len` (or omit when `-b` is given) |
| Barcode before UMI | `--bc-first` |
| Fixed or random bases after the barcode and UMI | `--spacer-len` |
| Bases between the read and the adapter | `--trim3` |
| 3' adapter | `-a <sequence>` |
| No UMI | `--no-dedup` |

The pipeline supports two 5' orders: `[UMI][BC][spacer][READ]` (default) and
`[BC][UMI][spacer][READ]` (`--bc-first`). If the layout fits neither order, tell the user.
Then suggest preprocessing outside the pipeline.
