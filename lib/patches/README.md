# STAR on macOS: the libc++ fix

Stock STAR 2.7.11b on macOS **reads zero reads from any input and still reports success**
(Linux is unaffected). `install_macos.sh` detects this with a one-second alignment test and,
only if STAR is broken, builds a patched copy from the official source.

| File | Role |
|------|------|
| `star-2.7.11b-libcxx-streambuf.patch` | The fix (3 files, +78/-8), vendored unmodified from upstream PR #2691 |
| `../star_selftest.sh` | Aligns 2,000 reads to a 100 kb synthetic genome; exit 0 only if STAR really aligns |
| `../build_patched_star.sh` | Downloads official STAR, verifies SHA-256, applies the patch, compiles, installs |
| `../modules.sh` `run_mapping_star` | Stops the pipeline if STAR reports 0 input reads on a non-empty FASTQ |

## Symptom

The run finishes without error but everything is empty: STAR's `*.Log.final.out` says
`Number of input reads | 0` (its `Log.out` ends `Thread #1 end of input stream, nextChar=-1`), then
Clink reports `Total input alignments: 0`, every pileup shows 0 reads, and peak calling finds
nothing. Re-running into an existing output folder overwrites good results with empty ones.

## Cause

STAR attaches its pre-allocated read/BAM buffers to string streams with
`rdbuf()->pubsetbuf(buf, n)`. That call is implementation-defined: libstdc++ (GCC, Linux) adopts
the buffer; **libc++ (clang, all macOS) silently ignores it**, so the streams stay empty.

STAR 2.7.10b worked only because it is an older GCC-era binary that bioconda ships for Intel
(`osx-64`) alone, so it ran under Rosetta. **No native Apple Silicon build of 2.7.10b exists**, so
pinning it cannot help there. Every 2.7.11b build we tried fails the same way (conda `_6`/`_7`/`_8`,
Homebrew `rna-star`). Upstream: STAR issue
[#2632](https://github.com/alexdobin/STAR/issues/2632) (the cause is in the comments); fix proposed
in [PR #2691](https://github.com/alexdobin/STAR/pull/2691) (closed, not merged).

## Fix and provenance

The patch adds a small `std::streambuf` (`ExtBufStreambuf.h`) that adopts the buffer explicitly, so
behavior is the same on every standard library. It changes how STAR *buffers* data, not how it
aligns.

* Patch: PR #2691 by `BenjaminDEMAILLE`, byte-for-byte (SHA-256 `ed4415b0...8091c9`; the builder
  refuses a modified file). STAR is MIT-licensed. The PR says it was written with an AI assistant,
  so we tested it ourselves rather than trust it.
* Source: official `STAR/archive/refs/tags/2.7.11b.tar.gz`, checked against SHA-256
  `3f65305e...c408e7` (the value Homebrew pins) before use.

## Verification (Oct 2026, Apple Silicon, macOS 27)

| Check | Result |
|-------|--------|
| Stock conda and Homebrew STAR, synthetic test | 0/2000 reads (fail) |
| Patched STAR, same test | 2000/2000 |
| 1,399,276 real reads (JL0380 CytoArs_2) through the pipeline's exact fastp + STAR arguments, GRCh38 | all reads read (equals fastp's count); 71.0% unique |
| Same reads: 300 kB buffers (many chunk boundaries) / 1 thread | Alignments identical to default, all 1,995,827 records (name, flag, chr, pos, MAPQ, CIGAR, NH) |
| `--readFilesCommand gzip -dc` | Works (fails on stock conda STAR) |
| Full pipeline subsample (`-s 300000`, `--run-clink --group-xlsite`) | Passed |

**Not yet done:** comparison against an *unpatched* STAR 2.7.11b on Linux. The above shows the
build is self-consistent and aligns correctly, not that it matches a reference read for read.
Before publishing, align one dataset on Linux and diff (`samtools view | cut -f1-6 | sort | md5`).

## How it is used

```
install_macos.sh -> conda env (STAR 2.7.10b on Intel, 2.7.11b on Apple Silicon)
                 -> star_selftest.sh        pass: done, nothing built
                                            fail: build_patched_star.sh (~30 s on 10 cores), then re-test;
                                                  the install stops if STAR still fails
```

* The previous STAR is kept as `$CONDA_PREFIX/bin/STAR.orig`; build details are in
  `$CONDA_PREFIX/bin/STAR.patch-info.txt`.
* Re-running the installer is safe and repairs a broken STAR in an existing environment.
  `mamba install star` / `update star` restores the broken binary: re-run the installer or the builder.
* Needs Xcode Command Line Tools, `llvm-openmp` (installed by the installer) and GitHub access.

Manual use:

```bash
lib/star_selftest.sh                                       # check the STAR on PATH (or pass a path)
conda activate clipittyclip && mamba install llvm-openmp
lib/build_patched_star.sh --dest "$CONDA_PREFIX/bin"       # only if the test fails
```

## Removing this

Delete the patch, builder, the installer's Step 4b and this folder once a STAR release or a
bioconda/conda-forge rebuild fixes the bug (a fresh environment passes `star_selftest.sh` with no
build) and you no longer support older environments. Until then the build runs only when the
self-test fails, so it retires itself once conda ships a fixed STAR.

## Methods wording

> Reads were aligned with STAR 2.7.11b built from source with a patch fixing stream-buffer handling
> under libc++ (alexdobin/STAR pull request 2691; CLIPittyClip `lib/patches/`), because the
> unpatched macOS build returns no alignments.

On Linux state the version only. On Intel macOS the installer uses STAR 2.7.10b.
