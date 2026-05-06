#!/usr/bin/env python3
"""
Clink cits.py — Crosslink-Induced Truncation Sites (strand-aware)

Calls statistically significant RT truncation sites from a pileup.npz or BAM.
Processes forward (+) and reverse (-) strand independently.

Algorithm:
  1. Load stranded pileup (pileup.py) — separate arrays per strand
  2. Estimate genome-wide background truncation rate λ from both strands combined
     (excluding top 1% signal positions to avoid inflating background)
  3. Binomial test at every position per strand: P(X ≥ observed | Binomial(cov, λ))
  4. Benjamini-Hochberg FDR correction
  5. Write BED6+ output with correct strand column (+/-)

Truncation site convention (matching CTK):
  Forward read: crosslink = reference_start - 1
  Reverse read: crosslink = reference_end

Output:
  {prefix}_truncations.bed
    #chrom  start  end  name  score  strand  signal  coverage  fraction  pvalue  qvalue
    score = min(-log10(qvalue) * 100, 1000)   BED-compatible 0-1000
    strand = '+' (fwd reads) or '-' (rev reads)

Usage:
    python cits.py --pileup sample_pileup.npz --prefix results/SAMPLE
    python cits.py --bam sample_dedup.bam --prefix results/SAMPLE --chrom chr1
"""

import sys
import argparse
from pathlib import Path

from pileup import scan_bam, to_arrays, load_pileup
from stats import estimate_background, test_signal, write_bed, HEADER

import numpy as np


def run_cits(pileup_path:  str   = None,
             bam_path:     str   = None,
             prefix:       str   = None,
             chrom:        str   = None,
             min_mapq:     int   = 255,
             max_nh:       int   = 1,
             min_coverage: int   = 5,
             min_fraction: float = 0.05,
             min_signal:   int   = 1,
             fdr:          float = 0.05,
             verbose:      bool  = True):
    """
    End-to-end CITS: pileup (.npz) or BAM → significant truncation sites BED.

    Strand-aware: fwd (+) and rev (-) strand signals are tested independently
    using a shared genome-wide background rate (computed across both strands).

    Args:
        pileup_path  : .npz file from clink pileup (preferred)
        bam_path     : sorted, deduplicated BAM (used if pileup_path is None)
        prefix       : output file prefix (default: derived from input name)
        chrom        : restrict to one chromosome (default: all in pileup)
        min_mapq     : minimum mapping quality — only used with --bam
        max_nh       : max NH tag — only used with --bam
        min_coverage : minimum coverage to test a position
        min_fraction : minimum raw truncation fraction pre-filter
        min_signal   : minimum raw truncation read count at a position
                       (default 1 = no extra filter)
        fdr          : Benjamini-Hochberg FDR threshold
        verbose      : print progress to stderr
    """
    if pileup_path is None and bam_path is None:
        print("ERROR: supply either --pileup or --bam", file=sys.stderr)
        sys.exit(1)

    src_name = Path(pileup_path or bam_path).name
    if prefix is None:
        stem = Path(pileup_path or bam_path).stem
        prefix = stem.replace('_pileup', '')

    print(f"\nClink cits  |  {src_name}", file=sys.stderr)
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  min_signal={min_signal}  fdr={fdr}",
          file=sys.stderr)

    # --- Load pileup data ---
    if pileup_path:
        chrom_data_full = load_pileup(pileup_path)
    else:
        print(f"  mapq>={min_mapq}  nh<={max_nh}", file=sys.stderr)
        pileups = scan_bam(bam_path, chrom=chrom, min_mapq=min_mapq,
                           max_nh=max_nh, verbose=verbose)
        chrom_data_full = {}
        for c, pileup in pileups.items():
            result = to_arrays(pileup)
            if result is None:
                continue
            chrom_data_full[c] = result

    # --- Filter to requested chrom if specified ---
    if chrom:
        chrom_data_full = {c: v for c, v in chrom_data_full.items() if c == chrom}

    # --- Build strand-separated data, accumulate genome-wide for background ---
    # chrom_fwd/chrom_rev: {chrom: (positions, coverage, truncations)}
    chrom_fwd = {}
    chrom_rev = {}
    g_cov   = []   # combined both strands for background estimation
    g_trunc = []

    for c, strands in chrom_data_full.items():
        fwd = strands.get('fwd')
        rev = strands.get('rev')
        if fwd is not None:
            positions, coverage, truncations = fwd[0], fwd[1], fwd[2]
            chrom_fwd[c] = (positions, coverage, truncations)
            g_cov.append(coverage)
            g_trunc.append(truncations)
        if rev is not None:
            positions, coverage, truncations = rev[0], rev[1], rev[2]
            chrom_rev[c] = (positions, coverage, truncations)
            g_cov.append(coverage)
            g_trunc.append(truncations)

    if not chrom_fwd and not chrom_rev:
        print("  No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background rate (both strands pooled) ---
    all_cov   = np.concatenate(g_cov)
    all_trunc = np.concatenate(g_trunc)
    lambda_trunc = estimate_background(all_trunc, all_cov, min_coverage)

    print(f"\n  Background truncation rate (genome-wide, both strands):",
          file=sys.stderr)
    print(f"    λ = {lambda_trunc:.6f}  "
          f"(1/{int(1/lambda_trunc) if lambda_trunc > 0 else '∞'} reads)",
          file=sys.stderr)

    if lambda_trunc <= 0:
        print("  ERROR: background rate is zero — no testable positions.",
              file=sys.stderr)
        return

    # --- Test each strand per chromosome, write to single output BED ---
    out_path = f"{prefix}_truncations.bed"
    total_fwd = total_rev = 0

    all_chroms = sorted(set(list(chrom_fwd.keys()) + list(chrom_rev.keys())))

    with open(out_path, 'w') as fh:
        fh.write(HEADER)

        for c in all_chroms:
            if c in chrom_fwd:
                positions, coverage, truncations = chrom_fwd[c]
                results = test_signal(
                    positions, truncations, coverage,
                    lambda_trunc, min_coverage, min_fraction, fdr,
                    min_signal=min_signal)
                write_bed(results, c, 'trunc', fh, strand='+')
                total_fwd += len(results)

            if c in chrom_rev:
                positions, coverage, truncations = chrom_rev[c]
                results = test_signal(
                    positions, truncations, coverage,
                    lambda_trunc, min_coverage, min_fraction, fdr,
                    min_signal=min_signal)
                write_bed(results, c, 'trunc', fh, strand='-')
                total_rev += len(results)

    total = total_fwd + total_rev
    print(f"\n  Significant truncation sites (FDR < {fdr}): {total:,}  "
          f"(+: {total_fwd:,}  -: {total_rev:,})", file=sys.stderr)
    print(f"  Output: {out_path}", file=sys.stderr)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink cits: call RT truncation sites (strand-aware) from a pileup or BAM')

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument('--pileup',
        help='Pre-computed pileup .npz (output of clink pileup) — preferred')
    src.add_argument('--bam',
        help='Sorted, deduplicated BAM — scanned on the fly (slower)')

    parser.add_argument('--prefix', default=None,
        help='Output file prefix (default: derived from input filename)')
    parser.add_argument('--chrom', default=None,
        help='Restrict to one chromosome (default: all)')
    parser.add_argument('--mapq', type=int, default=255,
        help='Minimum MAPQ — BAM mode only (default: 255)')
    parser.add_argument('--nh', type=int, default=1,
        help='Max NH tag — BAM mode only (default: 1)')
    parser.add_argument('--min-cov', type=int, default=5,
        help='Minimum coverage to test a position (default: 5)')
    parser.add_argument('--min-frac', type=float, default=0.05,
        help='Minimum raw truncation fraction pre-filter (default: 0.05)')
    parser.add_argument('--min-signal', type=int, default=1,
        help='Minimum raw truncation read count at a position (default: 1 = no extra filter)')
    parser.add_argument('--fdr', type=float, default=0.05,
        help='BH FDR threshold (default: 0.05)')

    args = parser.parse_args()

    run_cits(
        pileup_path  = args.pileup,
        bam_path     = args.bam,
        prefix       = args.prefix,
        chrom        = args.chrom,
        min_mapq     = args.mapq,
        max_nh       = args.nh,
        min_coverage = args.min_cov,
        min_fraction = args.min_frac,
        min_signal   = args.min_signal,
        fdr          = args.fdr,
    )
