#!/usr/bin/env python3
"""
Clink cits.py — Crosslink-Induced Truncation Sites

Calls statistically significant RT truncation sites from a sorted,
deduplicated BAM. Replaces CTK CITS.pl in the CLIPittyClip pipeline.

Algorithm:
  1. Scan BAM → per-position truncation counts and coverage (pileup.py)
  2. Estimate genome-wide background truncation rate λ (excluding top 1% sites)
  3. Binomial test at every position: P(X ≥ observed | Binomial(coverage, λ))
  4. Benjamini-Hochberg FDR correction
  5. Write BED6+ output

Truncation site convention (matching CTK):
  The crosslink nucleotide = last base synthesized by RT = one position
  upstream of the read's 5' end:
    Forward read: reference_start - 1
    Reverse read: reference_end

Output:
  {prefix}_truncations.bed
    #chrom  start  end  name  score  strand  signal  coverage  fraction  pvalue  qvalue
    score = min(-log10(qvalue) * 100, 1000)   BED-compatible 0-1000

Usage:
    python cits.py --bam sample_dedup.bam --prefix results/SAMPLE
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
             fdr:          float = 0.05,
             verbose:      bool  = True):
    """
    End-to-end CITS: pileup (.npz) or BAM → significant truncation sites BED.

    Prefer --pileup (pre-computed .npz from clink pileup) to avoid re-scanning
    the BAM. Falls back to --bam for standalone use.

    Args:
        pileup_path  : .npz file from clink pileup (preferred)
        bam_path     : sorted, deduplicated BAM (used if pileup_path is None)
        prefix       : output file prefix (default: derived from input name)
        chrom        : restrict to one chromosome (default: all in pileup)
        min_mapq     : minimum mapping quality — only used with --bam
        max_nh       : max NH tag — only used with --bam
        min_coverage : minimum coverage to test a position
        min_fraction : minimum raw truncation fraction pre-filter
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
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  fdr={fdr}",
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
            positions, coverage, truncations, deletions, subs = result
            chrom_data_full[c] = (positions, coverage, truncations, deletions, subs)

    # --- Filter to requested chrom if specified ---
    if chrom:
        chrom_data_full = {c: v for c, v in chrom_data_full.items() if c == chrom}

    # --- Extract truncation signal, accumulate genome-wide ---
    chrom_data = {}
    g_cov, g_trunc = [], []

    for c, arrays in chrom_data_full.items():
        positions, coverage, truncations = arrays[0], arrays[1], arrays[2]
        chrom_data[c] = (positions, coverage, truncations)
        g_cov.append(coverage)
        g_trunc.append(truncations)

    if not chrom_data:
        print("  No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background rate ---
    all_cov   = np.concatenate(g_cov)
    all_trunc = np.concatenate(g_trunc)
    lambda_trunc = estimate_background(all_trunc, all_cov, min_coverage)

    print(f"\n  Background truncation rate (genome-wide):", file=sys.stderr)
    print(f"    λ = {lambda_trunc:.6f}  "
          f"(1/{int(1/lambda_trunc) if lambda_trunc > 0 else '∞'} reads)",
          file=sys.stderr)

    if lambda_trunc <= 0:
        print("  ERROR: background rate is zero — no testable positions.", file=sys.stderr)
        return

    # --- Open output file ---
    out_path = f"{prefix}_truncations.bed"
    total = 0

    with open(out_path, 'w') as fh:
        fh.write(HEADER)

        for c, (positions, coverage, truncations) in chrom_data.items():
            results = test_signal(
                positions, truncations, coverage,
                lambda_trunc, min_coverage, min_fraction, fdr
            )
            write_bed(results, c, 'trunc', fh)
            total += len(results)

    print(f"\n  Significant truncation sites (FDR < {fdr}): {total:,}", file=sys.stderr)
    print(f"  Output: {out_path}", file=sys.stderr)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink cits: call RT truncation sites from a pileup or BAM')

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
        fdr          = args.fdr,
    )
