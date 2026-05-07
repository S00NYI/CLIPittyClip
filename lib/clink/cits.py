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


def build_local_rates(positions: np.ndarray,
                      clean_truncations: np.ndarray,
                      coverage: np.ndarray,
                      max_gap: int = 25) -> np.ndarray:
    """
    Compute a per-position local truncation rate from read clusters.

    A cluster is a maximal run of positions where consecutive entries in
    `positions` are ≤ max_gap bp apart.  Within each cluster:

        local_rate = Σ(clean_truncations) / Σ(coverage)

    This mirrors CTK's approach: truncation enrichment is tested relative to
    the local read density (cluster background) rather than a genome-wide rate.
    Positions in clusters with no truncations get rate=0 and are skipped by
    test_signal's pre-filter.

    Args:
        positions         : 0-based reference positions (sorted, sparse)
        clean_truncations : truncation counts from reads with no deletion in CIGAR
        coverage          : read coverage at each position
        max_gap           : maximum bp gap between positions in same cluster

    Returns:
        Array of per-position local rates, same length as positions.
    """
    n = len(positions)
    if n == 0:
        return np.array([], dtype=np.float64)

    local_rates = np.zeros(n, dtype=np.float64)

    cluster_start = 0
    for i in range(1, n + 1):
        # Cluster ends at gap or end of array
        if i == n or (positions[i] - positions[i - 1]) > max_gap:
            sl = slice(cluster_start, i)
            total_cov   = float(coverage[sl].sum())
            total_trunc = float(clean_truncations[sl].sum())
            if total_cov > 0 and total_trunc > 0:
                local_rates[sl] = total_trunc / total_cov
            # else leave as 0.0 — no signal in cluster, skip in test_signal
            cluster_start = i

    return local_rates


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
             cluster_gap:  int   = 25,
             verbose:      bool  = True):
    """
    End-to-end CITS: pileup (.npz) or BAM → significant truncation sites BED.

    Strand-aware: fwd (+) and rev (-) strand signals are tested independently
    using per-position local cluster background rates (matching CTK).

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
        cluster_gap  : max bp gap between positions in the same read cluster
                       (used for local background rate; default 25, matching CTK)
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
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  min_signal={min_signal}"
          f"  fdr={fdr}  cluster_gap={cluster_gap}bp",
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
            positions, coverage, truncations, deletions, clean_truncations, subs = fwd
            chrom_fwd[c] = (positions, coverage, clean_truncations)
            g_cov.append(coverage)
            g_trunc.append(clean_truncations)
        if rev is not None:
            positions, coverage, truncations, deletions, clean_truncations, subs = rev
            chrom_rev[c] = (positions, coverage, clean_truncations)
            g_cov.append(coverage)
            g_trunc.append(clean_truncations)

    if not chrom_fwd and not chrom_rev:
        print("  No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background rate (reference only — local rates used for testing) ---
    all_cov   = np.concatenate(g_cov)
    all_trunc = np.concatenate(g_trunc)
    lambda_global = estimate_background(all_trunc, all_cov, min_coverage)

    print(f"\n  Global truncation rate (both strands, reference only):",
          file=sys.stderr)
    print(f"    λ_global = {lambda_global:.6f}  "
          f"(1/{int(1/lambda_global) if lambda_global > 0 else '∞'} reads)",
          file=sys.stderr)
    print(f"  Local cluster background enabled  (gap≤{cluster_gap}bp)",
          file=sys.stderr)
    print(f"  Signal: clean truncations (reads without deletion in CIGAR)",
          file=sys.stderr)

    if lambda_global <= 0:
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
                positions, coverage, clean_truncations = chrom_fwd[c]
                local_rates = build_local_rates(
                    positions, clean_truncations, coverage, max_gap=cluster_gap)
                results = test_signal(
                    positions, clean_truncations, coverage,
                    local_rates, min_coverage, min_fraction, fdr,
                    min_signal=min_signal)
                write_bed(results, c, 'trunc', fh, strand='+')
                total_fwd += len(results)

            if c in chrom_rev:
                positions, coverage, clean_truncations = chrom_rev[c]
                local_rates = build_local_rates(
                    positions, clean_truncations, coverage, max_gap=cluster_gap)
                results = test_signal(
                    positions, clean_truncations, coverage,
                    local_rates, min_coverage, min_fraction, fdr,
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
    parser.add_argument('--cluster-gap', type=int, default=25,
        help='Max bp gap between positions in same read cluster for local '
             'background estimation (default: 25, matching CTK)')

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
        cluster_gap  = args.cluster_gap,
    )
