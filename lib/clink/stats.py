#!/usr/bin/env python3
"""
Clink stats.py — Binomial significance testing on pileup signals

Algorithm:
  1. Scan BAM → pileup arrays (via pileup.py)
  2. Estimate genome-wide background rates, excluding top 1% signal positions
       λ_trunc = Σ(truncations) / Σ(coverage)    across background positions
       λ_del   = Σ(deletions)   / Σ(coverage)
       λ_sub   = Σ(X>Y events)  / Σ(coverage)    per substitution type
  3. Vectorized Binomial test at every position with sufficient coverage:
       p = P(X ≥ observed | Binomial(coverage, λ))   via scipy binom.sf
  4. Benjamini-Hochberg FDR correction (no statsmodels needed)
  5. Write BED6+ files per signal type

Output columns:
    chrom  start  end  name  score  strand  signal  coverage  fraction  pvalue  qvalue
    score = min(-log10(qvalue) * 100, 1000)   BED-compatible 0–1000 range
    strand = '+' or '-'  (from read strand; '.' for legacy unstranded pileups)

Usage:
    python stats.py sample.bam --chrom chr1
    python stats.py sample.bam --prefix results/SRR26536505
"""

import sys
import argparse
import numpy as np
from pathlib import Path
from scipy.stats import binom

# ---------------------------------------------------------------------------
# Background estimation
# ---------------------------------------------------------------------------

def estimate_background(signal: np.ndarray,
                        coverage: np.ndarray,
                        min_coverage: int = 5,
                        exclude_top: float = 0.01) -> float:
    """
    Global background rate: λ = Σ(signal) / Σ(coverage).

    Computed only over positions with coverage >= min_coverage.
    The top `exclude_top` fraction of positions by signal fraction are excluded
    to prevent genuine crosslink sites from inflating the background estimate.

    Returns 0.0 if no valid positions exist.
    """
    mask = coverage >= min_coverage
    if mask.sum() == 0:
        return 0.0

    cov = coverage[mask].astype(np.float64)
    sig = signal[mask].astype(np.float64)

    fracs = np.where(cov > 0, sig / cov, 0.0)
    threshold = np.quantile(fracs, 1.0 - exclude_top)

    if threshold <= 0:
        # Sparse signal: >99% of positions have fraction=0 (typical for deletions
        # and rare substitution types). Quantile filter would exclude all signal
        # positions and return rate=0. Use raw rate over all covered positions instead.
        return sig.sum() / cov.sum()

    keep = fracs <= threshold
    if keep.sum() == 0:
        return sig.sum() / cov.sum()

    return sig[keep].sum() / cov[keep].sum()


# ---------------------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ---------------------------------------------------------------------------

def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg adjusted p-values (q-values).

    Input:  raw p-values for tested positions (1D array, no NaN)
    Output: BH-adjusted q-values, same shape

    A position is significant if qvalue <= desired FDR threshold.
    """
    m = len(pvalues)
    if m == 0:
        return np.array([], dtype=np.float64)

    order   = np.argsort(pvalues)
    ranks   = np.arange(1, m + 1, dtype=np.float64)
    adjusted = pvalues[order] * m / ranks

    # Enforce monotonicity right-to-left (take cumulative min from right)
    for i in range(m - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    adjusted = np.minimum(adjusted, 1.0)

    qvalues = np.empty(m, dtype=np.float64)
    qvalues[order] = adjusted
    return qvalues


# ---------------------------------------------------------------------------
# Per-signal significance testing
# ---------------------------------------------------------------------------

def test_signal(positions:        np.ndarray,
                signal:           np.ndarray,
                coverage:         np.ndarray,
                background_rate:  float,
                min_coverage:     int   = 5,
                min_fraction:     float = 0.05,
                fdr_threshold:    float = 0.05,
                min_signal:       int   = 1) -> list:
    """
    Binomial test + BH correction for one signal type.

    Pre-filters:
      - coverage >= min_coverage          (too few reads → unreliable)
      - signal >= min_signal              (minimum raw count of mutation/truncation reads)
      - signal / coverage >= min_fraction (must be above a raw fraction floor)

    Returns list of dicts for significant sites, sorted by q-value then
    fraction descending (best sites first).
    """
    # Accept either scalar or per-position array
    scalar_rate = np.isscalar(background_rate)
    if scalar_rate and background_rate <= 0:
        return []
    if not scalar_rate and not np.any(background_rate > 0):
        return []

    # Pre-filter: positions worth testing
    with np.errstate(divide='ignore', invalid='ignore'):
        fracs = np.where(coverage > 0, signal / coverage.astype(float), 0.0)

    testable = (
        (coverage >= min_coverage) &
        (signal >= max(1, min_signal)) &
        (fracs >= min_fraction)
    )

    n_test = testable.sum()
    if n_test == 0:
        return []

    # Vectorized Binomial survival function: P(X >= k) = binom.sf(k-1, n, p)
    k = signal[testable].astype(int)
    n = coverage[testable].astype(int)
    rate = background_rate if np.isscalar(background_rate) else background_rate[testable]
    # Mask out zero-rate positions (no local cluster signal) — skip test
    if not np.isscalar(rate):
        valid = rate > 0
        if not np.any(valid):
            return []
        raw_pvals      = np.ones(len(k))
        raw_pvals[valid] = binom.sf(k[valid] - 1, n[valid], rate[valid])
    else:
        raw_pvals = binom.sf(k - 1, n, rate)

    # BH correction over all tested positions
    qvals = bh_fdr(raw_pvals)

    # Collect significant sites
    sig_mask = qvals <= fdr_threshold
    testable_indices = np.where(testable)[0]

    results = []
    for j, i in enumerate(testable_indices):
        if not sig_mask[j]:
            continue
        results.append({
            'pos':      int(positions[i]),
            'signal':   int(signal[i]),
            'coverage': int(coverage[i]),
            'fraction': float(fracs[i]),
            'pvalue':   float(raw_pvals[j]),
            'qvalue':   float(qvals[j]),
        })

    results.sort(key=lambda x: (x['qvalue'], -x['fraction']))
    return results


# ---------------------------------------------------------------------------
# BED output
# ---------------------------------------------------------------------------

HEADER = (
    "#chrom\tstart\tend\tname\tscore\tstrand\t"
    "signal\tcoverage\tfraction\tpvalue\tqvalue\n"
)

def write_bed(results: list, chrom: str, signal_label: str, fh,
              strand: str = '.'):
    """
    Write significant sites in BED6 + extra columns.
    score = min(-log10(qvalue) * 100, 1000)  — BED-compatible 0–1000 integer
    strand: '+', '-', or '.' (default '.' for legacy unstranded pileups)
    """
    for r in results:
        pos   = r['pos']
        score = min(int(-np.log10(max(r['qvalue'], 1e-300)) * 100), 1000)
        name  = f"{signal_label}:{chrom}:{pos}:{strand}"
        fh.write(
            f"{chrom}\t{pos}\t{pos+1}\t{name}\t{score}\t{strand}\t"
            f"{r['signal']}\t{r['coverage']}\t{r['fraction']:.4f}\t"
            f"{r['pvalue']:.3e}\t{r['qvalue']:.3e}\n"
        )


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_analysis(bam_path:     str,
                 chrom:        str   = None,
                 min_mapq:     int   = 255,
                 max_nh:       int   = 1,
                 min_coverage: int   = 5,
                 min_fraction: float = 0.05,
                 fdr:          float = 0.05,
                 prefix:       str   = None,
                 verbose:      bool  = True):
    """
    End-to-end: BAM → pileup → background → test → BED files.

    Output files written to:
        {prefix}_truncations.bed
        {prefix}_deletions.bed
        {prefix}_{X}to{Y}.bed    (one per observed substitution type)
    """
    # Import here so stats.py can also be used as a library
    from pileup import scan_bam, to_arrays

    print(f"\nClink  |  {Path(bam_path).name}", file=sys.stderr)
    print(f"  mapq>={min_mapq}  nh<={max_nh}  min_cov={min_coverage}"
          f"  min_frac={min_fraction}  fdr={fdr}", file=sys.stderr)

    # --- Scan BAM ---
    pileups = scan_bam(bam_path, chrom=chrom, min_mapq=min_mapq,
                       max_nh=max_nh, verbose=verbose)

    # --- Convert to arrays, accumulate genome-wide ---
    chrom_data = {}
    g_cov, g_trunc, g_del = [], [], []
    g_subs     = {}   # sub_type -> [signal arrays per chrom]
    g_subs_cov = {}   # sub_type -> [coverage arrays for same chroms only]

    for c, pileup in pileups.items():
        result = to_arrays(pileup)
        if result is None:
            continue
        positions, coverage, truncations, deletions, clean_truncations, subs = result
        chrom_data[c] = (positions, coverage, truncations, deletions, clean_truncations, subs)

        g_cov.append(coverage)
        g_trunc.append(truncations)
        g_del.append(deletions)
        for sub_type, arr in subs.items():
            g_subs.setdefault(sub_type, []).append(arr)
            g_subs_cov.setdefault(sub_type, []).append(coverage)

    if not chrom_data:
        print("No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background estimation ---
    all_cov   = np.concatenate(g_cov)
    all_trunc = np.concatenate(g_trunc)
    all_del   = np.concatenate(g_del)

    λ_trunc = estimate_background(all_trunc, all_cov, min_coverage)
    λ_del   = estimate_background(all_del,   all_cov, min_coverage)
    λ_subs  = {
        st: estimate_background(
                np.concatenate(g_subs[st]),
                np.concatenate(g_subs_cov[st]),
                min_coverage)
        for st in g_subs
    }

    print(f"\n  Background rates (genome-wide):", file=sys.stderr)
    print(f"    truncation : {λ_trunc:.6f}  "
          f"(1/{int(1/λ_trunc) if λ_trunc > 0 else '∞'} reads)", file=sys.stderr)
    print(f"    deletion   : {λ_del:.6f}",   file=sys.stderr)
    for st, lam in sorted(λ_subs.items()):
        print(f"    {st[0]}>{st[1]}         : {lam:.6f}", file=sys.stderr)

    # --- Open output files ---
    if prefix is None:
        stem = Path(bam_path).stem
        prefix = stem

    def open_bed(suffix):
        path = f"{prefix}_{suffix}.bed"
        fh = open(path, 'w')
        fh.write(HEADER)
        return path, fh

    path_trunc, fh_trunc = open_bed("truncations")
    path_del,   fh_del   = open_bed("deletions")

    sub_files = {}
    for sub_type in λ_subs:
        label = f"{sub_type[0]}to{sub_type[1]}"
        sub_files[sub_type] = open_bed(label)

    totals = {'trunc': 0, 'del': 0}
    totals.update({st: 0 for st in λ_subs})

    # --- Per-chromosome testing ---
    for c, (positions, coverage, truncations, deletions, clean_truncations, subs) in chrom_data.items():

        # Truncation sites (CITS-equivalent)
        t_res = test_signal(positions, truncations, coverage,
                            λ_trunc, min_coverage, min_fraction, fdr)
        write_bed(t_res, c, "trunc", fh_trunc)
        totals['trunc'] += len(t_res)

        # Deletion sites (CIMS-equivalent)
        d_res = test_signal(positions, deletions, coverage,
                            λ_del, min_coverage, min_fraction, fdr)
        write_bed(d_res, c, "del", fh_del)
        totals['del'] += len(d_res)

        # Substitution sites (PAR-CLIP T>C, iCLIP any sub)
        for sub_type, sub_arr in subs.items():
            if sub_type not in λ_subs:
                continue
            s_res = test_signal(positions, sub_arr, coverage,
                                λ_subs[sub_type], min_coverage, min_fraction, fdr)
            _, fh_sub = sub_files[sub_type]
            write_bed(s_res, c, f"{sub_type[0]}to{sub_type[1]}", fh_sub)
            totals[sub_type] = totals.get(sub_type, 0) + len(s_res)

    # --- Close files ---
    fh_trunc.close()
    fh_del.close()
    for _, fh in sub_files.values():
        fh.close()

    # --- Summary ---
    print(f"\n  Significant sites (FDR < {fdr}):", file=sys.stderr)
    print(f"    Truncations : {totals['trunc']:>8,}  → {path_trunc}", file=sys.stderr)
    print(f"    Deletions   : {totals['del']:>8,}  → {path_del}",   file=sys.stderr)
    for st in sorted(λ_subs):
        label = f"{st[0]}>{st[1]}"
        path, _ = sub_files[st]
        print(f"    {label:<12}: {totals.get(st,0):>8,}  → {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink stats: Binomial site calling from BAM pileup')
    parser.add_argument('bam',
        help='Sorted, indexed BAM file')
    parser.add_argument('--chrom', default=None,
        help='Restrict to one chromosome (default: all standard chroms)')
    parser.add_argument('--mapq', type=int, default=255,
        help='Minimum MAPQ (default: 255 = STAR unique mappers)')
    parser.add_argument('--nh', type=int, default=1,
        help='Max NH tag — 1 = unique only (default: 1)')
    parser.add_argument('--min-cov', type=int, default=5,
        help='Minimum coverage to test a position (default: 5)')
    parser.add_argument('--min-frac', type=float, default=0.05,
        help='Minimum raw signal fraction pre-filter (default: 0.05)')
    parser.add_argument('--fdr', type=float, default=0.05,
        help='BH FDR threshold (default: 0.05)')
    parser.add_argument('--prefix', default=None,
        help='Output file prefix (default: derived from BAM filename)')
    args = parser.parse_args()

    run_analysis(
        bam_path     = args.bam,
        chrom        = args.chrom,
        min_mapq     = args.mapq,
        max_nh       = args.nh,
        min_coverage = args.min_cov,
        min_fraction = args.min_frac,
        fdr          = args.fdr,
        prefix       = args.prefix,
    )
