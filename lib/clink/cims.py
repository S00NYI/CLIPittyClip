#!/usr/bin/env python3
"""
Clink cims.py — Crosslink-Induced Mutation Sites

Calls statistically significant deletion and substitution sites from a
sorted, deduplicated BAM. Replaces CTK CIMS.pl in the CLIPittyClip pipeline.

Covers:
  - Deletions  : 1-bp CIGAR D operations (classic CIMS signal)
  - Substitutions: all 12 mismatch types (T>C for PAR-CLIP; any sub for iCLIP)
    Use --sub-types to restrict output (e.g. T>C only for PAR-CLIP)

Algorithm:
  1. Scan BAM → per-position deletion + substitution counts and coverage
  2. Estimate genome-wide background rates per signal type
  3. Binomial test + Benjamini-Hochberg FDR correction
  4. Write BED6+ output per signal type

Output:
  {prefix}_deletions.bed
  {prefix}_{X}to{Y}.bed    (one per observed/requested substitution type)

  Columns: chrom  start  end  name  score  strand  signal  coverage  fraction  pvalue  qvalue
  score = min(-log10(qvalue) * 100, 1000)

Usage:
    # All signal types (deletions + all 12 substitutions):
    python cims.py --bam sample_dedup.bam --prefix results/SAMPLE

    # Deletions only:
    python cims.py --bam sample_dedup.bam --prefix results/SAMPLE --no-subs

    # PAR-CLIP (T>C only):
    python cims.py --bam sample_dedup.bam --prefix results/SAMPLE --sub-types TC

    # iCLIP (deletions + all subs):
    python cims.py --bam sample_dedup.bam --prefix results/SAMPLE --chrom chr1
"""

import sys
import argparse
from pathlib import Path

from pileup import scan_bam, to_arrays, load_pileup
from stats import estimate_background, test_signal, write_bed, HEADER

import numpy as np


# Valid substitution type codes (ref+alt, no separator)
ALL_SUB_TYPES = {
    'AC', 'AG', 'AT',
    'CA', 'CG', 'CT',
    'GA', 'GC', 'GT',
    'TA', 'TC', 'TG',
}


def parse_sub_types(raw: str) -> set:
    """
    Parse a comma-separated list of substitution type codes.
    Accepts 'TC', 'T>C', 'TtoC' — all normalised to ('T','C') tuples.
    Returns None if raw is None (= all types).
    """
    if raw is None:
        return None
    types = set()
    for token in raw.split(','):
        token = token.strip().upper().replace('>', '').replace('TO', '')
        if len(token) == 2 and token in ALL_SUB_TYPES:
            types.add((token[0], token[1]))
        else:
            print(f"  WARNING: unrecognised substitution type '{token}' — skipped",
                  file=sys.stderr)
    return types if types else None


def run_cims(pileup_path:  str   = None,
             bam_path:     str   = None,
             prefix:       str   = None,
             chrom:        str   = None,
             min_mapq:     int   = 255,
             max_nh:       int   = 1,
             min_coverage: int   = 5,
             min_fraction: float = 0.05,
             fdr:          float = 0.05,
             no_subs:      bool  = False,
             sub_types:    set   = None,
             verbose:      bool  = True):
    """
    End-to-end CIMS: pileup (.npz) or BAM → significant deletion + substitution BEDs.

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
        min_fraction : minimum raw signal fraction pre-filter
        fdr          : Benjamini-Hochberg FDR threshold
        no_subs      : if True, skip substitution calling entirely
        sub_types    : set of (ref, alt) tuples to output; None = all observed
        verbose      : print progress to stderr
    """
    if pileup_path is None and bam_path is None:
        print("ERROR: supply either --pileup or --bam", file=sys.stderr)
        sys.exit(1)

    src_name = Path(pileup_path or bam_path).name
    if prefix is None:
        stem = Path(pileup_path or bam_path).stem
        prefix = stem.replace('_pileup', '')

    print(f"\nClink cims  |  {src_name}", file=sys.stderr)
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  fdr={fdr}",
          file=sys.stderr)
    if no_subs:
        print(f"  Mode: deletions only", file=sys.stderr)
    elif sub_types:
        labels = [f"{r}>{a}" for r, a in sorted(sub_types)]
        print(f"  Mode: deletions + substitutions {labels}", file=sys.stderr)
    else:
        print(f"  Mode: deletions + all substitution types", file=sys.stderr)

    # --- Load pileup data ---
    if pileup_path:
        chrom_data_full = load_pileup(pileup_path)
    else:
        print(f"  mapq>={min_mapq}  nh<={max_nh}", file=sys.stderr)
        pileups = scan_bam(bam_path, chrom=chrom, min_mapq=min_mapq,
                           max_nh=max_nh, verbose=verbose)
        chrom_data_full = {}
        for c, strand_pileups in pileups.items():
            chrom_data_full[c] = {}
            for s_label, pileup in strand_pileups.items():
                result = to_arrays(pileup)
                if result is None:
                    continue
                positions, coverage, truncations, deletions, subs = result
                chrom_data_full[c][s_label] = (positions, coverage, truncations, deletions, subs)
            if not chrom_data_full[c]:
                del chrom_data_full[c]

    # --- Filter to requested chrom if specified ---
    if chrom:
        chrom_data_full = {c: v for c, v in chrom_data_full.items() if c == chrom}

    # --- Extract deletion + sub signal, accumulate genome-wide (both strands) ---
    # chrom_data: {chrom: {strand: (positions, coverage, deletions, subs)}}
    chrom_data = {}
    g_cov, g_del = [], []
    g_subs     = {}   # sub_type -> [signal arrays per (chrom, strand)]
    g_subs_cov = {}   # sub_type -> [coverage arrays for same (chrom, strand)s]

    for c, strands in chrom_data_full.items():
        chrom_data[c] = {}
        for s_label, arrays in strands.items():
            positions, coverage, truncations, deletions, subs = arrays
            chrom_data[c][s_label] = (positions, coverage, deletions, subs)

            g_cov.append(coverage)
            g_del.append(deletions)
            for sub_type, arr in subs.items():
                if sub_types and sub_type not in sub_types:
                    continue
                g_subs.setdefault(sub_type, []).append(arr)
                g_subs_cov.setdefault(sub_type, []).append(coverage)
        if not chrom_data[c]:
            del chrom_data[c]

    if not chrom_data:
        print("  No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background rates ---
    all_cov = np.concatenate(g_cov)
    all_del = np.concatenate(g_del)

    lambda_del = estimate_background(all_del, all_cov, min_coverage)
    lambda_subs = {}
    if not no_subs:
        lambda_subs = {
            st: estimate_background(
                    np.concatenate(g_subs[st]),
                    np.concatenate(g_subs_cov[st]),
                    min_coverage)
            for st in g_subs
        }

    print(f"\n  Background rates (genome-wide):", file=sys.stderr)
    print(f"    deletion : {lambda_del:.6f}", file=sys.stderr)
    for st, lam in sorted(lambda_subs.items()):
        print(f"    {st[0]}>{st[1]}      : {lam:.6f}", file=sys.stderr)

    # --- Open output files ---
    out_del  = f"{prefix}_deletions.bed"
    sub_files = {}
    if not no_subs:
        for sub_type in lambda_subs:
            label = f"{sub_type[0]}to{sub_type[1]}"
            sub_files[sub_type] = f"{prefix}_{label}.bed"

    totals = {'del': 0}
    totals.update({st: 0 for st in lambda_subs})

    STRAND_CHAR = {'pos': '+', 'neg': '-'}
    with open(out_del, 'w') as fh_del:
        fh_del.write(HEADER)

        sub_fhs = {}
        for st, path in sub_files.items():
            sub_fhs[st] = open(path, 'w')
            sub_fhs[st].write(HEADER)

        for c, strands in chrom_data.items():
            for s_label, (positions, coverage, deletions, subs) in strands.items():
                strand_char = STRAND_CHAR.get(s_label, '.')

                # Deletions
                if lambda_del > 0:
                    d_res = test_signal(positions, deletions, coverage,
                                        lambda_del, min_coverage, min_fraction, fdr)
                    write_bed(d_res, c, 'del', fh_del, strand=strand_char)
                    totals['del'] += len(d_res)

                # Substitutions
                for sub_type, fh_sub in sub_fhs.items():
                    if sub_type not in subs:
                        continue
                    lam = lambda_subs.get(sub_type, 0)
                    if lam <= 0:
                        continue
                    s_res = test_signal(positions, subs[sub_type], coverage,
                                        lam, min_coverage, min_fraction, fdr)
                    write_bed(s_res, c, f"{sub_type[0]}to{sub_type[1]}", fh_sub,
                              strand=strand_char)
                    totals[sub_type] = totals.get(sub_type, 0) + len(s_res)

        for fh in sub_fhs.values():
            fh.close()

    # --- Summary ---
    print(f"\n  Significant sites (FDR < {fdr}):", file=sys.stderr)
    print(f"    Deletions : {totals['del']:>8,}  → {out_del}", file=sys.stderr)
    for st in sorted(lambda_subs):
        label = f"{st[0]}>{st[1]}"
        print(f"    {label:<12}: {totals.get(st,0):>8,}  → {sub_files[st]}",
              file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink cims: call deletion and substitution sites from a pileup or BAM')

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
        help='Minimum raw signal fraction pre-filter (default: 0.05)')
    parser.add_argument('--fdr', type=float, default=0.05,
        help='BH FDR threshold (default: 0.05)')
    parser.add_argument('--no-subs', action='store_true',
        help='Skip substitution calling — output deletions only')
    parser.add_argument('--sub-types', default=None,
        help='Comma-separated substitution types to output '
             '(e.g. "TC" or "T>C" for PAR-CLIP, "AG,TC" for multiple). '
             'Default: all observed types.')

    args = parser.parse_args()

    requested_subs = parse_sub_types(args.sub_types) if args.sub_types else None

    run_cims(
        pileup_path  = args.pileup,
        bam_path     = args.bam,
        prefix       = args.prefix,
        chrom        = args.chrom,
        min_mapq     = args.mapq,
        max_nh       = args.nh,
        min_coverage = args.min_cov,
        min_fraction = args.min_frac,
        fdr          = args.fdr,
        no_subs      = args.no_subs,
        sub_types    = requested_subs,
    )
