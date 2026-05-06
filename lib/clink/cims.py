#!/usr/bin/env python3
"""
Clink cims.py — Crosslink-Induced Mutation Sites (strand-aware)

Calls statistically significant deletion and substitution sites from a
pileup.npz or BAM. Processes forward (+) and reverse (-) strand independently.

Covers:
  - Deletions  : 1-bp CIGAR D operations (classic CIMS signal)
  - Substitutions: all 12 mismatch types (T>C for PAR-CLIP; any sub for iCLIP)
    Use --sub-types to restrict output (e.g. T>C only for PAR-CLIP)

Algorithm:
  1. Load stranded pileup (pileup.py) — separate arrays per strand
  2. Estimate genome-wide background rates per signal type (both strands pooled)
  3. Binomial test + Benjamini-Hochberg FDR correction per strand
  4. Write BED6+ output per signal type with correct strand column (+/-)

Output:
  {prefix}_deletions.bed
  {prefix}_{X}to{Y}.bed    (one per observed/requested substitution type)

  Columns: chrom  start  end  name  score  strand  signal  coverage  fraction  pvalue  qvalue
  score = min(-log10(qvalue) * 100, 1000)
  strand = '+' (fwd reads) or '-' (rev reads)

Usage:
    # All signal types (deletions + all 12 substitutions):
    python cims.py --pileup sample_pileup.npz --prefix results/SAMPLE

    # Deletions only:
    python cims.py --pileup sample_pileup.npz --prefix results/SAMPLE --no-subs

    # PAR-CLIP (T>C only):
    python cims.py --pileup sample_pileup.npz --prefix results/SAMPLE --sub-types TC

    # iCLIP (deletions + all subs):
    python cims.py --pileup sample_pileup.npz --prefix results/SAMPLE --chrom chr1
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

    Strand-aware: fwd (+) and rev (-) strand signals are tested independently
    using shared genome-wide background rates (computed across both strands).

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
        for c, pileup in pileups.items():
            result = to_arrays(pileup)
            if result is None:
                continue
            chrom_data_full[c] = result

    # --- Filter to requested chrom if specified ---
    if chrom:
        chrom_data_full = {c: v for c, v in chrom_data_full.items() if c == chrom}

    # --- Build strand-separated data, accumulate genome-wide for background ---
    # chrom_fwd/chrom_rev: {chrom: (positions, coverage, deletions, subs)}
    chrom_fwd = {}
    chrom_rev = {}
    g_cov = []
    g_del = []
    g_subs     = {}   # sub_type -> [signal arrays]
    g_subs_cov = {}   # sub_type -> [coverage arrays]

    for c, strands in chrom_data_full.items():
        for strand_key, strand_char in (('fwd', '+'), ('rev', '-')):
            strand_arr = strands.get(strand_key)
            if strand_arr is None:
                continue
            positions, coverage, truncations, deletions, subs = strand_arr

            store = chrom_fwd if strand_key == 'fwd' else chrom_rev
            store[c] = (positions, coverage, deletions, subs)

            g_cov.append(coverage)
            g_del.append(deletions)
            for sub_type, arr in subs.items():
                if sub_types and sub_type not in sub_types:
                    continue
                g_subs.setdefault(sub_type, []).append(arr)
                g_subs_cov.setdefault(sub_type, []).append(coverage)

    if not chrom_fwd and not chrom_rev:
        print("  No signal found.", file=sys.stderr)
        return

    # --- Genome-wide background rates (both strands pooled) ---
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

    print(f"\n  Background rates (genome-wide, both strands):", file=sys.stderr)
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

    all_chroms = sorted(set(list(chrom_fwd.keys()) + list(chrom_rev.keys())))

    with open(out_del, 'w') as fh_del:
        fh_del.write(HEADER)

        sub_fhs = {}
        for st, path in sub_files.items():
            sub_fhs[st] = open(path, 'w')
            sub_fhs[st].write(HEADER)

        for c in all_chroms:
            for strand_dict, strand_char in ((chrom_fwd, '+'), (chrom_rev, '-')):
                if c not in strand_dict:
                    continue
                positions, coverage, deletions, subs = strand_dict[c]

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
                    write_bed(s_res, c, f"{sub_type[0]}to{sub_type[1]}",
                              fh_sub, strand=strand_char)
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
        description='Clink cims: call deletion and substitution sites (strand-aware) from a pileup or BAM')

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
