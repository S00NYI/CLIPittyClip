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
  1. Read the stranded pileup (pileup.py) one chromosome and strand at a time
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
import math
import argparse
from pathlib import Path

from pileup import scan_bam, to_arrays
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
             min_signal:   int   = 1,
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
        min_signal   : minimum raw mutation/deletion read count at a position
                       (equivalent to CTK CIMS.pl -m; default 1 = no extra filter)
        fdr          : Benjamini-Hochberg FDR threshold
        no_subs      : if True, skip substitution calling entirely
        sub_types    : set of (ref, alt) tuples to output; None = all observed
        verbose      : print progress to stderr
    """
    if pileup_path is None and bam_path is None:
        print("ERROR: supply either --pileup or --bam", file=sys.stderr)
        sys.exit(1)

    # A pileup .npz is read one chromosome and strand at a time.
    if pileup_path:
        return run_cims_by_chrom(pileup_path, prefix=prefix, chrom=chrom,
                                 min_coverage=min_coverage,
                                 min_fraction=min_fraction,
                                 min_signal=min_signal, fdr=fdr,
                                 no_subs=no_subs, sub_types=sub_types)

    src_name = Path(pileup_path or bam_path).name
    if prefix is None:
        stem = Path(pileup_path or bam_path).stem
        prefix = stem.replace('_pileup', '')

    print(f"\nClink cims  |  {src_name}", file=sys.stderr)
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  min_signal={min_signal}  fdr={fdr}",
          file=sys.stderr)
    if no_subs:
        print(f"  Mode: deletions only", file=sys.stderr)
    elif sub_types:
        labels = [f"{r}>{a}" for r, a in sorted(sub_types)]
        print(f"  Mode: deletions + substitutions {labels}", file=sys.stderr)
    else:
        print(f"  Mode: deletions + all substitution types", file=sys.stderr)

    # --- Scan the BAM (pileup input returned above) ---
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
            positions, coverage, truncations, deletions, clean_truncations, subs = strand_arr

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
                                        lambda_del, min_coverage, min_fraction, fdr,
                                        min_signal=min_signal)
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
                                        lam, min_coverage, min_fraction, fdr,
                                        min_signal=min_signal)
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
# Pileup input: one chromosome and strand at a time
#
# Loading the whole pileup and calling estimate_background() on genome-wide
# arrays needs about 66 bytes per position (22.7 GB for a 352M-position pig
# group pileup), which ran out of RAM on a human group pileup. This code reads
# one chromosome and strand at a time from the .npz (0.8 GB on the same pig
# pileup, same run time) and gives byte-identical output:
#   - The FDR step already works per chromosome and strand.
#   - The background rate is rebuilt exactly from small per-chromosome parts:
#     integer sums, the count of zero-fraction positions, and the non-zero
#     fractions. The quantile uses the same interpolation as np.quantile
#     (method 'linear').
# ---------------------------------------------------------------------------

def _npz_layout(data):
    """
    Map the keys of an open pileup .npz, in the same order as load_pileup().

    Returns [(chrom, strand_key, key_prefix, {sub_type: key})], with chroms
    sorted and 'fwd' before 'rev'. An old unstranded file is read as fwd only.
    """
    files = data.files
    is_stranded = any(('__fwd__' in k or '__rev__' in k) for k in files)
    layout = []
    if not is_stranded:
        chroms = sorted(set(k.split('__')[0] for k in files))
        for chrom in chroms:
            pfx = chrom
            sub_prefix = f'{pfx}__sub__'
            subs = {(k[len(sub_prefix):][0], k[len(sub_prefix):][1]): k
                    for k in files if k.startswith(sub_prefix)}
            layout.append((chrom, 'fwd', pfx, subs))
        return layout

    chroms = set()
    for key in files:
        parts = key.split('__')
        if len(parts) >= 3 and parts[1] in ('fwd', 'rev'):
            chroms.add(parts[0])
    for chrom in sorted(chroms):
        for strand_key in ('fwd', 'rev'):
            pfx = f'{chrom}__{strand_key}'
            if f'{pfx}__positions' not in files:
                continue
            sub_prefix = f'{pfx}__sub__'
            subs = {(k[len(sub_prefix):][0], k[len(sub_prefix):][1]): k
                    for k in files if k.startswith(sub_prefix)}
            layout.append((chrom, strand_key, pfx, subs))
    return layout


def _linear_quantile(n: int, n_zero: int, positives: np.ndarray, q: float) -> float:
    """
    np.quantile(fracs, q) (method 'linear') for a virtual sorted array of
    n values: n_zero zeros, then the sorted values of `positives`.
    Follows numpy's virtual index, neighbour indexes and _lerp steps.
    """
    v = (n - 1) * q
    if v >= n - 1:
        lo = hi = n - 1
    elif v < 0:
        lo = hi = 0
    else:
        lo = int(math.floor(v))
        hi = lo + 1
    # numpy sets both indexes to -1 (the last value) above the upper bound,
    # so gamma is v - (-1). a == b there, so the result is the last value.
    gamma = v - (lo if v < n - 1 else -1)

    need = sorted({i - n_zero for i in (lo, hi) if i >= n_zero})
    part = {}
    if need:
        p = np.partition(positives, need)
        part = {k: float(p[k]) for k in need}
    a = np.float64(0.0 if lo < n_zero else part[lo - n_zero])
    b = np.float64(0.0 if hi < n_zero else part[hi - n_zero])

    diff_b_a = b - a
    result = a + diff_b_a * gamma
    if gamma >= 0.5:
        result = b - diff_b_a * (1 - gamma)
    return result


def _background_by_chrom(data, layout, signal_of, min_coverage: int,
                         exclude_top: float = 0.01):
    """
    Same value as estimate_background(concat(signal), concat(coverage)),
    computed one chromosome and strand at a time.

    signal_of(pfx, subs) returns the .npz key of the signal array for that
    chromosome and strand, or None to skip it (as run_cims() skips a strand
    that has no array for a substitution type).
    """
    n = 0
    n_zero = 0
    sum_sig = 0
    sum_cov = 0
    pos_frac, pos_sig, pos_cov = [], [], []

    for _chrom, _strand, pfx, subs in layout:
        sig_key = signal_of(pfx, subs)
        if sig_key is None:
            continue
        coverage = data[f'{pfx}__coverage']
        signal = data[sig_key]
        mask = coverage >= min_coverage
        m = int(mask.sum())
        if m == 0:
            continue
        cov = coverage[mask].astype(np.float64)
        sig = signal[mask].astype(np.float64)
        with np.errstate(divide='ignore', invalid='ignore'):
            fracs = np.where(cov > 0, sig / cov, 0.0)
        n += m
        sum_sig += int(signal[mask].astype(np.int64).sum())
        sum_cov += int(coverage[mask].astype(np.int64).sum())
        nz = fracs > 0
        n_zero += m - int(nz.sum())
        pos_frac.append(fracs[nz])
        pos_sig.append(signal[mask][nz].astype(np.int64))
        pos_cov.append(coverage[mask][nz].astype(np.int64))

    if n == 0:
        return 0.0

    positives = np.concatenate(pos_frac) if pos_frac else np.empty(0)
    threshold = _linear_quantile(n, n_zero, positives, 1.0 - exclude_top)

    raw_rate = np.float64(sum_sig) / np.float64(sum_cov)
    if threshold <= 0:
        return raw_rate

    # keep = fracs <= threshold. Zero-fraction positions are always kept,
    # so only the positive fractions above the threshold are removed.
    over = positives > threshold
    keep_count = n - int(over.sum())
    if keep_count == 0:
        return raw_rate
    keep_sig = sum_sig - int(np.concatenate(pos_sig)[over].sum())
    keep_cov = sum_cov - int(np.concatenate(pos_cov)[over].sum())
    return np.float64(keep_sig) / np.float64(keep_cov)


def run_cims_by_chrom(pileup_path:  str,
                      prefix:       str   = None,
                      chrom:        str   = None,
                      min_coverage: int   = 5,
                      min_fraction: float = 0.05,
                      min_signal:   int   = 1,
                      fdr:          float = 0.05,
                      no_subs:      bool  = False,
                      sub_types:    set   = None):
    """
    CIMS for a pileup .npz, one chromosome and strand at a time.
    run_cims() calls this for pileup input. Memory is about one chromosome
    and strand at a time.
    """
    src_name = Path(pileup_path).name
    if prefix is None:
        prefix = Path(pileup_path).stem.replace('_pileup', '')

    print(f"\nClink cims  |  {src_name}", file=sys.stderr)
    print(f"  min_cov={min_coverage}  min_frac={min_fraction}  min_signal={min_signal}  fdr={fdr}",
          file=sys.stderr)
    if no_subs:
        print(f"  Mode: deletions only", file=sys.stderr)
    elif sub_types:
        labels = [f"{r}>{a}" for r, a in sorted(sub_types)]
        print(f"  Mode: deletions + substitutions {labels}", file=sys.stderr)
    else:
        print(f"  Mode: deletions + all substitution types", file=sys.stderr)

    data = np.load(pileup_path)
    layout = _npz_layout(data)
    if chrom:
        layout = [entry for entry in layout if entry[0] == chrom]
    print(f"  Pileup opened ← {pileup_path}  "
          f"({len(set(e[0] for e in layout))} chromosomes)", file=sys.stderr)

    if not layout:
        print("  No signal found.", file=sys.stderr)
        return

    # Substitution types in first-seen order, as run_cims() builds g_subs.
    seen_subs = []
    for _c, _s, _pfx, subs in layout:
        for st in subs:
            if sub_types and st not in sub_types:
                continue
            if st not in seen_subs:
                seen_subs.append(st)

    # --- Genome-wide background rates (both strands pooled) ---
    lambda_del = _background_by_chrom(
        data, layout, lambda pfx, subs: f'{pfx}__deletions', min_coverage)
    lambda_subs = {}
    if not no_subs:
        for st in seen_subs:
            lambda_subs[st] = _background_by_chrom(
                data, layout, lambda pfx, subs, st=st: subs.get(st), min_coverage)

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

    with open(out_del, 'w') as fh_del:
        fh_del.write(HEADER)

        sub_fhs = {}
        for st, path in sub_files.items():
            sub_fhs[st] = open(path, 'w')
            sub_fhs[st].write(HEADER)

        for c, strand_key, pfx, subs in layout:
            strand_char = '+' if strand_key == 'fwd' else '-'
            positions = data[f'{pfx}__positions']
            coverage  = data[f'{pfx}__coverage']

            # Deletions
            if lambda_del > 0:
                deletions = data[f'{pfx}__deletions']
                d_res = test_signal(positions, deletions, coverage,
                                    lambda_del, min_coverage, min_fraction, fdr,
                                    min_signal=min_signal)
                write_bed(d_res, c, 'del', fh_del, strand=strand_char)
                totals['del'] += len(d_res)
                del deletions

            # Substitutions
            for sub_type, fh_sub in sub_fhs.items():
                if sub_type not in subs:
                    continue
                lam = lambda_subs.get(sub_type, 0)
                if lam <= 0:
                    continue
                s_res = test_signal(positions, data[subs[sub_type]], coverage,
                                    lam, min_coverage, min_fraction, fdr,
                                    min_signal=min_signal)
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
    parser.add_argument('--min-signal', type=int, default=1,
        help='Minimum raw mutation/deletion read count at a position '
             '(equivalent to CTK CIMS.pl -m; default: 1 = no extra filter)')
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
        min_signal   = args.min_signal,
        fdr          = args.fdr,
        no_subs      = args.no_subs,
        sub_types    = requested_subs,
    )
