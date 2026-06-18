#!/usr/bin/env python3
"""
Clink extract_crosslinks.py — Export all crosslink positions from a pileup

Writes a BED6 file of every position in the pileup with at least one
truncation event. No significance filtering is applied. This is the
unfiltered crosslink file required by tools such as PEKA.

Output:
  {prefix}_all_crosslinks.bed
    chrom  start  end  name  score  strand
    score = truncation count (capped at 1000 for BED compatibility)

Usage:
    python extract_crosslinks.py --pileup sample_pileup.npz --prefix results/SAMPLE
"""

import sys
import argparse
from pathlib import Path

from pileup import load_pileup

import numpy as np


def extract_crosslinks(pileup_path: str, prefix: str = None, min_truncations: int = 1):
    """
    Extract all crosslink positions from a pileup .npz.

    Args:
        pileup_path     : .npz file from clink pileup
        prefix          : output file prefix (default: derived from input name)
        min_truncations : minimum truncation count to include (default: 1)
    """
    src_name = Path(pileup_path).name
    if prefix is None:
        stem = Path(pileup_path).stem
        prefix = stem.replace('_pileup', '')

    print(f"\nClink extract_crosslinks  |  {src_name}", file=sys.stderr)

    chrom_data = load_pileup(pileup_path)

    out_path = f"{prefix}_all_crosslinks.bed"
    total_fwd = total_rev = 0

    with open(out_path, 'w') as fh:
        for c in sorted(chrom_data.keys()):
            strands = chrom_data[c]
            for strand_key, strand_char in [('fwd', '+'), ('rev', '-')]:
                strand_data = strands.get(strand_key)
                if strand_data is None:
                    continue
                positions, coverage, truncations, deletions, clean_truncations, subs = strand_data
                mask = truncations >= min_truncations
                sel_pos = positions[mask]
                sel_trunc = truncations[mask]
                for pos, trunc in zip(sel_pos, sel_trunc):
                    score = min(int(trunc), 1000)
                    name = f"{c}:{pos}:{strand_char}"
                    fh.write(f"{c}\t{pos}\t{pos + 1}\t{name}\t{score}\t{strand_char}\n")
                if strand_char == '+':
                    total_fwd += int(mask.sum())
                else:
                    total_rev += int(mask.sum())

    total = total_fwd + total_rev
    print(f"  Crosslink sites (truncations >= {min_truncations}): {total:,}  "
          f"(+: {total_fwd:,}  -: {total_rev:,})", file=sys.stderr)
    print(f"  Output: {out_path}", file=sys.stderr)
    return out_path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink extract_crosslinks: export all crosslink positions from a pileup')

    parser.add_argument('--pileup', required=True,
        help='Pre-computed pileup .npz (output of clink pileup)')
    parser.add_argument('--prefix', default=None,
        help='Output file prefix (default: derived from input filename)')
    parser.add_argument('--min-truncations', type=int, default=1,
        help='Minimum truncation count to include a position (default: 1)')

    args = parser.parse_args()

    extract_crosslinks(
        pileup_path=args.pileup,
        prefix=args.prefix,
        min_truncations=args.min_truncations,
    )
