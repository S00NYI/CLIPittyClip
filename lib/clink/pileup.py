#!/usr/bin/env python3
"""
Clink pileup.py — Single-pass BAM signal extractor

Scans a sorted, indexed BAM and accumulates four per-position signals:

  coverage     : reads spanning each position
  truncations  : 5' read ends (RT stop sites — CITS signal)
  deletions    : 1-bp CIGAR D operations (crosslink-induced deletions — CIMS signal)
  substitutions: mismatches by type, e.g. T>C (PAR-CLIP), any sub (iCLIP)

Key design decisions:
  - One forward pass per chromosome — no random BAM access
  - Sparse storage (defaultdict) — only positions with signal are stored
  - MD tag used to distinguish CIGAR D (deletion, ref_base present in MD)
    from CIGAR N (intron skip, ref_base=None) without needing a reference FASTA
  - Deletions count toward coverage (read spans the position, base is absent)
  - Intron-spanning positions (N) are not counted as coverage or deletions

Usage:
    python pileup.py sample.bam --chrom chr1
    python pileup.py sample.bam --chrom chrX --mapq 255
"""

import sys
import re
import argparse
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pysam


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASES = {'A', 'C', 'G', 'T'}

# Standard chromosomes: chr1-22, X, Y, M
def is_standard_chrom(name: str) -> bool:
    if name in ('chrX', 'chrY', 'chrM', 'chrMT'):
        return True
    if name.startswith('chr') and name[3:].isdigit():
        return 1 <= int(name[3:]) <= 22
    return False


# ---------------------------------------------------------------------------
# Core data structure
# ---------------------------------------------------------------------------

@dataclass
class ChromPileup:
    """
    Sparse per-position signal accumulator for one chromosome.

    All four dicts are keyed by 0-based reference position.
    subs is keyed by (ref_base, alt_base) tuples e.g. ('T', 'C').
    """
    chrom: str
    coverage:    Dict[int, int] = field(default_factory=lambda: defaultdict(int))
    truncations: Dict[int, int] = field(default_factory=lambda: defaultdict(int))
    deletions:   Dict[int, int] = field(default_factory=lambda: defaultdict(int))
    subs:        Dict[int, Counter] = field(
                     default_factory=lambda: defaultdict(Counter))
    n_reads:   int = 0
    n_skipped: int = 0


# ---------------------------------------------------------------------------
# Per-read signal extraction
# ---------------------------------------------------------------------------

def truncation_site(read: pysam.AlignedSegment) -> int:
    """
    Truncation site = last nucleotide synthesized by RT (0-based reference coords).

    iCLIP library structure:
      - RT reads 3'→5' on RNA template, stops AT the crosslinked nucleotide
      - The last base incorporated is the one immediately 3' of the crosslink on RNA
      - The cDNA 3' end = crosslink site; adapter is ligated here; Read 1 begins here
      - The crosslink nucleotide is therefore ONE POSITION UPSTREAM of the read start

    In reference coordinates (matching CTK/CITS convention):
      Forward (is_reverse=False): crosslink = reference_start - 1
      Reverse (is_reverse=True):  crosslink = reference_end
        (reference_end is pysam's exclusive end; the 5' start of the read is
         reference_end - 1, so one upstream = reference_end)

    This matches CTK's bedExt.pl -n up -l -1 -r -1 convention used by CITS.pl.
    Validated: Clink vs CTK CITS on SRR26536514 STAR BAM shows systematic +1 offset
    on - strand and -1 offset on + strand when using read-start directly; this
    1bp shift corrects it (7 → ~80 exact matches).

    frac = truncations[X] / coverage[X]:
      1.0 → every read spanning position X stopped there (all truncated, none read through)
      0.x → mix of truncated and read-through reads
    """
    return read.reference_end if read.is_reverse else read.reference_start - 1


def extract_signals(read: pysam.AlignedSegment, pileup: ChromPileup) -> None:
    """
    Accumulate all four signals from one read into pileup.

    Uses get_aligned_pairs(with_seq=True) which decodes the MD tag to
    supply reference bases. This lets us:
      - Distinguish D (deletion, MD supplies ref_base) from
        N (intron, MD is absent → ref_base=None)
      - Identify substitutions without a reference FASTA

    Signal rules:
      Intron (N):     query_pos=None, ref_base=None  → skip entirely
      Deletion (D):   query_pos=None, ref_base=set   → +deletion, +coverage
      Insertion (I):  ref_pos=None                   → skip (no ref position)
      Match/mismatch: both set                        → +coverage, +sub if mismatch
    """
    # --- Truncation ---
    # The truncation site is one position upstream of the read's 5' start
    # (matching CTK convention: the last base synthesized by RT = the crosslink site).
    # Because the truncated read starts just AFTER this position, it does not
    # contribute to coverage[trunc_pos] through the normal aligned-pairs loop.
    # We manually add +1 here so the binomial denominator (coverage) includes
    # both truncated reads AND read-through reads at the crosslink position.
    t_pos = truncation_site(read)
    pileup.truncations[t_pos] += 1
    pileup.coverage[t_pos]    += 1  # truncated read counts toward its own crosslink site

    # --- Coverage / deletions / substitutions via aligned pairs ---
    query_seq = read.query_sequence
    if query_seq is None:
        return  # hard-clipped read with no sequence stored

    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except Exception:
        # MD tag absent or malformed — skip substitution/deletion parsing
        # Still count truncation (already done above)
        pileup.n_skipped += 1
        return

    for query_pos, ref_pos, ref_base in pairs:

        # Insertion: read has extra base(s), consumes no reference position
        if ref_pos is None:
            continue

        # Intron skip (N): ref_pos set but ref_base is None (not in MD tag)
        if query_pos is None and ref_base is None:
            continue

        # Deletion (D): ref_pos set, ref_base set (from MD ^X), query_pos None
        if query_pos is None:
            pileup.deletions[ref_pos] += 1
            pileup.coverage[ref_pos]  += 1  # read spans this position
            continue

        # Aligned base: count coverage
        pileup.coverage[ref_pos] += 1

        # Substitution: ref_base (from MD) differs from read base
        ref_upper = ref_base.upper() if ref_base else None
        alt_base  = query_seq[query_pos].upper()

        if ref_upper in BASES and alt_base in BASES and ref_upper != alt_base:
            pileup.subs[ref_pos][(ref_upper, alt_base)] += 1


# ---------------------------------------------------------------------------
# BAM scan
# ---------------------------------------------------------------------------

def scan_bam(
    bam_path: str,
    chrom:    Optional[str] = None,
    min_mapq: int = 20,
    max_nh:   int = 1,
    verbose:  bool = True,
) -> Dict[str, Dict[str, 'ChromPileup']]:
    """
    One-pass BAM scan. Returns {chrom: {'pos': ChromPileup, 'neg': ChromPileup}}.

    Reads are routed to the appropriate strand-specific ChromPileup by
    read.is_reverse. Downstream callers iterate both strand keys to accumulate
    genome-wide background rates and write strand-labelled BED output.

    Args:
        bam_path : sorted, indexed BAM
        chrom    : scan only this chromosome (None = all standard chroms)
        min_mapq : minimum mapping quality filter (STAR unique = 255, BT2 varies)
        max_nh   : max NH tag value; 1 = unique mappers only
        verbose  : print per-chromosome summary to stderr
    """
    pileups: Dict[str, Dict[str, ChromPileup]] = {}

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chroms_to_scan = (
            [chrom] if chrom
            else [sq['SN'] for sq in bam.header['SQ']
                  if is_standard_chrom(sq['SN'])]
        )

        for c in chroms_to_scan:
            p_pos = ChromPileup(chrom=c)
            p_neg = ChromPileup(chrom=c)
            pileups[c] = {'pos': p_pos, 'neg': p_neg}

            for read in bam.fetch(c):

                # Skip unmapped / secondary / supplementary
                # Rejects go to p_pos.n_skipped (QC counter only, not signal)
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    p_pos.n_skipped += 1
                    continue

                # Mapping quality filter
                if read.mapping_quality < min_mapq:
                    p_pos.n_skipped += 1
                    continue

                # Multimapper filter via NH tag
                try:
                    if read.get_tag('NH') > max_nh:
                        p_pos.n_skipped += 1
                        continue
                except KeyError:
                    pass  # no NH tag — allow through

                # No CIGAR (shouldn't happen, but guard)
                if read.cigartuples is None:
                    p_pos.n_skipped += 1
                    continue

                # Route to strand-specific pileup
                strand_pileup = p_neg if read.is_reverse else p_pos
                strand_pileup.n_reads += 1
                extract_signals(read, strand_pileup)

            if verbose:
                for s_label, sp in (('pos', p_pos), ('neg', p_neg)):
                    print(
                        f"  {c} [{s_label}]: {sp.n_reads:>10,} reads | "
                        f"{len(sp.coverage):>10,} covered positions | "
                        f"{sum(sp.truncations.values()):>8,} truncation events | "
                        f"{sum(sp.deletions.values()):>8,} deletion events | "
                        f"{sp.n_skipped:>8,} skipped",
                        file=sys.stderr,
                    )

    return pileups


# ---------------------------------------------------------------------------
# Convert to numpy arrays for stats.py
# ---------------------------------------------------------------------------

def to_arrays(pileup: ChromPileup):
    """
    Convert sparse ChromPileup dicts to sorted numpy arrays.

    Returns:
        positions    : int32  [N]      — 0-based reference positions
        coverage     : uint32 [N]
        truncations  : uint32 [N]
        deletions    : uint32 [N]
        subs         : dict {(ref, alt): uint16 [N]} per substitution type
                       e.g. {('T','C'): array([0,0,3,0,...])}

    Returns None if the pileup is empty.
    """
    all_pos = (set(pileup.coverage) | set(pileup.truncations) |
               set(pileup.deletions) | set(pileup.subs))
    if not all_pos:
        return None

    positions = np.array(sorted(all_pos), dtype=np.int32)
    n = len(positions)
    idx = {p: i for i, p in enumerate(positions)}

    coverage    = np.zeros(n, dtype=np.uint32)
    truncations = np.zeros(n, dtype=np.uint32)
    deletions   = np.zeros(n, dtype=np.uint32)

    for pos, v in pileup.coverage.items():    coverage[idx[pos]]    = v
    for pos, v in pileup.truncations.items(): truncations[idx[pos]] = v
    for pos, v in pileup.deletions.items():   deletions[idx[pos]]   = v

    # Collect all observed substitution types across all positions
    all_sub_types = set()
    for counter in pileup.subs.values():
        all_sub_types.update(counter.keys())

    subs = {}
    for sub_type in sorted(all_sub_types):
        arr = np.zeros(n, dtype=np.uint16)
        for pos, counter in pileup.subs.items():
            if sub_type in counter:
                arr[idx[pos]] = min(counter[sub_type], 65535)
        subs[sub_type] = arr

    return positions, coverage, truncations, deletions, subs


# ---------------------------------------------------------------------------
# Save / load pileup arrays (intermediate .npz cache)
# ---------------------------------------------------------------------------

def save_pileup(chrom_data: dict, path: str) -> None:
    """
    Save pileup arrays for all chromosomes to a compressed .npz file.

    Key naming convention (triple-segment, double-underscore separator):
        {chrom}__{strand}__positions        e.g. chr1__pos__positions
        {chrom}__{strand}__coverage
        {chrom}__{strand}__truncations
        {chrom}__{strand}__deletions
        {chrom}__{strand}__sub__{ref}{alt}  e.g. chr1__pos__sub__TC

    strand is 'pos' (forward, +) or 'neg' (reverse, -).

    chrom_data format:
        {chrom: {'pos': (positions, coverage, truncations, deletions, subs),
                 'neg': (positions, coverage, truncations, deletions, subs)}}
    """
    arrays = {}
    for chrom, strands in chrom_data.items():
        for s_label, (positions, coverage, truncations, deletions, subs) in strands.items():
            pfx = f'{chrom}__{s_label}'
            arrays[f'{pfx}__positions']   = positions
            arrays[f'{pfx}__coverage']    = coverage
            arrays[f'{pfx}__truncations'] = truncations
            arrays[f'{pfx}__deletions']   = deletions
            for (ref, alt), arr in subs.items():
                arrays[f'{pfx}__sub__{ref}{alt}'] = arr

    np.savez_compressed(path, **arrays)
    print(f"  Pileup saved → {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)


def load_pileup(path: str) -> dict:
    """
    Load pileup arrays from a .npz file written by save_pileup().

    Returns:
        {chrom: {'pos': (positions, coverage, truncations, deletions, subs),
                 'neg': (positions, coverage, truncations, deletions, subs)}}
    """
    data = np.load(path)

    # Keys: {chrom}__{strand}__{field}  e.g. chr1__pos__positions
    # Discover (chrom, strand) pairs from first two segments
    chrom_strands: set = set()
    for key in data.files:
        parts = key.split('__')
        if len(parts) >= 2:
            chrom_strands.add((parts[0], parts[1]))

    chrom_data: dict = {}
    for chrom, s_label in sorted(chrom_strands):
        pfx = f'{chrom}__{s_label}'
        positions   = data[f'{pfx}__positions']
        coverage    = data[f'{pfx}__coverage']
        truncations = data[f'{pfx}__truncations']
        deletions   = data[f'{pfx}__deletions']

        subs = {}
        sub_pfx = f'{pfx}__sub__'
        for key in data.files:
            if key.startswith(sub_pfx):
                code = key[len(sub_pfx):]   # e.g. 'TC'
                subs[(code[0], code[1])] = data[key]

        if chrom not in chrom_data:
            chrom_data[chrom] = {}
        chrom_data[chrom][s_label] = (positions, coverage, truncations, deletions, subs)

    print(f"  Pileup loaded ← {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)
    return chrom_data


# ---------------------------------------------------------------------------
# Quick summary report (for validation / sanity-checking)
# ---------------------------------------------------------------------------

def print_summary(chrom: str, positions, coverage, truncations, deletions, subs,
                  top_n: int = 10):
    """Print a human-readable summary of the pileup arrays."""
    print(f"\n{'='*60}")
    print(f"  {chrom}  —  {len(positions):,} positions with signal")
    print(f"{'='*60}")
    print(f"  Max coverage:       {coverage.max():>10,}")
    print(f"  Total truncations:  {truncations.sum():>10,}")
    print(f"  Total deletions:    {deletions.sum():>10,}")
    print(f"  Substitution types: {[f'{r}>{a}' for r,a in sorted(subs.keys())]}")

    if len(subs) > 0:
        for (ref, alt), arr in sorted(subs.items()):
            print(f"    {ref}>{alt}: {arr.sum():,} total events")

    # Top truncation sites
    if truncations.sum() > 0:
        print(f"\n  Top {top_n} truncation sites (raw count):")
        top = np.argsort(truncations)[-top_n:][::-1]
        for i in top:
            if truncations[i] > 0:
                cov = coverage[i] if coverage[i] > 0 else 1
                frac = truncations[i] / cov
                print(f"    pos={positions[i]:>12,}  "
                      f"trunc={truncations[i]:>6,}  "
                      f"cov={cov:>8,}  "
                      f"frac={frac:.3f}")

    # Top deletion sites
    if deletions.sum() > 0:
        print(f"\n  Top {top_n} deletion sites (raw count):")
        top = np.argsort(deletions)[-top_n:][::-1]
        for i in top:
            if deletions[i] > 0:
                cov = coverage[i] if coverage[i] > 0 else 1
                frac = deletions[i] / cov
                print(f"    pos={positions[i]:>12,}  "
                      f"del={deletions[i]:>6,}  "
                      f"cov={cov:>8,}  "
                      f"frac={frac:.3f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink pileup: scan BAM and save per-position signal arrays')
    parser.add_argument('bam',
        help='Sorted, indexed BAM file (output of clink collapse)')
    parser.add_argument('--out', default=None,
        help='Output .npz file for cits/cims (default: <bam_stem>_pileup.npz)')
    parser.add_argument('--chrom', default=None,
        help='Restrict scan to one chromosome (e.g. chr1). Default: all standard chroms.')
    parser.add_argument('--mapq', type=int, default=255,
        help='Minimum mapping quality (default: 255 = STAR unique-only).')
    parser.add_argument('--nh', type=int, default=1,
        help='Maximum NH tag value — 1 = unique mappers only (default: 1).')
    parser.add_argument('--summary', action='store_true',
        help='Print per-chromosome signal summary (top sites etc).')
    parser.add_argument('--top', type=int, default=10,
        help='Top N sites to show in summary (default: 10). Requires --summary.')
    args = parser.parse_args()

    out_path = args.out or (Path(args.bam).stem + '_pileup.npz')

    print(f"\nClink pileup  |  {Path(args.bam).name}", file=sys.stderr)
    print(f"  mapq>={args.mapq}  NH<={args.nh}  chrom={args.chrom or 'all'}",
          file=sys.stderr)
    print(f"  Output: {out_path}\n", file=sys.stderr)

    pileups = scan_bam(
        args.bam,
        chrom=args.chrom,
        min_mapq=args.mapq,
        max_nh=args.nh,
    )

    chrom_data = {}
    for chrom, strand_pileups in pileups.items():
        chrom_data[chrom] = {}
        for s_label, pileup in strand_pileups.items():
            result = to_arrays(pileup)
            if result is None:
                print(f"  {chrom} [{s_label}]: no signal", file=sys.stderr)
                continue
            positions, coverage, truncations, deletions, subs = result
            chrom_data[chrom][s_label] = (positions, coverage, truncations, deletions, subs)
            if args.summary:
                print_summary(chrom, positions, coverage, truncations, deletions, subs,
                              top_n=args.top)
        # Drop chrom entry if neither strand has signal
        if not chrom_data[chrom]:
            del chrom_data[chrom]

    save_pileup(chrom_data, out_path)
