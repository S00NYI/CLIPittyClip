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
  - CIGAR-tuple walk + numpy bincount/cumsum for coverage and deletions:
      * Collect M/D interval endpoints as Python lists (O(n_reads × avg_cigar_ops))
      * np.bincount builds a difference array in C  (O(n_intervals))
      * np.cumsum converts difference → coverage    (O(span), C-level)
      * No per-base Python loop; eliminates get_aligned_pairs tuple allocation
  - Substitutions: MD tag pre-screened with a regex; get_aligned_pairs only
    called for reads that actually carry mismatches (~5-15% in typical CLIP)
  - Chromosome-level multiprocessing: each worker opens its own BAM handle
    so N chromosomes scan in parallel across N cores (--threads)
  - CIGAR D distinguishes deletion from intron (N) without a reference FASTA
  - Deletions count toward coverage (read spans the position, base is absent)
  - Intron-spanning positions (N) are not counted as coverage or deletions

Usage:
    python pileup.py sample.bam --out sample_pileup.npz
    python pileup.py sample.bam --out sample_pileup.npz --threads 8
    python pileup.py sample.bam --chrom chr1
"""

import re
import sys
import argparse
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pysam


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASES = frozenset('ACGT')

# CIGAR operation codes
_OP_M  = 0   # alignment match (can be match or mismatch)
_OP_I  = 1   # insertion to reference
_OP_D  = 2   # deletion from reference
_OP_N  = 3   # skipped region (intron)
_OP_S  = 4   # soft clip
_OP_H  = 5   # hard clip
_OP_P  = 6   # padding
_OP_EQ = 7   # sequence match
_OP_X  = 8   # sequence mismatch

# Pre-compiled regex: any uppercase letter NOT immediately preceded by '^'
# The '^' prefix marks reference deletions in MD (e.g. ^ACG); bare letters = mismatches
_MD_HAS_MISMATCH = re.compile(r'(?<!\^)[ACGT]')

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

    After the fast CIGAR path, _arrays is set directly so to_arrays()
    can skip re-conversion.  The raw dict fields are only populated when
    the legacy extract_signals() path is used (kept for compatibility).
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
# Truncation site helper
# ---------------------------------------------------------------------------

def truncation_site(read: pysam.AlignedSegment) -> int:
    """
    Truncation site = last nucleotide synthesized by RT (0-based reference coords).

    iCLIP library structure:
      - RT reads 3'→5' on RNA template, stops AT the crosslinked nucleotide
      - The cDNA 3' end = crosslink site; the read begins one position downstream
      - Crosslink nucleotide is ONE POSITION UPSTREAM of the read's 5' end

    In reference coordinates (matching CTK/CITS convention):
      Forward (is_reverse=False): crosslink = reference_start - 1
      Reverse (is_reverse=True):  crosslink = reference_end
        (reference_end is pysam's exclusive end; the 5' start of the read is
         reference_end - 1, so one upstream = reference_end)
    """
    return read.reference_end if read.is_reverse else read.reference_start - 1


# ---------------------------------------------------------------------------
# Substitution extraction (called only for reads with confirmed mismatches)
# ---------------------------------------------------------------------------

def _extract_subs_from_pairs(read: pysam.AlignedSegment,
                              subs_sparse: dict) -> None:
    """
    Extract per-position substitution counts using get_aligned_pairs(with_seq=True).

    Only called when the MD tag regex confirmed at least one mismatch, so the
    fraction of reads reaching this path is typically 5-15% for CLIP data.
    """
    query_seq = read.query_sequence
    if query_seq is None:
        return
    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except Exception:
        return
    for query_pos, ref_pos, ref_base in pairs:
        if query_pos is None or ref_pos is None or ref_base is None:
            continue
        ref_upper = ref_base.upper()
        alt_base  = query_seq[query_pos].upper()
        if ref_upper in BASES and alt_base in BASES and ref_upper != alt_base:
            subs_sparse[ref_pos][(ref_upper, alt_base)] += 1


# ---------------------------------------------------------------------------
# Per-chromosome scan (module-level — required for multiprocessing pickling)
# ---------------------------------------------------------------------------

# Maximum bp per chunk.  Each chunk allocates three int32 arrays of this size:
#   3 arrays × CHUNK_SIZE × 4 bytes = 60 MB per worker at default 5 Mbp.
# With 8 parallel workers the peak pileup footprint is ~480 MB regardless of
# how widely reads are spread across the chromosome.
_CHUNK_SIZE = 5_000_000


def _process_chunk(M_starts: list, M_ends: list,
                   D_starts: list, D_ends: list,
                   T_off: list, size: int,
                   subs_sparse: dict, sub_reads: list) -> tuple:
    """
    Run bincount/cumsum on one chunk's interval lists and return sparse arrays.

    Parameters
    ----------
    M_starts, M_ends : M/=/X block endpoints, already offset to [0, size)
    D_starts, D_ends : D block endpoints, already offset to [0, size)
    T_off            : truncation positions offset to [0, size)
    size             : chunk size in bp
    subs_sparse      : accumulate substitutions here (modified in place)
    sub_reads        : reads with mismatches owned by this chunk

    Returns
    -------
    (local_positions, coverage, truncations, deletions) as numpy arrays,
    or None if the chunk has no signal.
    """
    if not M_starts and not T_off:
        return None

    cov_diff = np.zeros(size + 1, dtype=np.int32)

    # M/=/X coverage
    if M_starts:
        Ms = np.array(M_starts, dtype=np.int64)
        Me = np.array(M_ends,   dtype=np.int64)
        cov_diff += np.bincount(Ms, minlength=size + 1).astype(np.int32)
        cov_diff -= np.bincount(Me, minlength=size + 1).astype(np.int32)

    # Truncation single-base coverage
    if T_off:
        T  = np.array(T_off,         dtype=np.int64)
        T1 = np.clip(T + 1, 0, size).astype(np.int64)
        cov_diff += np.bincount(T,  minlength=size + 1).astype(np.int32)
        cov_diff -= np.bincount(T1, minlength=size + 1).astype(np.int32)
        trunc_arr = np.bincount(T, minlength=size).astype(np.uint32)
    else:
        trunc_arr = np.zeros(size, dtype=np.uint32)

    # D block coverage + deletions
    del_arr = np.zeros(size, dtype=np.uint32)
    if D_starts:
        Ds = np.array(D_starts, dtype=np.int64)
        De = np.array(D_ends,   dtype=np.int64)
        d_diff = (np.bincount(Ds, minlength=size + 1).astype(np.int32) -
                  np.bincount(De, minlength=size + 1).astype(np.int32))
        cov_diff += d_diff
        del_arr   = np.maximum(0, np.cumsum(d_diff[:size])).astype(np.uint32)

    cov_arr = np.maximum(0, np.cumsum(cov_diff[:size])).astype(np.uint32)

    # Substitutions for reads owned by this chunk
    for read in sub_reads:
        _extract_subs_from_pairs(read, subs_sparse)

    # Sparse extraction
    has_signal = (cov_arr > 0) | (trunc_arr > 0) | (del_arr > 0)
    local_idx  = np.where(has_signal)[0]
    if len(local_idx) == 0:
        return None

    return (local_idx,
            cov_arr[local_idx],
            trunc_arr[local_idx],
            del_arr[local_idx])


def _scan_single_chrom(bam_path: str, chrom: str,
                       min_mapq: int, max_nh: int) -> 'ChromPileup':
    """
    Scan one chromosome in fixed-size chunks and return a ChromPileup.

    Algorithm
    ---------
    The chromosome is split into windows of _CHUNK_SIZE bp.  For each window:

    1. bam.fetch(chrom, start, end) returns only reads overlapping that window
       (O(1) BAM index lookup for empty windows — fast skip).

    2. CIGAR intervals (M/D blocks) are clipped to the window boundaries and
       collected as Python lists.  No per-base work in Python.

    3. np.bincount builds a difference array; np.cumsum converts it to
       per-position counts.  Both are C-level — no Python inner loop for
       coverage or deletion accumulation.

    4. Read ownership: a read is "owned" by the chunk containing its
       truncation site.  n_reads, subs, and truncation counts are only
       incremented for owned reads, preventing double-counting across chunks.
       Coverage intervals are always clipped to the current chunk window.

    5. Peak RAM per worker = 3 arrays × _CHUNK_SIZE × 4 bytes = 60 MB at 5 Mbp.
       With 8 parallel workers the total pileup footprint is ~480 MB regardless
       of chromosome length or how widely reads spread.
    """
    n_reads   = 0
    n_skipped = 0

    # Accumulate chunk results across the chromosome
    all_positions:   list = []
    all_coverage:    list = []
    all_truncations: list = []
    all_deletions:   list = []
    subs_sparse: dict = defaultdict(Counter)  # global across chunks

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chrom_len = bam.get_reference_length(chrom)

        for chunk_start in range(0, chrom_len, _CHUNK_SIZE):
            chunk_end = min(chunk_start + _CHUNK_SIZE, chrom_len)
            size      = chunk_end - chunk_start

            # Per-chunk interval lists (offset to [0, size))
            M_starts: list = []
            M_ends:   list = []
            D_starts: list = []
            D_ends:   list = []
            T_off:    list = []   # owned truncation offsets
            sub_reads: list = []  # owned reads with mismatches

            for read in bam.fetch(chrom, chunk_start, chunk_end):
                # ── Filters ──────────────────────────────────────────────
                if (read.is_unmapped or read.is_secondary or
                        read.is_supplementary or read.cigartuples is None):
                    n_skipped += 1
                    continue
                if read.mapping_quality < min_mapq:
                    n_skipped += 1
                    continue
                if min_mapq < 255:
                    try:
                        if read.get_tag('NH') > max_nh:
                            n_skipped += 1
                            continue
                    except KeyError:
                        pass

                # ── Read ownership via truncation site ───────────────────
                t = (read.reference_end if read.is_reverse
                     else read.reference_start - 1)
                owned = chunk_start <= t < chunk_end

                if owned:
                    n_reads += 1
                    t_off = t - chunk_start
                    if 0 <= t_off < size:
                        T_off.append(t_off)
                    # Substitutions counted once per read in owning chunk
                    try:
                        md = read.get_tag('MD')
                        if _MD_HAS_MISMATCH.search(md):
                            sub_reads.append(read)
                    except KeyError:
                        pass

                # ── CIGAR walk — clip intervals to chunk window ───────────
                ref = read.reference_start
                for op, length in read.cigartuples:
                    if op == _OP_M or op == _OP_EQ or op == _OP_X:
                        s = max(ref,          chunk_start) - chunk_start
                        e = min(ref + length, chunk_end)   - chunk_start
                        if s < e:
                            M_starts.append(s)
                            M_ends.append(e)
                        ref += length
                    elif op == _OP_D:
                        s = max(ref,          chunk_start) - chunk_start
                        e = min(ref + length, chunk_end)   - chunk_start
                        if s < e:
                            D_starts.append(s)
                            D_ends.append(e)
                        ref += length
                    elif op == _OP_N:
                        ref += length
                    # I, S, H, P: no reference advance

            # ── Process chunk ─────────────────────────────────────────────
            result = _process_chunk(M_starts, M_ends, D_starts, D_ends,
                                    T_off, size, subs_sparse, sub_reads)
            if result is not None:
                local_idx, cov, trunc, dels = result
                all_positions.append(local_idx + chunk_start)
                all_coverage.append(cov)
                all_truncations.append(trunc)
                all_deletions.append(dels)

    # ── Build ChromPileup ─────────────────────────────────────────────────────
    p = ChromPileup(chrom=chrom)
    p.n_reads   = n_reads
    p.n_skipped = n_skipped

    if not all_positions:
        return p

    # Concatenate chunk arrays
    global_pos  = np.concatenate(all_positions).astype(np.int32)
    coverage    = np.concatenate(all_coverage).astype(np.uint32)
    truncations = np.concatenate(all_truncations).astype(np.uint32)
    deletions   = np.concatenate(all_deletions).astype(np.uint32)

    # Build substitution arrays
    all_sub_types = set()
    for counter in subs_sparse.values():
        all_sub_types.update(counter.keys())

    pos_to_idx = {int(gp): i for i, gp in enumerate(global_pos)}
    subs_out: dict = {}
    for sub_type in sorted(all_sub_types):
        arr = np.zeros(len(global_pos), dtype=np.uint16)
        for pos, counter in subs_sparse.items():
            if sub_type in counter:
                idx = pos_to_idx.get(pos)
                if idx is not None:
                    arr[idx] = min(counter[sub_type], 65535)
        subs_out[sub_type] = arr

    p._arrays = (global_pos, coverage, truncations, deletions, subs_out)  # type: ignore[attr-defined]
    return p


def _worker(args: tuple):
    """
    Multiprocessing worker — must be module-level for pickling.
    Returns (chrom, n_reads, n_skipped, arrays).
    """
    bam_path, chrom, min_mapq, max_nh = args
    p = _scan_single_chrom(bam_path, chrom, min_mapq, max_nh)
    arrays = to_arrays(p)
    return chrom, p.n_reads, p.n_skipped, arrays


# ---------------------------------------------------------------------------
# BAM scan — sequential or parallel
# ---------------------------------------------------------------------------

def scan_bam(
    bam_path: str,
    chrom:    Optional[str] = None,
    min_mapq: int  = 20,
    max_nh:   int  = 1,
    threads:  int  = 1,
    verbose:  bool = True,
) -> Dict[str, ChromPileup]:
    """
    One-pass BAM scan. Returns {chrom: ChromPileup}.

    When threads > 1, chromosomes are scanned in parallel using
    multiprocessing.Pool. Each worker opens its own BAM file handle —
    pysam is safe for concurrent reads on the same file.

    Args:
        bam_path : sorted, indexed BAM
        chrom    : scan only this chromosome (None = all standard chroms)
        min_mapq : minimum MAPQ filter
        max_nh   : max NH tag; 1 = unique mappers only
        threads  : parallel workers (one per chromosome); 1 = sequential
        verbose  : print per-chromosome summary to stderr
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chroms_to_scan = (
            [chrom] if chrom
            else [sq['SN'] for sq in bam.header['SQ']
                  if is_standard_chrom(sq['SN'])]
        )

    pileups: Dict[str, ChromPileup] = {}

    def _log(p: ChromPileup, arrays) -> None:
        if not verbose or arrays is None:
            return
        positions, coverage, truncations, deletions, subs = arrays
        print(
            f"  {p.chrom}: {p.n_reads:>10,} reads | "
            f"{len(positions):>10,} covered positions | "
            f"{int(truncations.sum()):>8,} truncation events | "
            f"{int(deletions.sum()):>8,} deletion events | "
            f"{p.n_skipped:>8,} skipped",
            file=sys.stderr,
        )

    if threads > 1:
        # ── Parallel path ────────────────────────────────────────────────────
        from multiprocessing import Pool
        n_workers = min(threads, len(chroms_to_scan))
        if verbose:
            print(f"  Scanning {len(chroms_to_scan)} chromosomes "
                  f"with {n_workers} parallel workers ...", file=sys.stderr)

        jobs = [(bam_path, c, min_mapq, max_nh) for c in chroms_to_scan]
        with Pool(n_workers) as pool:
            for c, n_reads, n_skipped, arrays in pool.map(_worker, jobs):
                p = ChromPileup(chrom=c)
                p.n_reads   = n_reads
                p.n_skipped = n_skipped
                p._arrays   = arrays  # type: ignore[attr-defined]
                pileups[c]  = p
                _log(p, arrays)

    else:
        # ── Sequential path ──────────────────────────────────────────────────
        for c in chroms_to_scan:
            p = _scan_single_chrom(bam_path, c, min_mapq, max_nh)
            arrays = to_arrays(p)
            p._arrays = arrays  # type: ignore[attr-defined]
            pileups[c] = p
            _log(p, arrays)

    return pileups


# ---------------------------------------------------------------------------
# Convert to numpy arrays for cits.py / cims.py
# ---------------------------------------------------------------------------

def to_arrays(pileup: ChromPileup):
    """
    Return pileup arrays for one chromosome.

    Fast path: if _arrays was pre-computed by _scan_single_chrom (the CIGAR
    path), return it immediately.  Slow path: convert legacy defaultdict fields
    (only reached if extract_signals() was used directly).

    Returns:
        (positions, coverage, truncations, deletions, subs)
          positions    : int32  [N]  — 0-based reference positions with signal
          coverage     : uint32 [N]
          truncations  : uint32 [N]
          deletions    : uint32 [N]
          subs         : dict {(ref, alt): uint16 [N]}
        or None if the chromosome has no signal.
    """
    # Fast path
    if hasattr(pileup, '_arrays'):
        return pileup._arrays  # type: ignore[attr-defined]

    # Legacy path (defaultdict → arrays)
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

    for pos, v in pileup.coverage.items():
        coverage[idx[pos]] = v
    for pos, v in pileup.truncations.items():
        truncations[idx[pos]] = v
    for pos, v in pileup.deletions.items():
        deletions[idx[pos]] = v

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

    Key naming convention (double-underscore separator):
        {chrom}__positions
        {chrom}__coverage
        {chrom}__truncations
        {chrom}__deletions
        {chrom}__sub__{ref}{alt}   e.g. chr1__sub__TC
    """
    arrays = {}
    for chrom, (positions, coverage, truncations, deletions, subs) in chrom_data.items():
        arrays[f'{chrom}__positions']   = positions
        arrays[f'{chrom}__coverage']    = coverage
        arrays[f'{chrom}__truncations'] = truncations
        arrays[f'{chrom}__deletions']   = deletions
        for (ref, alt), arr in subs.items():
            arrays[f'{chrom}__sub__{ref}{alt}'] = arr

    np.savez_compressed(path, **arrays)
    print(f"  Pileup saved → {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)


def load_pileup(path: str) -> dict:
    """
    Load pileup arrays from a .npz file written by save_pileup().
    Returns: {chrom: (positions, coverage, truncations, deletions, subs)}
    """
    data = np.load(path)

    chroms = set()
    for key in data.files:
        chroms.add(key.split('__')[0])

    chrom_data = {}
    for chrom in sorted(chroms):
        positions   = data[f'{chrom}__positions']
        coverage    = data[f'{chrom}__coverage']
        truncations = data[f'{chrom}__truncations']
        deletions   = data[f'{chrom}__deletions']

        subs = {}
        prefix = f'{chrom}__sub__'
        for key in data.files:
            if key.startswith(prefix):
                code = key[len(prefix):]
                ref, alt = code[0], code[1]
                subs[(ref, alt)] = data[key]

        chrom_data[chrom] = (positions, coverage, truncations, deletions, subs)

    print(f"  Pileup loaded ← {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)
    return chrom_data


# ---------------------------------------------------------------------------
# Quick summary report
# ---------------------------------------------------------------------------

def print_summary(chrom: str, positions, coverage, truncations, deletions, subs,
                  top_n: int = 10):
    print(f"\n{'='*60}")
    print(f"  {chrom}  —  {len(positions):,} positions with signal")
    print(f"{'='*60}")
    print(f"  Max coverage:       {coverage.max():>10,}")
    print(f"  Total truncations:  {truncations.sum():>10,}")
    print(f"  Total deletions:    {deletions.sum():>10,}")
    print(f"  Substitution types: {[f'{r}>{a}' for r,a in sorted(subs.keys())]}")
    for (ref, alt), arr in sorted(subs.items()):
        print(f"    {ref}>{alt}: {arr.sum():,} total events")

    if truncations.sum() > 0:
        print(f"\n  Top {top_n} truncation sites:")
        top = np.argsort(truncations)[-top_n:][::-1]
        for i in top:
            if truncations[i] > 0:
                cov = coverage[i] if coverage[i] > 0 else 1
                print(f"    pos={positions[i]:>12,}  "
                      f"trunc={truncations[i]:>6,}  "
                      f"cov={cov:>8,}  "
                      f"frac={truncations[i]/cov:.3f}")

    if deletions.sum() > 0:
        print(f"\n  Top {top_n} deletion sites:")
        top = np.argsort(deletions)[-top_n:][::-1]
        for i in top:
            if deletions[i] > 0:
                cov = coverage[i] if coverage[i] > 0 else 1
                print(f"    pos={positions[i]:>12,}  "
                      f"del={deletions[i]:>6,}  "
                      f"cov={cov:>8,}  "
                      f"frac={deletions[i]/cov:.3f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink pileup: scan BAM → per-position signal arrays (.npz)')
    parser.add_argument('bam',
        help='Sorted, indexed BAM file (output of clink collapse)')
    parser.add_argument('--out', default=None,
        help='Output .npz file (default: <bam_stem>_pileup.npz)')
    parser.add_argument('--chrom', default=None,
        help='Restrict to one chromosome (default: all standard chroms)')
    parser.add_argument('--mapq', type=int, default=255,
        help='Minimum MAPQ (default: 255 = STAR unique-only)')
    parser.add_argument('--nh', type=int, default=1,
        help='Maximum NH tag — 1 = unique mappers only (default: 1)')
    parser.add_argument('--threads', type=int, default=1,
        help='Parallel workers — one per chromosome (default: 1 = sequential)')
    parser.add_argument('--summary', action='store_true',
        help='Print per-chromosome signal summary')
    parser.add_argument('--top', type=int, default=10,
        help='Top N sites in summary (default: 10)')
    args = parser.parse_args()

    out_path = args.out or (Path(args.bam).stem + '_pileup.npz')

    print(f"\nClink pileup  |  {Path(args.bam).name}", file=sys.stderr)
    print(f"  mapq>={args.mapq}  NH<={args.nh}  "
          f"threads={args.threads}  chrom={args.chrom or 'all'}",
          file=sys.stderr)
    print(f"  Output: {out_path}\n", file=sys.stderr)

    pileups = scan_bam(
        args.bam,
        chrom    = args.chrom,
        min_mapq = args.mapq,
        max_nh   = args.nh,
        threads  = args.threads,
    )

    chrom_data = {}
    for chrom, pileup in pileups.items():
        result = to_arrays(pileup)
        if result is None:
            print(f"  {chrom}: no signal", file=sys.stderr)
            continue
        positions, coverage, truncations, deletions, subs = result
        chrom_data[chrom] = (positions, coverage, truncations, deletions, subs)
        if args.summary:
            print_summary(chrom, positions, coverage, truncations, deletions, subs,
                          top_n=args.top)

    save_pileup(chrom_data, out_path)
