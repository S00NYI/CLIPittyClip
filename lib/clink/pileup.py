#!/usr/bin/env python3
"""
Clink pileup.py — Single-pass BAM signal extractor (strand-aware)

Scans a sorted, indexed BAM and accumulates four per-position signals
separately for forward (+) and reverse (-) strand reads:

  coverage     : reads spanning each position
  truncations  : 5' read ends (RT stop sites — CITS signal)
  deletions    : 1-bp CIGAR D operations (crosslink-induced deletions — CIMS signal)
  substitutions: mismatches by type, e.g. T>C (PAR-CLIP), any sub (iCLIP)

Key design decisions:
  - Strand separation: all signal arrays accumulated independently per strand
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

NPZ format (stranded):
  {chrom}__fwd__positions     {chrom}__rev__positions
  {chrom}__fwd__coverage      {chrom}__rev__coverage
  {chrom}__fwd__truncations   {chrom}__rev__truncations
  {chrom}__fwd__deletions     {chrom}__rev__deletions
  {chrom}__fwd__sub__{XY}     {chrom}__rev__sub__{XY}

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

# Complement map for strand-aware substitution recording.
# Reverse-strand reads are stored in reverse-complement orientation in the BAM,
# so pysam reports substitutions against the + strand reference.  To convert
# to the biological (gene-strand) context we complement both ref and alt:
#   e.g. A→G on a − strand read = T→C in the RNA (4SU crosslink signal)
_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

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

# Standard chromosomes — handles both UCSC (chr1-22, chrX, chrY, chrM)
# and Ensembl (1-22, X, Y, MT) naming conventions.
def is_standard_chrom(name: str) -> bool:
    # UCSC style
    if name in ('chrX', 'chrY', 'chrM', 'chrMT'):
        return True
    if name.startswith('chr') and name[3:].isdigit():
        return 1 <= int(name[3:]) <= 22
    # Ensembl style
    if name in ('X', 'Y', 'MT'):
        return True
    if name.isdigit():
        return 1 <= int(name) <= 22
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

    _arrays format (stranded):
        {'fwd': (positions, coverage, truncations, deletions, subs) or None,
         'rev': (positions, coverage, truncations, deletions, subs) or None}
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

    In reference coordinates (matching CTK/CITS convention):
      Forward (is_reverse=False): crosslink = reference_start - 1
      Reverse (is_reverse=True):  crosslink = reference_end
        (reference_end is pysam's exclusive end; the 5' start of the reverse read
         is at reference_end - 1 in reference coords; one upstream = reference_end)
    """
    return read.reference_end if read.is_reverse else read.reference_start - 1


# ---------------------------------------------------------------------------
# Substitution extraction (called only for reads with confirmed mismatches)
# ---------------------------------------------------------------------------

def _extract_subs_from_pairs(read: pysam.AlignedSegment,
                              subs_sparse: dict,
                              is_reverse: bool = False) -> None:
    """
    Extract per-position substitution counts using get_aligned_pairs(with_seq=True).

    Only called when the MD tag regex confirmed at least one mismatch, so the
    fraction of reads reaching this path is typically 5-15% for CLIP data.

    Strand convention
    -----------------
    pysam always reports ref_base in + strand orientation, and for reverse-strand
    reads the stored query sequence is the reverse complement of the RNA.  This
    means the same biological event (e.g. a 4SU T→C crosslink) appears as:
      - T→C on a + strand read  (ref=T, read=C)
      - A→G on a − strand read  (ref=A, read=G — reverse complement of T→C)

    When is_reverse=True we complement both ref and alt before recording, so
    all substitution types are stored in biological (gene-strand) context:
      A→G on − strand → _COMP('A')=T, _COMP('G')=C → T→C  ✓
    This ensures T→C arrays accumulate signal from both strand orientations,
    making background rates and site calls biologically coherent.
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
            if is_reverse:
                # Convert to biological (gene-strand) context
                ref_upper = _COMP[ref_upper]
                alt_base  = _COMP[alt_base]
            subs_sparse[ref_pos][(ref_upper, alt_base)] += 1


# ---------------------------------------------------------------------------
# Per-chromosome scan (module-level — required for multiprocessing pickling)
# ---------------------------------------------------------------------------

# Maximum bp per chunk.  Each chunk allocates six int32 arrays of this size
# (three per strand: cov_diff, trunc, del):
#   6 arrays × CHUNK_SIZE × 4 bytes = 120 MB per worker at default 5 Mbp.
# With 8 parallel workers the peak pileup footprint is ~960 MB.
_CHUNK_SIZE = 5_000_000


def _process_chunk(M_starts: list, M_ends: list,
                   D_starts: list, D_ends: list,
                   T_off: list, T_clean: list, size: int,
                   subs_sparse: dict, sub_reads: list,
                   is_reverse: bool = False) -> tuple:
    """
    Run bincount/cumsum on one chunk's interval lists and return sparse arrays.

    Parameters
    ----------
    M_starts, M_ends : M/=/X block endpoints, already offset to [0, size)
    D_starts, D_ends : D block endpoints, already offset to [0, size)
    T_off            : truncation positions offset to [0, size)
    T_clean          : truncation positions for reads with no deletion in CIGAR
    size             : chunk size in bp
    subs_sparse      : accumulate substitutions here (modified in place)
    sub_reads        : reads with mismatches owned by this chunk

    Returns
    -------
    (local_positions, coverage, truncations, deletions, clean_truncations)
    as numpy arrays, or None if the chunk has no signal.
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

    # Clean truncations: reads without any deletion in their CIGAR
    if T_clean:
        Tc = np.array(T_clean, dtype=np.int64)
        clean_trunc_arr = np.bincount(Tc, minlength=size).astype(np.uint32)
    else:
        clean_trunc_arr = np.zeros(size, dtype=np.uint32)

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
        _extract_subs_from_pairs(read, subs_sparse, is_reverse=is_reverse)

    # Sparse extraction
    has_signal = (cov_arr > 0) | (trunc_arr > 0) | (del_arr > 0)
    local_idx  = np.where(has_signal)[0]
    if len(local_idx) == 0:
        return None

    return (local_idx,
            cov_arr[local_idx],
            trunc_arr[local_idx],
            del_arr[local_idx],
            clean_trunc_arr[local_idx])


def _build_strand_arrays(all_positions: list, all_coverage: list,
                          all_truncations: list, all_deletions: list,
                          all_clean_truncations: list,
                          subs_sparse: dict):
    """
    Concatenate chunk accumulation lists into final numpy arrays for one strand.
    Returns (positions, coverage, truncations, deletions, clean_truncations, subs_out)
    or None.
    """
    if not all_positions:
        return None

    global_pos        = np.concatenate(all_positions).astype(np.int32)
    coverage          = np.concatenate(all_coverage).astype(np.uint32)
    truncations       = np.concatenate(all_truncations).astype(np.uint32)
    deletions         = np.concatenate(all_deletions).astype(np.uint32)
    clean_truncations = np.concatenate(all_clean_truncations).astype(np.uint32)

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

    return (global_pos, coverage, truncations, deletions, clean_truncations, subs_out)


def _scan_single_chrom(bam_path: str, chrom: str,
                       min_mapq: int, max_nh: int) -> 'ChromPileup':
    """
    Scan one chromosome in fixed-size chunks and return a ChromPileup.

    Reads are split by strand at the CIGAR-walk level so that fwd (+) and
    rev (-) signals accumulate into independent arrays throughout.

    Algorithm
    ---------
    The chromosome is split into windows of _CHUNK_SIZE bp.  For each window:

    1. bam.fetch(chrom, start, end) returns only reads overlapping that window.

    2. Reads are split by strand; CIGAR intervals (M/D blocks) are clipped to
       the window boundaries and collected into fwd and rev Python lists.

    3. np.bincount + np.cumsum (C-level) convert interval lists to per-position
       coverage/deletion arrays for each strand.

    4. Read ownership: a read is "owned" by the chunk containing its truncation
       site.  n_reads, subs, and truncation counts are only incremented for owned
       reads.  Coverage intervals are accumulated in every overlapping chunk.

    5. Peak RAM per worker ≈ 6 arrays × _CHUNK_SIZE × 4 bytes = 120 MB at 5 Mbp.
    """
    n_reads   = 0
    n_skipped = 0

    # Accumulate chunk results separately per strand
    all_positions_fwd:         list = []
    all_coverage_fwd:          list = []
    all_truncations_fwd:       list = []
    all_deletions_fwd:         list = []
    all_clean_truncations_fwd: list = []
    subs_sparse_fwd: dict = defaultdict(Counter)

    all_positions_rev:         list = []
    all_coverage_rev:          list = []
    all_truncations_rev:       list = []
    all_deletions_rev:         list = []
    all_clean_truncations_rev: list = []
    subs_sparse_rev: dict = defaultdict(Counter)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        chrom_len = bam.get_reference_length(chrom)

        for chunk_start in range(0, chrom_len, _CHUNK_SIZE):
            chunk_end = min(chunk_start + _CHUNK_SIZE, chrom_len)
            size      = chunk_end - chunk_start

            # Per-chunk interval lists, split by strand
            M_starts_fwd: list = []; M_ends_fwd: list = []
            M_starts_rev: list = []; M_ends_rev: list = []
            D_starts_fwd: list = []; D_ends_fwd: list = []
            D_starts_rev: list = []; D_ends_rev: list = []
            T_off_fwd:       list = []; T_off_rev:       list = []
            T_off_fwd_clean: list = []; T_off_rev_clean: list = []
            sub_reads_fwd: list = []; sub_reads_rev: list = []

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

                is_rev = read.is_reverse

                # ── Pre-check: does this read have any deletion in CIGAR? ─
                # Used to populate clean_truncations (CITS signal excluding
                # reads that also carry a deletion event, matching CTK's
                # removeRow step before CITS calling).
                has_del = any(op == _OP_D for op, _ in read.cigartuples)

                # ── Read ownership via truncation site ───────────────────
                t = (read.reference_end if is_rev
                     else read.reference_start - 1)
                owned = chunk_start <= t < chunk_end

                if owned:
                    n_reads += 1
                    t_off = t - chunk_start
                    if 0 <= t_off < size:
                        if is_rev:
                            T_off_rev.append(t_off)
                            if not has_del:
                                T_off_rev_clean.append(t_off)
                        else:
                            T_off_fwd.append(t_off)
                            if not has_del:
                                T_off_fwd_clean.append(t_off)
                    # Substitutions counted once per read in owning chunk
                    try:
                        md = read.get_tag('MD')
                        if _MD_HAS_MISMATCH.search(md):
                            if is_rev:
                                sub_reads_rev.append(read)
                            else:
                                sub_reads_fwd.append(read)
                    except KeyError:
                        pass

                # ── CIGAR walk — clip intervals to chunk window ───────────
                ref = read.reference_start
                for op, length in read.cigartuples:
                    if op == _OP_M or op == _OP_EQ or op == _OP_X:
                        s = max(ref,          chunk_start) - chunk_start
                        e = min(ref + length, chunk_end)   - chunk_start
                        if s < e:
                            if is_rev:
                                M_starts_rev.append(s)
                                M_ends_rev.append(e)
                            else:
                                M_starts_fwd.append(s)
                                M_ends_fwd.append(e)
                        ref += length
                    elif op == _OP_D:
                        s = max(ref,          chunk_start) - chunk_start
                        e = min(ref + length, chunk_end)   - chunk_start
                        if s < e:
                            if is_rev:
                                D_starts_rev.append(s)
                                D_ends_rev.append(e)
                            else:
                                D_starts_fwd.append(s)
                                D_ends_fwd.append(e)
                        ref += length
                    elif op == _OP_N:
                        ref += length
                    # I, S, H, P: no reference advance

            # ── Process chunk — both strands independently ────────────────
            result_fwd = _process_chunk(
                M_starts_fwd, M_ends_fwd, D_starts_fwd, D_ends_fwd,
                T_off_fwd, T_off_fwd_clean, size, subs_sparse_fwd, sub_reads_fwd)
            if result_fwd is not None:
                local_idx, cov, trunc, dels, clean_trunc = result_fwd
                all_positions_fwd.append(local_idx + chunk_start)
                all_coverage_fwd.append(cov)
                all_truncations_fwd.append(trunc)
                all_deletions_fwd.append(dels)
                all_clean_truncations_fwd.append(clean_trunc)

            result_rev = _process_chunk(
                M_starts_rev, M_ends_rev, D_starts_rev, D_ends_rev,
                T_off_rev, T_off_rev_clean, size, subs_sparse_rev, sub_reads_rev,
                is_reverse=True)
            if result_rev is not None:
                local_idx, cov, trunc, dels, clean_trunc = result_rev
                all_positions_rev.append(local_idx + chunk_start)
                all_coverage_rev.append(cov)
                all_truncations_rev.append(trunc)
                all_deletions_rev.append(dels)
                all_clean_truncations_rev.append(clean_trunc)

    # ── Build ChromPileup ─────────────────────────────────────────────────────
    p = ChromPileup(chrom=chrom)
    p.n_reads   = n_reads
    p.n_skipped = n_skipped

    fwd_arrays = _build_strand_arrays(
        all_positions_fwd, all_coverage_fwd, all_truncations_fwd,
        all_deletions_fwd, all_clean_truncations_fwd, subs_sparse_fwd)
    rev_arrays = _build_strand_arrays(
        all_positions_rev, all_coverage_rev, all_truncations_rev,
        all_deletions_rev, all_clean_truncations_rev, subs_sparse_rev)

    p._arrays = {'fwd': fwd_arrays, 'rev': rev_arrays}  # type: ignore[attr-defined]
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
        n_pos = n_trunc = n_del = 0
        for strand_arr in (arrays.get('fwd'), arrays.get('rev')):
            if strand_arr is None:
                continue
            positions, coverage, truncations, deletions, clean_truncations, subs = strand_arr
            n_pos   += len(positions)
            n_trunc += int(truncations.sum())
            n_del   += int(deletions.sum())
        print(
            f"  {p.chrom}: {p.n_reads:>10,} reads | "
            f"{n_pos:>10,} covered positions (fwd+rev) | "
            f"{n_trunc:>8,} truncation events | "
            f"{n_del:>8,} deletion events | "
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
    path), return it immediately.

    Returns:
        dict {'fwd': (positions, coverage, truncations, deletions, subs) or None,
              'rev': (positions, coverage, truncations, deletions, subs) or None}
        or None if both strands have no signal.

        Each strand tuple:
          positions    : int32  [N]  — 0-based reference positions with signal
          coverage     : uint32 [N]
          truncations  : uint32 [N]
          deletions    : uint32 [N]
          subs         : dict {(ref, alt): uint16 [N]}
    """
    # Fast path
    if hasattr(pileup, '_arrays'):
        return pileup._arrays  # type: ignore[attr-defined]

    # Legacy path (defaultdict → arrays) — wrap everything as fwd (unstranded)
    all_pos = (set(pileup.coverage) | set(pileup.truncations) |
               set(pileup.deletions) | set(pileup.subs))
    if not all_pos:
        return None

    positions = np.array(sorted(all_pos), dtype=np.int32)
    n = len(positions)
    idx = {p: i for i, p in enumerate(positions)}

    coverage          = np.zeros(n, dtype=np.uint32)
    truncations       = np.zeros(n, dtype=np.uint32)
    deletions         = np.zeros(n, dtype=np.uint32)
    clean_truncations = np.zeros(n, dtype=np.uint32)  # legacy: assume all clean

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

    legacy_arrays = (positions, coverage, truncations, deletions, clean_truncations, subs)
    return {'fwd': legacy_arrays, 'rev': None}


# ---------------------------------------------------------------------------
# Save / load pileup arrays (intermediate .npz cache)
# ---------------------------------------------------------------------------

def save_pileup(chrom_data: dict, path: str) -> None:
    """
    Save stranded pileup arrays for all chromosomes to a compressed .npz file.

    chrom_data: {chrom: {'fwd': arrays_or_None, 'rev': arrays_or_None}}

    Key naming convention (triple double-underscore segments):
        {chrom}__fwd__positions       {chrom}__rev__positions
        {chrom}__fwd__coverage        {chrom}__rev__coverage
        {chrom}__fwd__truncations     {chrom}__rev__truncations
        {chrom}__fwd__deletions       {chrom}__rev__deletions
        {chrom}__fwd__sub__{XY}       {chrom}__rev__sub__{XY}
    """
    arrays = {}
    for chrom, strands in chrom_data.items():
        for strand_key in ('fwd', 'rev'):
            strand_arr = strands.get(strand_key)
            if strand_arr is None:
                continue
            positions, coverage, truncations, deletions, clean_truncations, subs = strand_arr
            pfx = f'{chrom}__{strand_key}'
            arrays[f'{pfx}__positions']         = positions
            arrays[f'{pfx}__coverage']          = coverage
            arrays[f'{pfx}__truncations']       = truncations
            arrays[f'{pfx}__deletions']         = deletions
            arrays[f'{pfx}__clean_truncations'] = clean_truncations
            for (ref, alt), arr in subs.items():
                arrays[f'{pfx}__sub__{ref}{alt}'] = arr

    np.savez_compressed(path, **arrays)
    print(f"  Pileup saved → {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)


def load_pileup(path: str) -> dict:
    """
    Load pileup arrays from a .npz file written by save_pileup().

    Returns: {chrom: {'fwd': arrays_or_None, 'rev': arrays_or_None}}

    Backward compatible: if the file uses the old unstranded format
    ({chrom}__positions etc.), wraps all data as fwd-only (rev=None) with
    a warning.
    """
    data = np.load(path)

    # Detect format: stranded files have '__fwd__' or '__rev__' in keys
    is_stranded = any(('__fwd__' in k or '__rev__' in k) for k in data.files)

    if not is_stranded:
        # --- Backward-compatible: old unstranded format → wrap as fwd ---
        print("  WARNING: pileup.npz is unstranded (old format). "
              "Re-run pileup.py to get strand-aware results.", file=sys.stderr)
        chroms = set(k.split('__')[0] for k in data.files)
        chrom_data = {}
        for chrom in sorted(chroms):
            positions   = data[f'{chrom}__positions']
            coverage    = data[f'{chrom}__coverage']
            truncations = data[f'{chrom}__truncations']
            deletions   = data[f'{chrom}__deletions']
            # Backward compat: old npz has no clean_truncations → assume all clean
            ct_key = f'{chrom}__clean_truncations'
            clean_truncations = data[ct_key] if ct_key in data.files else truncations.copy()
            subs = {}
            sub_prefix = f'{chrom}__sub__'
            for key in data.files:
                if key.startswith(sub_prefix):
                    code = key[len(sub_prefix):]
                    ref, alt = code[0], code[1]
                    subs[(ref, alt)] = data[key]
            chrom_data[chrom] = {
                'fwd': (positions, coverage, truncations, deletions, clean_truncations, subs),
                'rev': None,
            }
        print(f"  Pileup loaded (unstranded→fwd) ← {path}  "
              f"({len(chrom_data)} chromosomes)", file=sys.stderr)
        return chrom_data

    # --- Stranded format ---
    chroms = set()
    for key in data.files:
        parts = key.split('__')
        if len(parts) >= 3 and parts[1] in ('fwd', 'rev'):
            chroms.add(parts[0])

    chrom_data = {}
    for chrom in sorted(chroms):
        chrom_data[chrom] = {}
        for strand_key in ('fwd', 'rev'):
            pos_key = f'{chrom}__{strand_key}__positions'
            if pos_key not in data.files:
                chrom_data[chrom][strand_key] = None
                continue
            pfx = f'{chrom}__{strand_key}'
            positions   = data[f'{pfx}__positions']
            coverage    = data[f'{pfx}__coverage']
            truncations = data[f'{pfx}__truncations']
            deletions   = data[f'{pfx}__deletions']
            # Backward compat: old npz has no clean_truncations → assume all clean
            ct_key = f'{pfx}__clean_truncations'
            clean_truncations = data[ct_key] if ct_key in data.files else truncations.copy()
            subs = {}
            sub_prefix = f'{pfx}__sub__'
            for key in data.files:
                if key.startswith(sub_prefix):
                    code = key[len(sub_prefix):]
                    ref, alt = code[0], code[1]
                    subs[(ref, alt)] = data[key]
            chrom_data[chrom][strand_key] = (
                positions, coverage, truncations, deletions, clean_truncations, subs)

    print(f"  Pileup loaded ← {path}  ({len(chrom_data)} chromosomes)", file=sys.stderr)
    return chrom_data


# ---------------------------------------------------------------------------
# Quick summary report
# ---------------------------------------------------------------------------

def print_summary(chrom: str, strand: str, positions, coverage, truncations,
                  deletions, clean_truncations, subs, top_n: int = 10):
    print(f"\n{'='*60}")
    print(f"  {chrom} [{strand}]  —  {len(positions):,} positions with signal")
    print(f"{'='*60}")
    print(f"  Max coverage:            {coverage.max():>10,}")
    print(f"  Total truncations:       {truncations.sum():>10,}")
    print(f"  Clean truncations (no D):{clean_truncations.sum():>10,}")
    print(f"  Total deletions:         {deletions.sum():>10,}")
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
        description='Clink pileup: scan BAM → per-position strand-aware signal arrays (.npz)')
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
        chrom_data[chrom] = result
        if args.summary:
            for strand_key in ('fwd', 'rev'):
                strand_arr = result.get(strand_key)
                if strand_arr is None:
                    continue
                positions, coverage, truncations, deletions, clean_truncations, subs = strand_arr
                strand_char = '+' if strand_key == 'fwd' else '-'
                print_summary(chrom, strand_char, positions, coverage,
                              truncations, deletions, clean_truncations, subs, top_n=args.top)

    save_pileup(chrom_data, out_path)
