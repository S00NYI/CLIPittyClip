#!/usr/bin/env python3
"""
Clink collapse.py — PCR duplicate removal via UMI-aware deduplication

Wraps umi_tools dedup with sensible iCLIP defaults. Reads UMIs directly
from the BAM read names (format: {readID}#{count}#{UMI} or {readID}#{UMI}).
Outputs a deduplicated BAM ready for clink pileup.

UMI extraction:
    Read names from CLIPittyClip are formatted as:
        SRR26536514.3066536#1#GGAACCT   (after fastq2collapse: #count#UMI)
        SRR26536514.3066536#GGAACCT     (without fastq2collapse: #UMI)
    In both cases the UMI is the last '#'-delimited token.

Deduplication method:
    --method directional (umi_tools default):
        Builds a UMI graph per genomic position. Two UMIs are considered
        related if they differ by ≤1 edit AND one has ≥2× the count of the
        other. Collapses related UMIs into the most abundant representative.
        More conservative than 'unique' (keeps more reads at saturated sites).

No-UMI mode:
    When --umi-len 0 is passed, falls back to position-only deduplication
    (umi_tools dedup --method unique --ignore-umi).

Usage:
    python collapse.py \\
        --bam   /path/to/sorted.bam \\
        --out   /path/to/dedup.bam \\
        --umi-len 7                          # UMI length (0 = no UMI)
        --threads 4

    # Position-only (no UMI):
    python collapse.py --bam input.bam --out dedup.bam --umi-len 0

    # Rescue multi-mapped reads via EM assignment (requires pysam):
    python collapse.py --bam input.bam --out dedup.bam --multi-map
"""

import sys
import os
import argparse
import subprocess
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    import pysam
    import numpy as np
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


# ---------------------------------------------------------------------------
# umi_tools environment
# ---------------------------------------------------------------------------

def find_umi_tools() -> str:
    """
    Return path to umi_tools binary.
    Checks PATH first, then the dedicated umi_tools conda env.
    """
    # Check PATH
    path = shutil.which('umi_tools')
    if path:
        return path

    # Try dedicated conda env (installed separately due to Python 3.13 incompatibility)
    conda_candidates = [
        os.path.expanduser('~/anaconda3/envs/umi_tools/bin/umi_tools'),
        os.path.expanduser('~/miniconda3/envs/umi_tools/bin/umi_tools'),
        '/opt/anaconda3/envs/umi_tools/bin/umi_tools',
        '/opt/miniconda3/envs/umi_tools/bin/umi_tools',
    ]
    for p in conda_candidates:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p

    # Try conda run as fallback
    result = subprocess.run(
        ['conda', 'run', '-n', 'umi_tools', 'which', 'umi_tools'],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        found = result.stdout.strip()
        if found:
            return found

    return None


# ---------------------------------------------------------------------------
# BAM utilities
# ---------------------------------------------------------------------------

def check_bam_index(bam_path: str) -> None:
    """Ensure BAM is indexed; create index if missing."""
    bai = bam_path + '.bai'
    if not os.path.exists(bai):
        print(f"  Indexing {Path(bam_path).name} ...", file=sys.stderr)
        subprocess.run(['samtools', 'index', bam_path], check=True)


def detect_umi_in_bam(bam_path: str, n_reads: int = 200) -> tuple[bool, int]:
    """
    Peek at the first n_reads read names to detect UMI presence and length.

    Returns:
        (has_umi: bool, umi_len: int)
    """
    result = subprocess.run(
        ['samtools', 'view', bam_path],
        capture_output=True, text=True
    )
    lines = result.stdout.strip().split('\n')[:n_reads]

    umi_lengths = []
    for line in lines:
        if not line:
            continue
        qname = line.split('\t')[0]
        parts = qname.split('#')
        if len(parts) >= 2:
            umi_candidate = parts[-1]
            if all(c in 'ACGTN' for c in umi_candidate) and len(umi_candidate) >= 3:
                umi_lengths.append(len(umi_candidate))

    if not umi_lengths:
        return False, 0

    # Use most common length
    from collections import Counter
    most_common_len = Counter(umi_lengths).most_common(1)[0][0]
    return True, most_common_len


def bam_read_count(bam_path: str) -> int:
    """Count mapped reads in BAM."""
    result = subprocess.run(
        ['samtools', 'view', '-c', '-F', '4', bam_path],
        capture_output=True, text=True
    )
    try:
        return int(result.stdout.strip())
    except ValueError:
        return 0


def filter_multimappers(bam_path: str, out_path: str,
                        max_nh: int = 1, threads: int = 4) -> int:
    """
    Filter a BAM to reads with NH tag <= max_nh.

    Reads with no NH tag are kept (treated as uniquely mapped — some aligners
    omit the tag for unique mappers). Unmapped reads (FLAG 4) are excluded.

    Uses: samtools view | awk NH filter | samtools view -bS

    Args:
        bam_path : input BAM
        out_path : filtered output BAM (sorted, indexed)
        max_nh   : maximum NH value to keep (default: 1 = unique mappers only)
        threads  : samtools threads

    Returns:
        Number of reads passing filter.
    """
    # Build awk pattern matching NH:i:0 through NH:i:max_nh
    nh_pattern = '|'.join(str(i) for i in range(0, max_nh + 1))
    # Keep: header lines  OR  reads with NH <= max_nh  OR  reads with no NH tag
    awk_expr = f'/^@/ || /\\tNH:i:({nh_pattern})(\\t|$)/ || !/\\tNH:i:/'

    cmd = (
        f"samtools view -h -F 4 '{bam_path}' | "
        f"awk '{awk_expr}' | "
        f"samtools view -bS -@ {threads} -o '{out_path}'"
    )
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(['samtools', 'index', out_path], check=True)
    return bam_read_count(out_path)


# ---------------------------------------------------------------------------
# umi_tools dedup helper
# ---------------------------------------------------------------------------

def _run_umi_tools_dedup(umi_tools: str,
                         in_bam:    str,
                         out_bam:   str,
                         log_path:  str,
                         umi_len:   int) -> None:
    """
    Run umi_tools dedup and exit on failure.

    Internal helper used by both standard mode and multi-map mode.
    The caller is responsible for sorting/indexing inputs; this function
    indexes the output BAM.

    Args:
        umi_tools : path to umi_tools binary
        in_bam    : sorted, indexed input BAM
        out_bam   : output deduplicated BAM
        log_path  : umi_tools log file path
        umi_len   : UMI length (0 → --ignore-umi / position-only)
    """
    cmd = [
        umi_tools, 'dedup',
        '-I', in_bam,
        '-S', out_bam,
        '-L', log_path,
        '--method', 'directional',
        '--spliced-is-unique',
    ]
    if umi_len > 0:
        cmd += [
            '--extract-umi-method', 'read_id',
            '--umi-separator', '#',
        ]
    else:
        cmd += ['--ignore-umi']

    print(f"\n  Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"\n  ERROR: umi_tools dedup failed:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)

    subprocess.run(['samtools', 'index', out_bam], check=True)


# ---------------------------------------------------------------------------
# Multi-mapper EM assignment
# ---------------------------------------------------------------------------

def extract_multimappers(bam_path: str, out_path: str, threads: int = 4) -> int:
    """
    Write reads with NH:i: > 1 to out_path (excludes unmapped).

    Returns number of multi-mapped reads written.
    """
    cmd = (
        f"samtools view -h -F 4 '{bam_path}' | "
        f"awk '/^@/ || /\\tNH:i:[2-9][0-9]*/ || /\\tNH:i:[1-9][0-9]+/' | "
        f"samtools view -bS -@ {threads} -o '{out_path}'"
    )
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(['samtools', 'index', out_path], check=True)
    return bam_read_count(out_path)


def build_pileup_from_bam(bam_path: str) -> dict:
    """
    Build a per-position read-count dictionary from a BAM.

    Returns:
        pileup: dict keyed by (chrom, pos, strand) -> int count
            pos = 0-based leftmost position of alignment
            strand: '+' if not reverse complement, '-' if reverse complement

    Uses samtools view to avoid requiring pysam at runtime.
    """
    pileup = defaultdict(int)
    result = subprocess.run(
        ['samtools', 'view', '-F', '4', bam_path],
        capture_output=True, text=True
    )
    for line in result.stdout.splitlines():
        if not line:
            continue
        fields = line.split('\t')
        chrom = fields[2]
        pos   = int(fields[3]) - 1   # SAM is 1-based → 0-based
        flag  = int(fields[1])
        strand = '-' if (flag & 16) else '+'
        pileup[(chrom, pos, strand)] += 1
    return pileup


def em_assign_multimappers(multi_bam:   str,
                           pileup:      dict,
                           out_bam:     str,
                           pseudocount: float = 0.1,
                           window:      int   = 50,
                           threads:     int   = 4,
                           **kwargs) -> int:
    """
    Single-pass positional rescue of multi-mapped reads.

    For each read name, scores every candidate locus by summing the
    unique-dedup pileup within a ±window bp neighbourhood plus a
    pseudocount, then hard-assigns the read to the highest-scoring locus.

    Ties are broken by preferring the primary alignment (FLAG & 256 == 0).
    The winning alignment is re-tagged NH:i:1 so downstream umi_tools dedup
    treats it as a unique mapper.

    Requires pysam. Returns number of reads written.

    Args:
        multi_bam   : input BAM of multi-mapped reads (NH:i:>1)
        pileup      : dict (chrom, pos, strand) -> int — unique dedup pileup
        out_bam     : output BAM (one hard-assigned alignment per read name)
        pseudocount : floor score added to every locus (avoid zero weights)
        window      : half-width (bp) of pileup scoring window
        threads     : pysam / samtools threads
        **kwargs    : accepts (and ignores) max_iter / epsilon for API compat
    """
    if not HAS_PYSAM:
        print(
            "  ERROR: pysam is required for --multi-map mode.\n"
            "  Install with: conda install -c bioconda pysam",
            file=sys.stderr
        )
        sys.exit(1)

    import pysam as _pysam

    # Group all alignments per query name
    read_groups: dict = defaultdict(list)
    with _pysam.AlignmentFile(multi_bam, 'rb', threads=threads) as bam_in:
        header = bam_in.header.to_dict()
        for read in bam_in:
            if read.is_unmapped:
                continue
            read_groups[read.query_name].append(read)

    n_written = 0
    with _pysam.AlignmentFile(out_bam, 'wb',
                               header=header, threads=threads) as bam_out:
        for qname, alignments in read_groups.items():
            best_score = -1.0
            best_aln   = None

            for aln in alignments:
                chrom  = aln.reference_name
                pos    = aln.reference_start        # 0-based
                strand = '-' if aln.is_reverse else '+'

                score = pseudocount
                for p in range(max(0, pos - window), pos + window + 1):
                    score += pileup.get((chrom, p, strand), 0)

                # Prefer primary alignment on tie
                is_primary = not (aln.flag & 256)
                if (score > best_score) or \
                   (score == best_score and is_primary and
                    best_aln is not None and (best_aln.flag & 256)):
                    best_score = score
                    best_aln   = aln

            if best_aln is None:
                continue

            best_aln.set_tag('NH', 1)
            best_aln.flag &= ~256   # clear not-primary-alignment
            best_aln.flag &= ~2048  # clear supplementary
            bam_out.write(best_aln)
            n_written += 1

    # Sort and index
    sorted_tmp = out_bam + '.sorted.bam'
    subprocess.run(
        ['samtools', 'sort', '-@', str(threads), '-o', sorted_tmp, out_bam],
        check=True
    )
    os.replace(sorted_tmp, out_bam)
    subprocess.run(['samtools', 'index', out_bam], check=True)

    return n_written


# ---------------------------------------------------------------------------
# Core deduplication
# ---------------------------------------------------------------------------

def run_dedup(bam_path:   str,
              out_path:   str,
              umi_len:    int,
              max_nh:     int   = 1,
              multi_map:  bool  = False,
              threads:    int   = 4,
              log_path:   str   = None,
              tmp_dir:    str   = None,
              umi_tools:  str   = None) -> None:
    """
    Run umi_tools dedup on a sorted, indexed BAM.

    Multimapper filtering (max_nh) is applied before deduplication so that
    reads with NH > max_nh never enter the UMI graph. Default max_nh=1 keeps
    only uniquely mapped reads, matching CTK's parseAlignment.pl --map-qual 1
    behaviour.

    When multi_map=True, multi-mapped reads (NH:i:>1) are rescued via an
    EM-style assignment step (CLAM-inspired):
      1. Unique reads (NH:i:1) → umi_tools dedup → unique.dedup.bam
      2. Multi-mapped reads (NH:i:>1) → EM assignment using unique.dedup.bam
         pileup as prior → hard-assign to best locus → umi_tools dedup
         → multi.dedup.bam
      3. samtools merge → combined.dedup.bam (written to out_path)

    Args:
        bam_path  : input sorted BAM
        out_path  : output deduplicated BAM
        umi_len   : UMI length in bases (0 = position-only, no UMI)
        max_nh    : maximum NH tag value to keep when multi_map=False
                    (1 = unique only, 0 = disable; default: 1).
                    Ignored when multi_map=True (unique reads are always
                    extracted first; all multi-mappers are EM-assigned).
        multi_map : if True, rescue multi-mapped reads via single-pass
                    positional assignment (default: False)
        threads   : samtools threads for sort/index steps
        log_path  : path for umi_tools log (default: <out>.log)
        tmp_dir   : temp directory for intermediate files
        umi_tools : path to umi_tools binary (auto-detected if None)
    """
    if umi_tools is None:
        umi_tools = find_umi_tools()
        if umi_tools is None:
            print(
                "  ERROR: umi_tools not found.\n"
                "  Install with: conda create -n umi_tools -c bioconda umi_tools python=3.11",
                file=sys.stderr
            )
            sys.exit(1)

    if log_path is None:
        stem = Path(out_path).stem
        # Strip trailing _dedup to avoid SRR_dedup_dedup.log
        base = stem[:-6] if stem.endswith('_dedup') else stem
        log_path = str(Path(out_path).parent / base) + '_dedup.log'

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Clink collapse", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"  Input:   {bam_path}", file=sys.stderr)
    print(f"  Output:  {out_path}", file=sys.stderr)
    print(f"  UMI len: {umi_len} {'(position-only mode)' if umi_len == 0 else 'nt'}", file=sys.stderr)
    if multi_map:
        print(f"  Mode:    multi-map rescue (EM assignment)", file=sys.stderr)
    else:
        print(f"  max NH:  {max_nh} {'(unique mappers only)' if max_nh == 1 else '(filter disabled)' if max_nh == 0 else ''}", file=sys.stderr)
    print(f"  Log:     {log_path}", file=sys.stderr)

    # Ensure indexed
    check_bam_index(bam_path)

    # Count input alignments (a multi-mapped read contributes one record per locus)
    n_in = bam_read_count(bam_path)
    print(f"\n  Total input alignments: {n_in:,}", file=sys.stderr)

    use_tmp = tmp_dir or str(Path(out_path).parent)
    stem    = Path(out_path).stem

    # =========================================================================
    # MULTI-MAP MODE: EM-based rescue of multi-mapped reads
    # =========================================================================
    if multi_map:
        # --- Step 1: extract unique reads (NH:i:1) ---------------------------
        tmp_unique_bam = str(Path(use_tmp) / f"{stem}_unique.bam")
        print(f"\n  [multi-map] Extracting unique reads (NH:i:1)...", file=sys.stderr)
        n_unique = filter_multimappers(bam_path, tmp_unique_bam,
                                       max_nh=1, threads=threads)
        print(f"  [multi-map] Unique reads: {n_unique:,}", file=sys.stderr)

        # --- Step 2: dedup unique reads → unique.dedup.bam -------------------
        tmp_unique_dedup = str(Path(use_tmp) / f"{stem}_unique.dedup.bam")
        tmp_unique_log   = str(Path(use_tmp) / f"{stem}_unique_dedup.log")
        _run_umi_tools_dedup(umi_tools, tmp_unique_bam, tmp_unique_dedup,
                              tmp_unique_log, umi_len)
        n_unique_dedup = bam_read_count(tmp_unique_dedup)
        print(f"  [multi-map] After dedup (unique): {n_unique_dedup:,}", file=sys.stderr)

        # --- Step 3: extract multi-mapped reads (NH:i:>1) --------------------
        tmp_multi_bam = str(Path(use_tmp) / f"{stem}_multi.bam")
        print(f"\n  [multi-map] Extracting multi-mapped reads (NH:i:>1)...", file=sys.stderr)
        n_multi = extract_multimappers(bam_path, tmp_multi_bam, threads=threads)
        print(f"  [multi-map] Multi-mapped reads: {n_multi:,}", file=sys.stderr)

        n_out = n_unique_dedup   # default if no multimappers

        if n_multi > 0:
            # --- Step 4: build pileup from unique.dedup.bam ------------------
            print(f"\n  [multi-map] Building pileup from unique dedup reads...", file=sys.stderr)
            pileup = build_pileup_from_bam(tmp_unique_dedup)
            print(f"  [multi-map] Pileup sites: {len(pileup):,}", file=sys.stderr)

            # --- Step 5: EM-assign multi-mappers → assigned.bam --------------
            tmp_assigned_bam = str(Path(use_tmp) / f"{stem}_multi.assigned.bam")
            print(f"\n  [multi-map] EM-assigning multi-mapped reads...", file=sys.stderr)
            n_assigned = em_assign_multimappers(
                tmp_multi_bam, pileup, tmp_assigned_bam, threads=threads
            )
            print(f"  [multi-map] Assigned reads: {n_assigned:,}", file=sys.stderr)

            # --- Step 6: dedup assigned reads --------------------------------
            tmp_multi_dedup = str(Path(use_tmp) / f"{stem}_multi.dedup.bam")
            tmp_multi_log   = str(Path(use_tmp) / f"{stem}_multi_dedup.log")
            _run_umi_tools_dedup(umi_tools, tmp_assigned_bam, tmp_multi_dedup,
                                  tmp_multi_log, umi_len)
            n_multi_dedup = bam_read_count(tmp_multi_dedup)
            print(f"  [multi-map] After dedup (multi): {n_multi_dedup:,}", file=sys.stderr)

            # --- Step 7: merge unique.dedup + multi.dedup → out_path ---------
            print(f"\n  [multi-map] Merging unique + multi dedup BAMs...", file=sys.stderr)
            tmp_merged = str(Path(use_tmp) / f"{stem}_merged.bam")
            subprocess.run(
                ['samtools', 'merge', '-f', '-@', str(threads),
                 tmp_merged, tmp_unique_dedup, tmp_multi_dedup],
                check=True
            )
            subprocess.run(
                ['samtools', 'sort', '-@', str(threads), '-o', out_path, tmp_merged],
                check=True
            )
            subprocess.run(['samtools', 'index', out_path], check=True)

            n_out = bam_read_count(out_path)

            # Clean up intermediates
            for f in [tmp_multi_bam, tmp_multi_bam + '.bai',
                      tmp_assigned_bam, tmp_assigned_bam + '.bai',
                      tmp_multi_dedup, tmp_multi_dedup + '.bai',
                      tmp_multi_log, tmp_merged]:
                try:
                    os.remove(f)
                except FileNotFoundError:
                    pass
        else:
            # No multi-mappers — just rename unique.dedup to out_path
            print(f"  [multi-map] No multi-mapped reads found; using unique dedup only.",
                  file=sys.stderr)
            import shutil as _shutil
            _shutil.copy2(tmp_unique_dedup, out_path)
            _shutil.copy2(tmp_unique_dedup + '.bai', out_path + '.bai')

        # Clean up unique intermediates
        for f in [tmp_unique_bam, tmp_unique_bam + '.bai',
                  tmp_unique_dedup, tmp_unique_dedup + '.bai',
                  tmp_unique_log]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

        pct_kept = 100 * n_out / n_in if n_in > 0 else 0
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Deduplication complete (multi-map mode)", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)
        print(f"  Total input alignments:  {n_in:>10,}", file=sys.stderr)
        print(f"    Unique (NH:i:1):       {n_unique:>10,}", file=sys.stderr)
        print(f"    Multi  (NH:i:>1):      {n_multi:>10,}", file=sys.stderr)
        print(f"  After dedup (total):     {n_out:>10,}  ({pct_kept:.1f}% of input retained)", file=sys.stderr)
        print(f"  Log:                     {log_path}", file=sys.stderr)
        return   # done

    # =========================================================================
    # STANDARD MODE: unique-only (NH:i:<=max_nh) deduplication
    # =========================================================================

    # --- Multimapper filter (NH > max_nh) ------------------------------------
    dedup_input = bam_path   # may be replaced by filtered temp BAM

    if max_nh >= 1:
        tmp_filtered = str(Path(use_tmp) / (stem + '_nhfiltered.bam'))
        print(f"  Filtering NH > {max_nh} reads...", file=sys.stderr)
        n_filtered = filter_multimappers(bam_path, tmp_filtered,
                                         max_nh=max_nh, threads=threads)
        n_removed_nh = n_in - n_filtered
        print(f"  After NH filter: {n_filtered:,}  ({n_removed_nh:,} multimappers removed)",
              file=sys.stderr)
        dedup_input = tmp_filtered
    else:
        print(f"  NH filter disabled (max_nh=0) — all reads passed to dedup",
              file=sys.stderr)
        n_filtered = n_in

    _run_umi_tools_dedup(umi_tools, dedup_input, out_path, log_path, umi_len)

    # Clean up temp filtered BAM
    if max_nh >= 1 and dedup_input != bam_path:
        for f in [dedup_input, dedup_input + '.bai']:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    # Count output reads and report
    n_out = bam_read_count(out_path)
    n_pcr_dups = n_filtered - n_out
    pct_kept = 100 * n_out / n_in if n_in > 0 else 0

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Deduplication complete", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"  Total input alignments:  {n_in:>10,}", file=sys.stderr)
    if max_nh >= 1:
        print(f"  After NH filter:         {n_filtered:>10,}  ({n_in - n_filtered:,} multimappers removed)", file=sys.stderr)
    print(f"  After dedup:             {n_out:>10,}  ({pct_kept:.1f}% of input retained)", file=sys.stderr)
    print(f"  PCR dups removed:        {n_pcr_dups:>10,}", file=sys.stderr)
    print(f"  Log:                     {log_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Clink collapse: UMI-aware PCR deduplication via umi_tools dedup')

    parser.add_argument('--bam', required=True,
        help='Sorted, indexed input BAM')
    parser.add_argument('--out', required=True,
        help='Output deduplicated BAM')
    parser.add_argument('--umi-len', type=int, default=-1,
        help='UMI length in bases. 0 = position-only (no UMI). '
             'Default: auto-detect from read names.')
    parser.add_argument('--max-nh', type=int, default=1,
        help='Maximum NH tag value to keep (1 = unique mappers only, '
             '0 = disable filter and keep all reads; default: 1). '
             'Ignored when --multi-map is set.')
    parser.add_argument('--multi-map', action='store_true', default=False,
        help='Rescue multi-mapped reads (NH:i:>1) via single-pass positional '
             'assignment. Unique reads are deduped first; multi-mappers are '
             'hard-assigned to the candidate locus with the greatest unique-read '
             'pileup support (±50 bp window), then deduped and merged. '
             'Requires pysam. Default: off (unique-only).')
    parser.add_argument('--threads', type=int, default=4,
        help='samtools threads (default: 4)')
    parser.add_argument('--log', default=None,
        help='Path for umi_tools log (default: <out>.log)')
    parser.add_argument('--umi-tools', default=None,
        help='Path to umi_tools binary (default: auto-detect)')

    args = parser.parse_args()

    # Auto-detect UMI length if not specified
    umi_len = args.umi_len
    if umi_len < 0:
        print("  Auto-detecting UMI from read names...", file=sys.stderr)
        has_umi, umi_len = detect_umi_in_bam(args.bam)
        if has_umi:
            print(f"  Detected UMI length: {umi_len}nt", file=sys.stderr)
        else:
            print(f"  No UMI detected — using position-only mode", file=sys.stderr)
            umi_len = 0

    run_dedup(
        bam_path  = args.bam,
        out_path  = args.out,
        umi_len   = umi_len,
        max_nh    = args.max_nh,
        multi_map = args.multi_map,
        threads   = args.threads,
        log_path  = args.log,
        umi_tools = args.umi_tools,
    )
