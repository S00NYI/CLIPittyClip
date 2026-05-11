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
"""

import sys
import os
import argparse
import subprocess
import shutil
import tempfile
from pathlib import Path


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
# Core deduplication
# ---------------------------------------------------------------------------

def run_dedup(bam_path:   str,
              out_path:   str,
              umi_len:    int,
              max_nh:     int  = 1,
              threads:    int  = 4,
              log_path:   str  = None,
              tmp_dir:    str  = None,
              umi_tools:  str  = None) -> None:
    """
    Run umi_tools dedup on a sorted, indexed BAM.

    Multimapper filtering (max_nh) is applied before deduplication so that
    reads with NH > max_nh never enter the UMI graph. Default max_nh=1 keeps
    only uniquely mapped reads, matching CTK's parseAlignment.pl --map-qual 1
    behaviour.

    Args:
        bam_path  : input sorted BAM
        out_path  : output deduplicated BAM
        umi_len   : UMI length in bases (0 = position-only, no UMI)
        max_nh    : maximum NH tag value to keep (1 = unique only, 0 = disable filter)
        threads   : samtools threads for sort/index steps
        log_path  : path for umi_tools log (default: out_path + '.log')
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
    print(f"  max NH:  {max_nh} {'(unique mappers only)' if max_nh == 1 else '(filter disabled)' if max_nh == 0 else ''}", file=sys.stderr)
    print(f"  Log:     {log_path}", file=sys.stderr)

    # Ensure indexed
    check_bam_index(bam_path)

    # Count input reads
    n_in = bam_read_count(bam_path)
    print(f"\n  Input reads: {n_in:,}", file=sys.stderr)

    # --- Multimapper filter (NH > max_nh) ------------------------------------
    use_tmp = tmp_dir or str(Path(out_path).parent)
    dedup_input = bam_path   # may be replaced by filtered temp BAM

    if max_nh >= 1:
        tmp_filtered = str(Path(use_tmp) / (Path(out_path).stem + '_nhfiltered.bam'))
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

    # Build umi_tools command
    cmd = [
        umi_tools, 'dedup',
        '-I', dedup_input,
        '-S', out_path,
        '-L', log_path,
        '--method', 'directional',
        '--spliced-is-unique',          # reads with different splicing = distinct molecules
    ]

    if umi_len > 0:
        cmd += [
            '--extract-umi-method', 'read_id',
            '--umi-separator', '#',     # UMI is last '#'-delimited token
        ]
    else:
        cmd += ['--ignore-umi']         # position-only dedup

    print(f"\n  Running: {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Clean up temp filtered BAM
    if max_nh >= 1 and dedup_input != bam_path:
        for f in [dedup_input, dedup_input + '.bai']:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    if result.returncode != 0:
        print(f"\n  ERROR: umi_tools dedup failed:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)

    # Index output
    subprocess.run(['samtools', 'index', out_path], check=True)

    # Count output reads and report
    n_out = bam_read_count(out_path)
    n_pcr_dups = n_filtered - n_out
    pct_kept = 100 * n_out / n_in if n_in > 0 else 0

    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Deduplication complete", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"  Input reads:       {n_in:>10,}", file=sys.stderr)
    if max_nh >= 1:
        print(f"  After NH filter:   {n_filtered:>10,}  ({n_in - n_filtered:,} multimappers removed)", file=sys.stderr)
    print(f"  After dedup:       {n_out:>10,}  ({pct_kept:.1f}% of input retained)", file=sys.stderr)
    print(f"  PCR dups removed:  {n_pcr_dups:>10,}", file=sys.stderr)
    print(f"  Log:               {log_path}", file=sys.stderr)


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
             '0 = disable filter and keep all reads; default: 1)')
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
        threads   = args.threads,
        log_path  = args.log,
        umi_tools = args.umi_tools,
    )
