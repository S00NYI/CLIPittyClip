#!/usr/bin/env python3
"""
fastq_collapse_hash.py - Hash-based FASTQ exact-duplicate collapser
Part of CLIPittyClip lib/dedup.sh module

Replaces: fastq2collapse.pl (CTK) sort | uniq pipeline
Method:   Two-pass collections.Counter — O(n) time, no disk spill

Output format: identical to fastq2collapse.pl
  @readname#COUNT
  SEQUENCE
  +
  QUALITY

Usage:
  python3 fastq_collapse_hash.py <input.fastq> <output.fastq>

  input.fastq  - plain (uncompressed) FASTQ file; gzip handled by caller
  output.fastq - plain (uncompressed) FASTQ file; gzip handled by caller

Memory estimate:
  ~9 GB for 100M unique 30-nt sequences (Counter only)
  ~12 GB peak (Counter + emitted set)

Exit codes:
  0 - success
  1 - argument / file error
"""

import sys
import os
from collections import Counter


def fastq_records(path):
    """Yield (header_line, seq, plus_line, qual) tuples from a plain FASTQ file."""
    with open(path, "r") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq    = fh.readline().rstrip("\n")
            plus   = fh.readline()
            qual   = fh.readline().rstrip("\n")
            if not qual:
                print(f"[ERROR] Truncated FASTQ record near: {header.strip()}", file=sys.stderr)
                sys.exit(1)
            yield header.rstrip("\n"), seq, plus.rstrip("\n"), qual


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} <input.fastq> <output.fastq>", file=sys.stderr)
        sys.exit(1)

    input_path  = sys.argv[1]
    output_path = sys.argv[2]

    if not os.path.isfile(input_path):
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # ── Pass 1: count sequence occurrences ────────────────────────────────────
    print("[INFO] fastq_collapse_hash.py: Pass 1 — counting sequences...", file=sys.stderr)
    counts = Counter()
    for header, seq, plus, qual in fastq_records(input_path):
        counts[seq] += 1

    total_reads  = sum(counts.values())
    unique_seqs  = len(counts)
    print(f"[INFO] fastq_collapse_hash.py: {total_reads:,} reads, {unique_seqs:,} unique sequences",
          file=sys.stderr)

    if total_reads == 0:
        print("[ERROR] fastq_collapse_hash.py: input file is empty", file=sys.stderr)
        sys.exit(1)

    # ── Pass 2: emit first occurrence with count appended to read ID ──────────
    # Output format mirrors fastq2collapse.pl exactly:
    #   header stored by awk '{print $1}' includes the leading '@'
    #   output: $2"#"$1 = @original_id#COUNT
    print("[INFO] fastq_collapse_hash.py: Pass 2 — writing deduplicated output...", file=sys.stderr)
    emitted = set()

    with open(output_path, "w") as out:
        for header, seq, plus, qual in fastq_records(input_path):
            if seq in emitted:
                continue
            emitted.add(seq)

            # header is e.g. "@READ1 comment" — take first whitespace-delimited token
            id_token = header.split()[0]          # "@READ1"
            count    = counts[seq]
            out.write(f"{id_token}#{count}\n{seq}\n+\n{qual}\n")

    written = len(emitted)
    print(f"[INFO] fastq_collapse_hash.py: wrote {written:,} deduplicated records to {output_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
