#!/usr/bin/env bash
# fastq_collapse_sort.sh - FASTQ exact-duplicate collapser using only
# stock Unix tools (paste, awk, sort, uniq — no installs, no Python hash table)
#
# Same trick as CTK's fastq2collapse.pl (paste | sort | uniq -c), but merges
# the FASTQ into one-line-per-record with a single `paste - - - -` pass
# instead of three separate awk column-extraction scans, cutting full-file
# reads from 4 down to 1 before the sort/uniq stage.
#
# RAM is bounded by the system `sort`'s own default buffering behavior — no
# explicit --buffer-size/-S is passed, which lets sort fall back to its
# built-in external-merge threshold instead of forcing a large in-memory sort.
#
# Output format: identical to fastq2collapse.pl
#   @readname#COUNT
#   SEQUENCE
#   +
#   QUALITY
#
# Usage:
#   fastq_collapse_sort.sh <input.fastq> <output.fastq>

set -euo pipefail

[[ $# -ne 2 ]] && { echo "Usage: $0 <input.fastq> <output.fastq>" >&2; exit 1; }

INPUT="$1"
OUTPUT="$2"

[[ ! -f "$INPUT" ]] && { echo "[ERROR] Input not found: $INPUT" >&2; exit 1; }

TMP_DIR=$(mktemp -d "$(dirname "$OUTPUT")/.collapse_sort.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

# One pass: merge every 4 physical lines into one tab-separated record,
# trim header to its first whitespace token (Illumina headers carry a
# space-separated index/barcode comment that would otherwise be miscounted
# as an extra field by sort/uniq's default blank-splitting — BSD uniq has
# no -t/delimiter flag, so fields must be blank-free tokens instead), and
# reorder to (header, qual, seq) so seq is the LAST field for `uniq -f2`.
# Qual and seq are guaranteed blank-free (Phred33 excludes space/tab; ACGTN
# has none), so once the header is trimmed, default blank-splitting is safe.
echo "[INFO] fastq_collapse_sort.sh: merging + sorting by sequence..." >&2
paste - - - - < "$INPUT" \
    | awk -F'\t' -v OFS='\t' '{split($1,h," "); print h[1], $4, $2}' \
    | sort -k3 -T "$TMP_DIR" \
    | uniq -f2 -c \
    | awk '{print $2"#"$1"\n"$4"\n+\n"$3}' \
    > "$OUTPUT"

written=$(( $(wc -l < "$OUTPUT") / 4 ))
echo "[INFO] fastq_collapse_sort.sh: wrote ${written} deduplicated records" >&2
