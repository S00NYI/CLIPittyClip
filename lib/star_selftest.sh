#!/usr/bin/env bash
#===============================================================================
# star_selftest.sh -- does this STAR binary actually align reads?
#
# Usage:  lib/star_selftest.sh [path/to/STAR]      (default: STAR on PATH)
# Exit:   0 = aligns correctly, 1 = does not, 2 = STAR not found
#
# Why this exists: STAR 2.7.11b built with libc++ (every macOS build) exits
# "successfully" but reads 0 reads from any input (upstream STAR issue #2632).
# `STAR --version` passes, so version checks cannot catch it. This runs a real,
# tiny alignment (a 100 kb random genome and 2,000 error-free 30-mers, about a
# second) and checks every read was read and mapped. See lib/patches/README.md.
#
# Needs only STAR and awk. Works in a temp directory; touches nothing else.
#===============================================================================
set -u

STAR_BIN="${1:-$(command -v STAR || true)}"
if [[ -z "$STAR_BIN" || ! -x "$STAR_BIN" ]]; then
    echo "star_selftest: STAR not found (${STAR_BIN:-not on PATH})" >&2
    exit 2
fi

N_READS=2000
READ_LEN=30
GENOME_LEN=100000
MIN_UNIQUE=$(( N_READS * 95 / 100 ))   # error-free reads from a random genome: expect ~100%

WORK="$(mktemp -d "${TMPDIR:-/tmp}/star_selftest.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK" || exit 2

# Synthetic genome + reads cut from it (fixed seed => identical every run).
awk -v n="$N_READS" -v rl="$READ_LEN" -v gl="$GENOME_LEN" 'BEGIN {
    srand(42); split("A C G T", b, " ")
    g = ""
    for (c = 0; c < gl / 1000; c++) {          # build in chunks: avoids slow 1-char appends
        chunk = ""
        for (i = 0; i < 1000; i++) chunk = chunk b[int(rand() * 4) + 1]
        g = g chunk
    }
    print ">chr1" > "genome.fa"
    for (i = 1; i <= gl; i += 60) print substr(g, i, 60) > "genome.fa"
    q = ""; for (i = 0; i < rl; i++) q = q "I"
    for (r = 0; r < n; r++) {
        p = int(rand() * (gl - rl)) + 1
        printf "@read%d\n%s\n+\n%s\n", r, substr(g, p, rl), q > "reads.fq"
    }
}'

mkdir idx
if ! "$STAR_BIN" --runMode genomeGenerate --genomeDir idx --genomeFastaFiles genome.fa \
        --genomeSAindexNbases 6 --genomeChrBinNbits 12 --outFileNamePrefix gen_ > gen.out 2>&1; then
    echo "star_selftest: FAIL ($STAR_BIN could not even build a tiny genome index)" >&2
    tail -5 gen.out >&2
    exit 1
fi

if ! "$STAR_BIN" --runThreadN 2 --genomeDir idx --readFilesIn reads.fq \
        --outFileNamePrefix aln_ --outTmpDir tmp_aln > aln.out 2>&1; then
    echo "star_selftest: FAIL ($STAR_BIN exited with an error while aligning)" >&2
    tail -5 aln.out >&2
    exit 1
fi

field() { awk -F'|' -v k="$1" 'index($0, k) { gsub(/[ \t]/, "", $2); print $2; exit }' aln_Log.final.out; }
n_in="$(field 'Number of input reads')"
n_uniq="$(field 'Uniquely mapped reads number')"
version="$("$STAR_BIN" --version 2>/dev/null)"

if [[ "${n_in:-0}" -eq "$N_READS" && "${n_uniq:-0}" -ge "$MIN_UNIQUE" ]]; then
    echo "star_selftest: PASS (STAR $version: read $n_in/$N_READS reads, $n_uniq mapped uniquely)"
    exit 0
fi

echo "star_selftest: FAIL (STAR $version: read ${n_in:-?}/$N_READS reads, ${n_uniq:-?} mapped uniquely)" >&2
if [[ "${n_in:-0}" -eq 0 ]]; then
    echo "  STAR read zero reads. This is the known STAR 2.7.11b / libc++ defect (macOS); see lib/patches/README.md" >&2
fi
exit 1
