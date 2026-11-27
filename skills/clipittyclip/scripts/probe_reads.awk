# probe_reads.awk: measure the read layout from the first reads of a FASTQ file.
#
# Usage:
#   zcat -f READS.fastq.gz | head -400000 | awk -f probe_reads.awk
#   zcat -f READS.fastq.gz | head -400000 | awk -v bc=barcodes.txt -f probe_reads.awk
#
# The script prints measurements only. It does not pick a protocol.
#   1. Reads with a UMI in the header (eclipdemux format, for eCLIP pe).
#   2. Reads that contain each known 3' adapter (12-nt core, left-most match).
#   3. Adapter dimers: reads where the adapter core sits at the smallest common
#      position. In a dimer, all bases before the adapter come from the primers.
#   4. For each base before the adapter in dimers: the top base and its fraction.
#      Mask: N = random (top base < 80%), a letter = fixed base.
#   5. With -v bc=FILE (Name<TAB>Sequence): reads that match a barcode at each
#      offset 0-20, and the per-barcode counts at the best offset.

BEGIN {
    FIXED = 0.80
    MAXPOS = 60
    nad = split("L32 L3 L3MOD L7 SRA3 FPAR", an, " ")
    seq["L32"] = "GTGTCAGTCACTTCCAGCGG";   core["L32"] = "CAGTCACTTCCA"
    seq["L3"] = "TGAGATCGGAAGAGCGGTTCAG";  core["L3"] = "GGAAGAGCGGTT"
    seq["L3MOD"] = "AGATCGGAAGAGCGGTTCAG"; core["L3MOD"] = "GGAAGAGCGGTT"
    seq["L7"] = "AGATCGGAAGAGCACACGTC";    core["L7"] = "GGAAGAGCACAC"
    seq["SRA3"] = "TCGTATGCCGTCTTCTGCTTG"; core["SRA3"] = "TATGCCGTCTTC"
    seq["FPAR"] = "TGGAATTCTCGGGTGCCAAGG"; core["FPAR"] = "GAATTCTCGGGT"
    nbc = 0
    if (bc != "") {
        while ((getline line < bc) > 0) {
            sub(/\r$/, "", line)
            if (line == "" || line ~ /^#/) continue
            split(line, f, "\t")
            if (f[2] == "") continue
            nbc++; bname[f[2]] = f[1]; blen = length(f[2])
        }
        close(bc)
    }
}

NR % 4 == 1 {
    h = substr($1, 2); split(h, hp, ":")
    if (length(hp[1]) >= 5 && length(hp[1]) <= 10 && hp[1] ~ /^[ACGTN]+$/) hdr_umi++
}

NR % 4 == 2 {
    n++; s[n] = $0
    best = ""; bp = 0
    for (i = 1; i <= nad; i++) {
        a = an[i]; if (a == "L3MOD") continue    # L3 and L3MOD share a core
        p = index($0, core[a])
        if (p > 0 && (bp == 0 || p < bp)) { best = a; bp = p }
    }
    if (best != "") { hit[best]++; cpos[best, bp]++; rbest[n] = best; rpos[n] = bp }
    if (nbc > 0)
        for (o = 0; o <= 20; o++) {
            b = substr($0, o + 1, blen)
            if (b in bname) { bco[o]++; bcount[o, bname[b]]++ }
        }
}

END {
    if (n == 0) { print "No reads found."; exit 1 }
    printf "Reads sampled: %d\n", n
    printf "Reads with a UMI in the header: %d (%.1f%%)\n\n", hdr_umi, 100 * hdr_umi / n

    print "3' adapter core found (left-most match per read):"
    top = ""
    for (i = 1; i <= nad; i++) {
        a = an[i]; if (a == "L3MOD") continue
        name = (a == "L3") ? "L3/L3MOD" : a
        printf "  %-9s %8d  (%5.1f%%)\n", name, hit[a], 100 * hit[a] / n
        if (hit[a] > hit[top] || top == "") top = a
    }
    if (hit[top] == 0) { print "\nNo known adapter found. Ask the user for the adapter sequence."; exit 0 }

    # Dimer core position: the smallest position that holds >= 1% of the reads.
    dp = 0
    for (p = 1; p <= MAXPOS; p++) if (cpos[top, p] >= 0.01 * n) { dp = p; break }
    printf "\nAdapter core positions for %s (base where the core starts, top counts):\n", top
    for (p = 1; p <= MAXPOS; p++) if (cpos[top, p] >= 0.005 * n)
        printf "  base %3d: %8d  (%5.1f%%)\n", p, cpos[top, p], 100 * cpos[top, p] / n
    if (dp == 0) { print "\nNo adapter-dimer peak (no position holds 1% of the reads)."; exit 0 }

    nd = 0
    for (r = 1; r <= n; r++) if (rbest[r] == top && rpos[r] == dp) d[++nd] = s[r]

    # Walk left from the core while >= 80% of the dimers match the adapter.
    ad = (top == "L3") ? "L3MOD" : top    # walk L3MOD, then check for the L3 "TG"
    off = index(seq[ad], core[ad]) - 1
    k = 0
    while (k < off && dp - k - 1 >= 1) {
        want = substr(seq[ad], off - k, 1); m = 0
        for (j = 1; j <= nd; j++) if (substr(d[j], dp - k - 1, 1) == want) m++
        if (m / nd < FIXED) break
        k++
    }
    start = dp - k
    if (top == "L3") {
        # L3 starts with TG. L3MOD does not.
        tg = 0
        for (j = 1; j <= nd; j++) if (substr(d[j], start - 2, 2) == "TG") tg++
        ad = (start > 2 && tg / nd >= FIXED) ? "L3" : "L3MOD"
        if (ad == "L3") start -= 2
    }
    plen = start - 1
    printf "\nAdapter: %s (%s)\n", ad, seq[ad]
    printf "Adapter dimers: %d reads (%.1f%%). The adapter starts at base %d.\n", nd, 100 * nd / n, start
    printf "Bases before the adapter in dimers: %d\n\n", plen

    if (plen > 0) {
        print "Base  top  fraction  mask"
        mask = ""
        for (q = 1; q <= plen; q++) {
            split("", bcnt); tb = ""
            for (j = 1; j <= nd; j++) bcnt[substr(d[j], q, 1)]++
            for (b in bcnt) if (tb == "" || bcnt[b] > bcnt[tb]) tb = b
            fr = bcnt[tb] / nd
            mk = (fr >= FIXED) ? tb : "N"
            mask = mask mk
            printf "%4d  %3s  %8.2f  %4s\n", q, tb, fr, mk
        }
        printf "\nMask (5' to 3'): %s\n", mask
    }

    if (nbc > 0) {
        printf "\nBarcode matches by offset (barcode length %d, %d barcodes):\n", blen, nbc
        bo = 0
        for (o = 0; o <= 20; o++) {
            if (bco[o] > 0) printf "  offset %2d: %8d  (%5.1f%%)\n", o, bco[o], 100 * bco[o] / n
            if (bco[o] > bco[bo]) bo = o
        }
        printf "Best offset: %d\n", bo
        for (b in bname) printf "  %-20s %8d  (%5.1f%%)\n", bname[b], bcount[bo, bname[b]], 100 * bcount[bo, bname[b]] / n
    }
}
