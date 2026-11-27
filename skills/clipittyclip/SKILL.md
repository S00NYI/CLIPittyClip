---
name: clipittyclip
description: Configure and run a CLIP-seq analysis with the CLIPittyClip pipeline (iCLIP, iCLIP2, iCLIP3, irCLIP, coCLIP, BrdU-CLIP, eCLIP, PAR-CLIP, HITS-CLIP). Interviews the user about input files, protocol and read layout, genome index, analysis tracks (HOMER peaks, Clink, CTK CIMS/CITS), grouping, and runtime options, checks the files, shows the exact CLIPittyClip.sh command, and runs it only after the user confirms. Use when the user wants to analyze CLIP-seq FASTQ files, call peaks or crosslink sites from CLIP data, or run CLIPittyClip.
---

# CLIPittyClip analysis skill

This skill is the chat version of `CLIPittyClip.sh -w` (the wizard in `lib/wizard.sh`).
You interview the user, check their files, build one `CLIPittyClip.sh` command, and run it
after the user confirms.

Reference files in this skill:

- `references/read_architectures.md`: protocol-to-arguments table (UMI, barcode, spacer,
  adapter) from the CLIP read architecture index.
- `references/options.md`: every pipeline option with its default and when to ask about it.
- `scripts/probe_reads.awk`: measures the read layout from the first reads (Round 3).

The pipeline glossary is in the install folder, not in this skill. Find it in Step 0.

## Rules

1. Ask before you act. Do not install, start, download, write, or delete anything until the
   user gives a clear yes for that action. One yes covers one action only.
   Read-only checks (`ls`, `head`, `zcat | head`, the read probe) do not need a yes.
   Do not start the pipeline before the user confirms the final command.
2. Ask questions in rounds. Use `AskUserQuestion` for choices. Ask for file paths in plain
   chat text, because paths do not fit a multiple-choice question.
3. If a value has a pipeline default, offer the default as the first option. Mark it
   "(Recommended)" only when the default is correct for the user's protocol.
4. If you can find an answer from the files, find it and show the evidence. Do not ask the
   user for a fact the data shows (read length, file count, barcode length, index type).
5. Pass the full adapter sequence to `-a`, so the command shows the exact sequence. Older
   CLIPittyClip versions (3.5.2 and earlier) do not translate adapter names. If the adapter
   is L32, omit `-a`.
6. Use the user's words for sample and group names. Do not rename files.
7. After the run, report facts from the logs (read counts, mapping rates, peak counts,
   file paths). Do not interpret the biology.
8. Trust but verify, in every step. Treat each fact the user gives as a claim, and check
   it against the files before you use it. If you cannot check a claim, say what you could
   not check. If the files and the claim differ, show both and let the user decide.
   Examples:
   - Protocol and read layout: labs change protocols, so an "eCLIP" or "iCLIP2" library
     can differ from the published layout. Check it against the reads (Round 3).
   - Paths: check that each file or folder exists and has the expected contents.
   - Genome index: check the index type (STAR or Bowtie2) from its files.
   - Barcode and groups files: check the format, and that the names match the samples.
   - Genome FASTA and GTF: check that the chromosome names match the index.
   - Install: run the checks after the installer, and do not trust its exit message alone.
   - Run results: read the logs, and do not trust a success message alone.
9. Use the terms and explanations from the CLIPittyClip glossary (`docs/GLOSSARY.md` in
   the install folder) as the preferred wording. Read the entry before you explain a term,
   and adapt it to the user's question. The glossary guides the wording, but it does not
   limit what you can explain. Do not tell the user whether a term is in the glossary.

## Step 0: Preflight

Run these checks before the first question:

```bash
command -v CLIPittyClip.sh || conda run -n clipittyclip which CLIPittyClip.sh
nproc 2>/dev/null || sysctl -n hw.ncpu
free -g 2>/dev/null | awk 'NR==2{print $2" GB RAM"}' || echo "$(( $(sysctl -n hw.memsize) / 1073741824 )) GB RAM"
```

- If `CLIPittyClip.sh` is on PATH, call it directly.
- Find the glossary. The installer puts the CLIPittyClip folder on PATH, so the glossary is
  at `$(dirname "$(command -v CLIPittyClip.sh)")/docs/GLOSSARY.md`.
- If it is only in the conda environment, put `conda run -n <env> --no-capture-output` in front
  of every pipeline command.
- If neither check finds it, go to "Install" below.
- On macOS, if the user has not run it before, offer `lib/star_selftest.sh`. Stock STAR
  2.7.11b on macOS maps zero reads and still reports success.

Keep the CPU count and RAM. You use them in Round 5.

### Install

Do not install without a clear yes. First tell the user that CLIPittyClip is not installed.
Then ask with one `AskUserQuestion` call:

- "Install CLIPittyClip now"
- "Show me the install steps, I will install it"
- "Stop"

If the user picks "Show me the install steps", show the steps below and stop.
If the user picks "Install CLIPittyClip now":

1. Check for conda: `command -v mamba || command -v conda`. If neither exists, tell the
   user to install Miniforge first. Then stop.
2. Check for an existing environment: `conda env list | grep -w clipittyclip`. If it
   exists, ask whether to use another name (`--env <name>`). The installer asks
   interactively to remove an existing environment, so do not run it on that name.
3. Ask for the clone folder and the tools folder (default `~/Tools`, for CTK and HOMER).
4. Show the exact commands and what they change. The installer creates a conda
   environment, downloads CTK and HOMER into the tools folder, and can edit the shell
   startup file. Ask: "Run these install commands?" Continue only on a clear yes.
5. Run the commands:
   ```bash
   git clone https://github.com/LunaRNALab/CLIPittyClip.git <clone_folder>/CLIPittyClip
   cd <clone_folder>/CLIPittyClip
   ./install_linux.sh --env clipittyclip --tools-dir <tools_folder>   # macOS: ./install_macos.sh
   ```
   The installer can run for 10 minutes or more. Run it in the background and read its
   output. Do not poll in a tight loop.
6. If the installer fails, quote the error line and stop. Do not retry without the user.
7. Run the README checks: `conda run -n clipittyclip which CLIPittyClip.sh`, `which
   findPeaks`, `which parseAlignment.pl`. Report each result. Then go back to Step 0.

## Round 1: Input and protocol

Ask with one `AskUserQuestion` call:

- **Input type**: "Pooled FASTQ + barcode file (demultiplex first)", "One demultiplexed
  FASTQ", "Folder of demultiplexed FASTQs".
- **Protocol family**: "iCLIP family (iCLIP, iCLIP2, iCLIP3, irCLIP, irCLIP2)",
  "coCLIP or BrdU-CLIP", "eCLIP (seCLIP or peCLIP)", "PAR-CLIP or HITS-CLIP".
  "Other" covers a custom layout and "I do not know".

Then ask one follow-up `AskUserQuestion` for the exact variant inside the family. Use the
variant list in `references/read_architectures.md`.

- If the user names a variant, keep it as the claimed layout. Round 3 checks it.
- If the user has a custom layout, do the part-by-part interview in Round 3.
- If the user does not know the protocol, do not guess. Round 3 measures the layout
  from the reads first.

## Round 2: File paths

Ask in plain text for the paths that apply:

| Path | Ask when |
|---|---|
| FASTQ file or folder (`-i` or `-d`) | Always |
| Barcode file (`-b`), format `Name<TAB>Sequence` | Input type is pooled |
| Genome index folder (`-x`), STAR or Bowtie2 | Always |
| Output name or path (`-o`) | Always. A blank answer means `<input>_output` next to the input. |

Check each path after the user gives it:

```bash
ls -la <path>
```

- Index type: a STAR index folder has `SA`, `Genome`, and `chrNameLength.txt`. A Bowtie2
  folder has `*.1.bt2` or `*.1.bt2l`. If you find Bowtie2 files, set `-m bowtie2`.
- Folder input: list the FASTQ files (`*.fastq.gz`, `*.fq.gz`, `*.fastq`, `*.fq`) and show
  the sample names to the user. CLIPittyClip 3.5.2 and earlier accept only `*.fastq.gz` and
  `*.fq.gz` with `-d`. On those versions, tell the user to compress plain files with `gzip`.
- Barcode file: show the rows. Run `check_barcodes.sh -f <file> -m 1` from the
  CLIPittyClip repo to confirm that barcodes stay unique with 1 mismatch.
- If a path does not exist or has the wrong contents, say so and ask again.

## Round 3: Verify the read layout

Do this round for every run, also when the user names a known protocol.

### 3a. Confirm the claimed layout part by part

Skip this step if the user does not know the protocol.

Show the claimed layout as a table, one row per part, 5' to 3'. Take the values from
`references/read_architectures.md`, or from the interview below for a custom layout.

| Part | Length (nt) | Sequence | Argument |
|---|---|---|---|
| UMI | 7 | random (`N`) | `-u 7` |
| Barcode | 6 | from the barcode file | `--bc-len 6` |
| Spacer | 3 | random (`DDD`) | `--spacer-len 3` |
| Insert (RNA) | varies | | |
| Bases before the adapter | 0 | | `--trim3 0` |
| 3' adapter | 20 | `GTGTCAGTCACTTCCAGCGG` (L32) | default |

Ask with one `AskUserQuestion` call: "Does your library match this layout?"
Use these options:
"Yes, it matches", "Some parts differ", "I am not sure". If the user picks "Some parts
differ", do the part-by-part interview. "I am not sure" is a valid answer. Step 3b
measures the layout.

**Part-by-part interview.** Ask about one part at a time, 5' to 3'. Use
`AskUserQuestion` for choices and plain text for numbers and sequences:

1. What comes first in the read: the UMI, the barcode, or the RNA insert (no 5' tag)?
2. The UMI length.
3. The barcode length, and whether the barcode is before or after the UMI.
4. The spacer after the UMI and barcode: its length, and whether it is random or fixed
   (give the sequence if fixed).
5. The bases between the insert and the 3' adapter (for example, a second UMI): the length.
6. The 3' adapter: a name from the adapter table or the sequence.

If the user does not know a part, write "unknown" and continue. Step 3b measures it.

### 3b. Measure the layout from the reads

Show the first two reads:

```bash
zcat -f <fastq> | head -8
```

Run the read probe on the first 100,000 reads. It is in this skill, at
`scripts/probe_reads.awk`. For pooled input, add `-v bc=<barcode_file>`:

```bash
zcat -f <fastq> | head -400000 | awk -v bc=<barcode_file> -f <skill_dir>/scripts/probe_reads.awk
```

For folder input, run the probe on two or three files.

The probe reports:

- The reads with a UMI in the header (eCLIP pe after `eclipdemux`).
- The reads that contain each known 3' adapter.
- The adapter dimers: reads where the adapter starts at the smallest common position.
- For each base before the adapter in the dimers, the top base and its fraction. The mask
  shows `N` for a random base and the letter for a fixed base.
- With a barcode file: the reads that match a barcode at each offset from 0 to 20.

How to read the probe output:

- The dimer prefix length is the total of all 5' parts plus the bases before the
  adapter: UMI + barcode + spacer + `--trim3` bases.
- In pooled data, the barcode positions show as `N`, because the barcode changes from
  read to read. The best barcode offset shows where the barcode starts. In a
  demultiplexed file, the barcode shows as fixed letters.
- A fixed spacer shows as letters. A random spacer shows as `N`.
- The probe cannot split a block of `N` bases into UMI and spacer, or into 5' UMI and
  bases before the adapter. Only the protocol or the user can say where one part ends.

### 3c. Compare and decide

Show a table with three columns: the part, the claimed value, and the measured value.
Compare the adapter, the dimer prefix length, the fixed bases, the barcode offset, and
(for eCLIP pe) the header UMI.

- If all parts agree, say so and give the numbers. Then continue.
- If a part differs, show the difference with the numbers. Ask with one
  `AskUserQuestion` call: "Use the measured layout", "Keep my layout", "Let me describe
  the layout". Do not change the layout without the user.
- If the user did not know the protocol, compare the measured values with every row of
  `references/read_architectures.md`. Show the rows that fit, and the evidence for each.
  If one row fits, ask the user to confirm it. If more than one row fits, or none, show
  what the reads cannot tell (for example, the UMI and spacer split), and ask the user.

The probe needs adapter dimers. If no position holds at least 1% of the reads, the probe
cannot measure the prefix. Tell the user, and ask for the primer sequences or the
protocol paper. Do not guess the layout.

## Round 4: Analysis plan

Ask with one `AskUserQuestion` call:

- **Crosslink-site tracks** (multiSelect): "Clink (Python, BH FDR, recommended)",
  "CTK CITS (truncation sites)", "CTK CIMS (mutation sites)", "None (peaks only)".
  HOMER or CTK peaks always run.
- **Peak caller**: "HOMER (default)", "CTK tag2peak".
- **Grouping** (only for pooled or folder input): "No groups", "Groups file for averaged
  tracks and peak matrix", "Groups file and pooled crosslink-site calls (--group-xlsite)".

If the user picked any crosslink-site track, ask in plain text for the genome FASTA path
(`--genome-fasta`). The pipeline uses it for `samtools calmd`. The user can skip it, but
say that the help text calls it "strongly recommended".

Check that the FASTA chromosome names match the index. Compare the first names from
`grep '>' <fasta> | head` with `head <index>/chrNameLength.txt` (STAR). If they differ,
for example `chr1` against `1`, show both and ask the user.

If the user picked a groups file, ask for its path. The format is `SampleName<TAB>GroupName`.
Check that each sample name matches a sample from Round 2. The pipeline accepts names with
or without the `.fastq.gz` extension. Samples that are not in the file run alone.

Protocol notes for this round:

- PAR-CLIP: recommend Clink. The primary signal is T>C substitution. CITS gives few or no
  sites for PAR-CLIP, because truncated cDNAs are not cloned.
- `--group-xlsite` with CTK can need more than 64 GB RAM. If the RAM from Step 0 is below
  that, tell the user before they pick it.

## Round 5: Settings

Ask with one `AskUserQuestion` call:

- **Threads**: offer the CPU count from Step 0, half of it, and the default (1).
- **Settings depth**: "Use defaults (Recommended)", "Review one category",
  "Walk through all options".
- **Test run first**: "Yes, 100,000 reads with -s 100000 (Recommended)", "No, run all reads".
- **Where to run**: "Run here now", "Give me the command (for a cluster job or tmux)".

If the user picks "Review one category" or "Walk through all options", use
`references/options.md`. Show only the categories that apply to this run. Ask at most four
options per `AskUserQuestion` call.

Ask about these options only when they apply:

- More than about 30 million reads per sample, or low RAM: offer `--low-memory`.
- Clink is on: offer `--xl-bigwig` and `--clink-multi-map`.
- The user mentions rRNA, snRNA, or other abundant RNA: explain the two options.
  `--filter-repeat` removes reads before alignment and needs a repeat index in the genome
  index folder. `--report-abundant-rna` only reports fractions and needs `--gtf`.
- Scaffold-level genome with no `chr1`-style names: the run needs `--no-chr-filter`.

## Final confirmation

Show three things:

1. A summary table: input, protocol, genome index, aligner, tracks, groups, output path.
2. The options that differ from the defaults, with the default next to each.
3. The full command in a `bash` code block, one option per line with `\`.

Offer to write the command to `clipittyclip_command_<YYYYMMDD>.sh` in the current folder.
With this file, the user can run the command again without the interview. Write the file
only on a yes.

Then ask: "Start the analysis with this command?" Start only on a clear yes.

## Run

If the user picked "Give me the command", stop after the confirmation and give the file path.

If the user picked "Run here now":

1. Start the pipeline with `nohup`, so the run continues when the session closes:
   ```bash
   nohup <command> > clipittyclip_<YYYYMMDD_HHMM>.log 2>&1 &
   echo $!
   ```
2. Tell the user the PID and the log path.
3. Check progress with `tail -n 30 <log>`. Do not poll in a tight loop. A full run can take
   hours.
4. If the log shows `[ERROR]`, quote the error line. Find the cause in the log, then
   propose one fix. Do not restart the pipeline without the user.

If the user picked a test run, report the test results first. Then ask whether to run the
same command without `-s`.

## Report

After the run ends, report facts only:

- Output folder path and the folders it contains (see the README "Output Structure").
- Per-sample input reads, reads after trimming, and uniquely mapped reads, from
  `00_REPORTS/ALIGNER_LOGS/` and `00_REPORTS/SAMPLES/`.
- Peak count from `04_PEAKS/COMBINED_PEAKS/`.
- Crosslink-site counts per sample from the Clink or CTK folders, if those tracks ran.
- Warnings from the log.

Do not judge whether the experiment worked. The user draws the conclusions.

## Cite

At the end, give the citation from the README:

Yi S, Luna JM. CLIPittyClip: a single-command framework for CLIP sequencing data analysis.
NAR Genomics and Bioinformatics. 2026;8(4):lqag138. https://doi.org/10.1093/nargab/lqag138
