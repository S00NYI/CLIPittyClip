# Manual Installation Guide

This guide provides step-by-step instructions for manually installing CLIPittyClip and its dependencies if the automated installation scripts don't work for your system.

## Prerequisites

- **Conda/Mamba**: [Miniforge](https://github.com/conda-forge/miniforge) recommended
- **Git**: For cloning repositories
- **curl** or **wget**: For downloading files
- **Build tools**: `gcc`, `make` (Linux: `build-essential`; macOS: Xcode Command Line Tools)

## Step 1: Clone Repository

```bash
git clone -b v3-development https://github.com/S00NYI/CLIPittyClip.git
cd CLIPittyClip
```

## Step 2: Create Conda Environment

### Linux

```bash
mamba create -n clipittyclip -c conda-forge -c bioconda \
  perl>=5.32 bedtools samtools>=1.15 star>=2.7 bowtie2 \
  cutadapt fastp seqkit python>=3.10 pandas numpy scipy
```

### macOS (Apple Silicon and Intel)

Use a native conda, e.g. [Miniforge](https://github.com/conda-forge/miniforge) (arm64 build on Apple Silicon). Rosetta is not needed; do not set `CONDA_SUBDIR`.

```bash
mamba create -n clipittyclip --override-channels -c conda-forge -c bioconda \
  perl bedtools "samtools>=1.15" star=2.7.11b bowtie2 bwa \
  cutadapt fastp seqkit "python>=3.10,<3.12" pandas numpy scipy pysam umi_tools \
  perl-bioperl-core perl-math-cdf
```

Use `perl-bioperl-core`, not `perl-bioperl`: CTK only needs `Bio::SeqIO`, and the full package pins samtools to 0.1.19. STAR 2.7.10b has no arm64 build; 2.7.11b reads indices built with 2.7.10b.

## Step 3: Install Perl Modules via CPAN (Linux)

> **macOS:** skip this step. `Bio::SeqIO` and `Math::CDF` were installed from conda in Step 2 (nothing to compile). Verify with the two `perl -M` commands below.

```bash
conda activate clipittyclip

# Install cpanminus if not available
curl -L https://cpanmin.us | perl - App::cpanminus

# Install required modules
cpanm --notest Math::CDF
cpanm --notest --force XML::LibXML  # May need --force due to test failures
cpanm --notest Bio::SeqIO

# Verify
perl -MMath::CDF -e 'print "Math::CDF OK\n"'
perl -MBio::SeqIO -e 'print "Bio::SeqIO OK\n"'
```

## Step 4: Install CTK

```bash
mkdir -p ~/Tools

# Clone CTK
git clone https://github.com/chaolinzhanglab/ctk.git ~/Tools/ctk

# Clone czplib (required Perl library)
git clone https://github.com/chaolinzhanglab/czplib.git ~/Tools/ctk/czplib

# Make scripts executable
chmod +x ~/Tools/ctk/*.pl
```

## Step 5: Install HOMER

```bash
mkdir -p ~/Tools/homer && cd ~/Tools/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install homer
```

HOMER compiles from source, so it is built for the architecture it was installed under. If it was installed by an older Intel/Rosetta setup, rebuild in place (keeps downloaded genomes): `perl ~/Tools/homer/configureHomer.pl -make`.

## Step 6: Configure Shell Environment

Add to `~/.zshrc` (macOS) or `~/.bashrc` (Linux):

```bash
# CLIPittyClip
export PATH="$PATH:/path/to/CLIPittyClip"

# CTK (CLIP Tool Kit)
export PATH="$PATH:$HOME/Tools/ctk"
export PERL5LIB="$PERL5LIB:$HOME/Tools/ctk/czplib"

# HOMER
export PATH="$PATH:$HOME/Tools/homer/bin"
```

Then reload:
```bash
source ~/.zshrc  # or ~/.bashrc
```

## Step 7: Verify Installation

```bash
conda activate clipittyclip

which CLIPittyClip.sh
which parseAlignment.pl
which findPeaks
perl -MBio::SeqIO -e 'print "Bio::SeqIO OK\n"'
perl -MMath::CDF -e 'print "Math::CDF OK\n"'
```

## Common Issues

### "Can't locate MyConfig.pm"
CTK requires the czplib library. Make sure you cloned it:
```bash
git clone https://github.com/chaolinzhanglab/czplib.git ~/Tools/ctk/czplib
```

### XML::LibXML test failures
Use `--force` flag: `cpanm --notest --force XML::LibXML`
