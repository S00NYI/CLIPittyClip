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

### macOS (Apple Silicon)

macOS requires x86 emulation via Rosetta 2:

```bash
# Create x86 environment
CONDA_SUBDIR=osx-64 mamba create -n clipittyclip
conda activate clipittyclip
conda config --env --set subdir osx-64

# Install packages
mamba install -c conda-forge -c bioconda \
  perl>=5.32 bedtools samtools>=1.15 star>=2.7 bowtie2 \
  cutadapt fastp seqkit python>=3.10 pandas numpy scipy \
  libxml2 openssl
```

## Step 3: Install Perl Modules via CPAN

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

### Troubleshooting macOS CPAN Installation

If CPAN modules fail to compile on macOS:

1. **Install Xcode Command Line Tools:**
   ```bash
   xcode-select --install
   ```

2. **Install Homebrew dependencies:**
   ```bash
   brew install libxml2 openssl
   ```

3. **Create compiler wrappers** (for conda's Perl):
   ```bash
   CONDA_BIN=$(conda info --base)/envs/clipittyclip/bin
   
   # Remove existing symlinks
   rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang"
   rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++"
   
   # Create wrapper scripts
   cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang" << 'EOF'
   #!/bin/bash
   exec /usr/bin/clang -Wno-incompatible-pointer-types "$@"
   EOF
   
   cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++" << 'EOF'
   #!/bin/bash
   exec /usr/bin/clang++ -Wno-incompatible-pointer-types "$@"
   EOF
   
   chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang"
   chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++"
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

### "x86_64-apple-darwin13.4.0-clang not found"
See Step 3 troubleshooting section for creating compiler wrappers.

### XML::LibXML test failures
Use `--force` flag: `cpanm --notest --force XML::LibXML`
