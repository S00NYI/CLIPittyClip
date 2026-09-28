#!/usr/bin/env bash
#===============================================================================
# CLIPittyClip macOS Installation Script
# 
# Self-contained installation script for macOS (Apple Silicon and Intel).
# Builds a NATIVE conda environment (no Rosetta), installs the Perl
# dependencies CTK needs from conda, and configures CTK and HOMER.
#
# Usage:
#   ./install_macos.sh [OPTIONS]
#
# Options:
#   --env <name>         Conda environment name (default: clipittyclip)
#   --tools-dir <path>   Directory for CTK/HOMER (default: ~/Tools)
#   --help               Show this help message
#===============================================================================

set -e

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
ENV_NAME="clipittyclip"
TOOLS_DIR="$HOME/Tools"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors (with fallback for non-interactive terminals)
if [[ -t 1 ]]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    NC='\033[0m' # No Color
    BOLD='\033[1m'
else
    RED='' GREEN='' YELLOW='' BLUE='' CYAN='' NC='' BOLD=''
fi

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------
print_header() {
    echo -e "\n${BOLD}${CYAN}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║${NC}        ${BOLD}CLIPittyClip macOS Installation Script${NC}               ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}╚════════════════════════════════════════════════════════════════╝${NC}\n"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_success() {
    echo -e "${GREEN}  ✓${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}  ⚠${NC} $1"
}

print_error() {
    echo -e "${RED}  ✗${NC} $1"
}

print_info() {
    echo -e "${CYAN}  ℹ${NC} $1"
}

show_help() {
    cat << EOF
CLIPittyClip macOS Installation Script

Usage: ./install_macos.sh [OPTIONS]

Options:
    --env <name>         Conda environment name (default: clipittyclip)
    --tools-dir <path>   Directory for CTK/HOMER installation (default: ~/Tools)
    --help               Show this help message

Examples:
    ./install_macos.sh
    ./install_macos.sh --env myenv --tools-dir ~/Software

Notes:
    - Requires a native conda or mamba (Miniforge: arm64 on Apple Silicon)
    - Builds a native environment; Rosetta 2 is not needed or used
    - Installs CTK and HOMER from source (not available via conda on macOS)
    - Installs Perl dependencies (BioPerl core, Math::CDF) from conda

EOF
    exit 0
}

#-------------------------------------------------------------------------------
# Parse Arguments
#-------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --env)
            ENV_NAME="$2"
            shift 2
            ;;
        --tools-dir)
            TOOLS_DIR="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Expand ~ in TOOLS_DIR if present
TOOLS_DIR="${TOOLS_DIR/#\~/$HOME}"

#-------------------------------------------------------------------------------
# Main Installation
#-------------------------------------------------------------------------------
# Check if running on Linux
if [[ "$(uname)" != "Darwin" ]]; then
    print_error "This script is for macOS. Please run ./install_linux.sh instead."
    exit 1
fi

# arm64 on Apple Silicon, x86_64 on Intel. Everything is built for this, natively.
ARCH="$(uname -m)"

print_header

echo -e "${BOLD}Configuration:${NC}"
echo -e "  Environment name: ${CYAN}${ENV_NAME}${NC}"
echo -e "  Tools directory:  ${CYAN}${TOOLS_DIR}${NC}"
echo -e "  Script directory: ${CYAN}${SCRIPT_DIR}${NC}"
echo ""

#-------------------------------------------------------------------------------
# Step 1: Check Prerequisites
#-------------------------------------------------------------------------------
print_step "Checking prerequisites..."

# Find a conda/mamba that actually RUNS. `command -v` is not enough: an Intel-only
# Anaconda earlier in PATH fails with "Bad CPU type in executable" on Macs without
# Rosetta, so probe each candidate (PATH first, then common Miniforge locations).
CONDA_CMD=""
for cand in mamba conda \
            "$HOME/miniforge3/bin/mamba" "$HOME/miniforge3/bin/conda" \
            "$HOME/mambaforge/bin/mamba" "$HOME/mambaforge/bin/conda" \
            "$HOME/miniconda3/bin/conda" \
            /opt/homebrew/Caskroom/miniforge/base/bin/mamba; do
    exe="$(command -v "$cand" 2>/dev/null)" || continue
    if "$exe" --version &> /dev/null; then
        CONDA_CMD="$exe"
        break
    fi
    print_warning "Skipping $exe (cannot run on this Mac; Intel-only install?)"
done

if [[ -z "$CONDA_CMD" ]]; then
    print_error "No working conda or mamba found."
    echo -e "\n  Install Miniforge (native ${ARCH} build): https://github.com/conda-forge/miniforge"
    echo -e "  or: brew install --cask miniforge"
    exit 1
fi
print_success "Using $(basename "$CONDA_CMD"): $CONDA_CMD"

# Make later plain `conda ...` calls resolve to the working install too.
export PATH="$(dirname "$CONDA_CMD"):$PATH"

# Check for git
if ! command -v git &> /dev/null; then
    print_error "git is required but not found. Please install git first."
    exit 1
fi
print_success "Found git"

# Check for curl
if ! command -v curl &> /dev/null; then
    print_error "curl is required but not found. Please install curl first."
    exit 1
fi
print_success "Found curl"

# Check for Xcode Command Line Tools (provides git, and the compiler HOMER is built with)
if ! xcode-select -p &> /dev/null; then
    print_warning "Xcode Command Line Tools not found."
    print_info "Installing Xcode Command Line Tools..."
    xcode-select --install
    echo ""
    print_info "Please complete the Xcode installation popup, then re-run this script."
    exit 1
fi
print_success "Found Xcode Command Line Tools"

#-------------------------------------------------------------------------------
# Step 2: Check if environment exists
#-------------------------------------------------------------------------------
print_step "Checking for existing environment..."

SKIP_CONDA=""
if conda env list | grep -q "^${ENV_NAME} "; then
    print_warning "Environment '${ENV_NAME}' already exists."
    read -p "  Do you want to remove and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
    else
        print_info "Keeping existing environment. Skipping conda installation."
        SKIP_CONDA=true
    fi
fi

#-------------------------------------------------------------------------------
# Step 3: Create native conda environment
#-------------------------------------------------------------------------------
if [[ -z "$SKIP_CONDA" ]]; then
    print_step "Creating conda environment '${ENV_NAME}' (native ${ARCH})..."
    print_info "This may take several minutes..."

    # Never force a platform: conda picks osx-arm64 / osx-64 to match this Mac.
    # --override-channels keeps ~/.condarc (e.g. the Anaconda 'defaults' channel,
    # which also triggers its terms-of-service prompt) out of the solve.
    #
    # star=2.7.11b: 2.7.10b has no osx-arm64 build. 2.7.11b reads indices built
    # by 2.7.10b (index format is compatible since 2.7.4a).
    #
    # perl-bioperl-core, not perl-bioperl: CTK only needs Bio::SeqIO. The full
    # metapackage pulls in ~120 extra packages and pins samtools to 0.1.19.
    # perl-math-cdf is needed for CIMS/CITS. Both ship prebuilt, so nothing is
    # compiled from CPAN.
    print_info "Installing conda packages..."
    if ! $CONDA_CMD create -n "$ENV_NAME" -y --override-channels \
        -c conda-forge -c bioconda \
        wget \
        "python>=3.10,<3.12" \
        perl \
        perl-threaded \
        perl-yaml \
        perl-bioperl-core \
        perl-math-cdf \
        bedtools \
        ucsc-bedgraphtobigwig \
        samtools \
        htslib \
        bowtie2 \
        bwa \
        "star=2.7.11b" \
        cutadapt \
        fastp \
        seqkit \
        trim-galore \
        pandas \
        numpy \
        scipy \
        setuptools \
        seaborn \
        matplotlib \
        pysam \
        umi_tools \
        ca-certificates \
        openssl \
        certifi; then
        print_error "Failed to create conda environment"
        exit 1
    fi
    print_success "Conda environment created"
fi

#-------------------------------------------------------------------------------
# Step 4: Verify the environment is native and its Perl modules load
#-------------------------------------------------------------------------------
print_step "Verifying environment..."

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# An env built by an older (Intel/Rosetta) version of this installer cannot run on
# a Mac without Rosetta. Catch it here rather than as a cryptic failure mid-run.
if ! file "$CONDA_PREFIX/bin/python" | grep -q "$ARCH"; then
    print_error "Environment '${ENV_NAME}' is not native to this Mac (${ARCH})."
    print_info "Re-run this script and answer 'y' to remove and recreate it."
    exit 1
fi
print_success "Environment is native (${ARCH})"

# Loading each module is the only check that proves the installed build works
# with THIS perl.
PERL_BROKEN=()
for probe in Math::CDF Bio::SeqIO; do
    if perl -M"$probe" -e '1' 2>/dev/null; then
        print_success "  $probe loads"
    else
        print_warning "  $probe does NOT load"
        PERL_BROKEN+=("$probe")
    fi
done

if [[ ${#PERL_BROKEN[@]} -gt 0 ]]; then
    print_warning "Perl modules unavailable: ${PERL_BROKEN[*]}"
    print_warning "CTK CIMS/CITS (--run-cims-cits) will not work. Clink (--run-clink) is unaffected."
fi

#-------------------------------------------------------------------------------
# Step 5: Install CTK
#-------------------------------------------------------------------------------
print_step "Installing CTK (CLIP Tool Kit)..."

mkdir -p "$TOOLS_DIR"
CTK_DIR="$TOOLS_DIR/ctk"

if [[ -d "$CTK_DIR" ]]; then
    print_warning "CTK already exists at $CTK_DIR"
    read -p "  Do you want to remove and reinstall? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$CTK_DIR"
    else
        print_info "Keeping existing CTK installation"
        SKIP_CTK=true
    fi
fi

if [[ -z "$SKIP_CTK" ]]; then
    print_info "Cloning CTK repository..."
    git clone https://github.com/chaolinzhanglab/ctk.git "$CTK_DIR" 2>/dev/null
    
    print_info "Cloning czplib (Perl library for CTK)..."
    git clone https://github.com/chaolinzhanglab/czplib.git "$CTK_DIR/czplib" 2>/dev/null
    
    # Ensure MyConfig.pm exists (sometimes missing from czplib repo)
    if [[ ! -f "$CTK_DIR/czplib/MyConfig.pm" ]]; then
        print_info "Downloading MyConfig.pm..."
        curl -s -o "$CTK_DIR/czplib/MyConfig.pm" \
            https://raw.githubusercontent.com/chaolinzhanglab/czplib/master/MyConfig.pm
    fi
    
    print_success "CTK installed to $CTK_DIR"
fi

#-------------------------------------------------------------------------------
# Step 6: Install HOMER
#-------------------------------------------------------------------------------
print_step "Installing HOMER..."

HOMER_DIR="$TOOLS_DIR/homer"

if [[ -d "$HOMER_DIR" ]]; then
    print_warning "HOMER already exists at $HOMER_DIR"
    read -p "  Do you want to remove and reinstall? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$HOMER_DIR"
    else
        print_info "Keeping existing HOMER installation"
        SKIP_HOMER=true
    fi
fi

if [[ -z "$SKIP_HOMER" ]]; then
    print_info "Downloading HOMER..."
    mkdir -p "$HOMER_DIR"
    cd "$HOMER_DIR"
    curl -s -O http://homer.ucsd.edu/homer/configureHomer.pl
    
    print_info "Installing HOMER (minimal configuration)..."
    perl configureHomer.pl -install 2>/dev/null
    
    cd "$SCRIPT_DIR"
    print_success "HOMER installed to $HOMER_DIR"
fi

# HOMER compiles its binaries from source, so they match the architecture of the
# shell that built them. A copy built by an older (Intel/Rosetta) run cannot
# execute on a Mac without Rosetta. Rebuild in place (~10 s) rather than
# reinstalling, which would delete any genomes already downloaded into data/.
if ! file "$HOMER_DIR/bin/findPeaks" 2>/dev/null | grep -q "$ARCH"; then
    print_info "HOMER binaries are not native (${ARCH}); rebuilding..."
    perl "$HOMER_DIR/configureHomer.pl" -make > /dev/null 2>&1 || true
    if file "$HOMER_DIR/bin/findPeaks" 2>/dev/null | grep -q "$ARCH"; then
        print_success "HOMER binaries rebuilt natively"
    else
        print_warning "HOMER rebuild failed; --peak-caller homer will not work."
        print_info "Retry manually: perl $HOMER_DIR/configureHomer.pl -make"
    fi
fi

#-------------------------------------------------------------------------------
# Step 7: Configure PATH and PERL5LIB
#-------------------------------------------------------------------------------
print_step "Configuring shell environment..."

SHELL_RC="$HOME/.zshrc"
if [[ ! -f "$SHELL_RC" ]]; then
    SHELL_RC="$HOME/.bash_profile"
fi

# FIX: Remove existing entries before re-adding to prevent duplicate PATH
# entries accumulating across reinstalls.
print_info "Removing any existing PATH entries from $SHELL_RC..."
sed -i '' '/# CLIPittyClip$/,/^$/d' "$SHELL_RC" 2>/dev/null || true
sed -i '' '/# CTK (CLIP Tool Kit)$/,/^$/d' "$SHELL_RC" 2>/dev/null || true
sed -i '' '/# HOMER$/,/^$/d' "$SHELL_RC" 2>/dev/null || true

# Add CLIPittyClip to PATH
echo "" >> "$SHELL_RC"
echo "# CLIPittyClip" >> "$SHELL_RC"
echo "export PATH=\"\$PATH:${SCRIPT_DIR}\"" >> "$SHELL_RC"
print_success "Added CLIPittyClip to PATH"

# Add CTK to PATH and PERL5LIB
echo "" >> "$SHELL_RC"
echo "# CTK (CLIP Tool Kit)" >> "$SHELL_RC"
echo "export PATH=\"\$PATH:${CTK_DIR}\"" >> "$SHELL_RC"
echo "export PERL5LIB=\"\$PERL5LIB:${CTK_DIR}/czplib\"" >> "$SHELL_RC"
print_success "Added CTK to PATH and PERL5LIB"

# Add HOMER to PATH
echo "" >> "$SHELL_RC"
echo "# HOMER" >> "$SHELL_RC"
echo "export PATH=\"\$PATH:${HOMER_DIR}/bin\"" >> "$SHELL_RC"
print_success "Added HOMER to PATH"

# Set execute permissions on CLIPittyClip scripts
chmod +x "$SCRIPT_DIR/CLIPittyClip.sh" 2>/dev/null || true
chmod +x "$SCRIPT_DIR/MAPittyMap.sh" 2>/dev/null || true
chmod +x "$SCRIPT_DIR/PEAKittyPeak.sh" 2>/dev/null || true
chmod +x "$SCRIPT_DIR/check_barcodes.sh" 2>/dev/null || true
print_success "Set execute permissions on scripts"

#-------------------------------------------------------------------------------
# Step 8: Verification
#-------------------------------------------------------------------------------
print_step "Verifying installation..."

echo ""
echo -e "${BOLD}${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${GREEN}║${NC}                  ${BOLD}Installation Complete!${NC}                       ${BOLD}${GREEN}║${NC}"
echo -e "${BOLD}${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${BOLD}Installation Summary:${NC}"
echo -e "  Conda environment: ${CYAN}${ENV_NAME}${NC}"
echo -e "  CTK location:      ${CYAN}${CTK_DIR}${NC}"
echo -e "  HOMER location:    ${CYAN}${HOMER_DIR}${NC}"
echo -e "  Shell config:      ${CYAN}${SHELL_RC}${NC}"
echo ""
echo -e "${BOLD}Next steps:${NC}"
echo -e "  1. Restart your terminal or run: ${CYAN}source ${SHELL_RC}${NC}"
echo -e "  2. Activate environment: ${CYAN}conda activate ${ENV_NAME}${NC}"
echo -e "  3. Verify installation: ${CYAN}CLIPittyClip.sh --help${NC}"
echo ""
echo -e "${BOLD}Quick verification after activation:${NC}"
echo -e "  ${CYAN}which CLIPittyClip.sh${NC}"
echo -e "  ${CYAN}which parseAlignment.pl${NC}"
echo -e "  ${CYAN}which findPeaks${NC}"
echo -e "  ${CYAN}perl -MBio::Seq -e 'print \"BioPerl OK\\n\"'${NC}"
echo -e "  ${CYAN}perl -MMath::CDF -e 'print \"Math::CDF OK\\n\"'${NC}"
echo -e "  ${CYAN}python3 -c \"import pysam; print('pysam OK')\"${NC}"
echo -e "  ${CYAN}umi_tools --version${NC}"
echo ""