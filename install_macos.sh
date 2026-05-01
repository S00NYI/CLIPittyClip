#!/usr/bin/env bash
#===============================================================================
# CLIPittyClip macOS Installation Script
# 
# Self-contained installation script for macOS (Intel and Apple Silicon).
# Creates conda environment, installs Perl dependencies via CPAN, and
# configures CTK and HOMER from source.
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
    - Requires conda or mamba to be installed
    - Uses x86 emulation via Rosetta 2 for Apple Silicon compatibility
    - Installs CTK and HOMER from source (not available via conda on macOS)
    - Installs Perl dependencies (BioPerl, Math::CDF) via CPAN

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

# Check for conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    print_success "Found mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    print_success "Found conda"
else
    print_error "Neither conda nor mamba found. Please install Miniconda or Mambaforge first."
    echo -e "\n  Installation instructions: https://github.com/conda-forge/miniforge"
    exit 1
fi

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

# Check for Xcode Command Line Tools (needed for CPAN compilation)
if ! xcode-select -p &> /dev/null; then
    print_warning "Xcode Command Line Tools not found."
    print_info "Installing Xcode Command Line Tools..."
    xcode-select --install
    echo ""
    print_info "Please complete the Xcode installation popup, then re-run this script."
    exit 1
fi
print_success "Found Xcode Command Line Tools"

# Check for Homebrew and install required libraries
if command -v brew &> /dev/null; then
    print_success "Found Homebrew"
    
    # Check and install libxml2 (needed for XML::LibXML)
    if ! brew list libxml2 &> /dev/null; then
        print_info "Installing libxml2 via Homebrew..."
        brew install libxml2 || print_warning "Failed to install libxml2"
    else
        print_success "libxml2 already installed"
    fi
    
    # Check and install openssl (needed for Net::SSLeay)
    if ! brew list openssl &> /dev/null; then
        print_info "Installing openssl via Homebrew..."
        brew install openssl || print_warning "Failed to install openssl"
    else
        print_success "openssl already installed"
    fi

    # Check and install expat (needed for XML::Parser)
    if ! brew list expat &> /dev/null; then
        print_info "Installing expat via Homebrew..."
        brew install expat || print_warning "Failed to install expat"
    else
        print_success "expat already installed"
    fi

    # Check and install berkeley-db (needed for DB_File)
    if ! brew list berkeley-db &> /dev/null; then
        print_info "Installing berkeley-db via Homebrew..."
        brew install berkeley-db || print_warning "Failed to install berkeley-db"
    else
        print_success "berkeley-db already installed"
    fi
else
    print_warning "Homebrew not found. Some Perl modules may fail to compile."
    print_info "Install Homebrew from: https://brew.sh"
fi

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
# Step 3: Create conda environment with x86 architecture
#-------------------------------------------------------------------------------
if [[ -z "$SKIP_CONDA" ]]; then
    print_step "Creating conda environment '${ENV_NAME}' (x86 architecture)..."
    print_info "Using x86 emulation for Rosetta 2 compatibility"
    print_info "This may take several minutes..."

    # Create environment with x86 architecture
    CONDA_SUBDIR=osx-64 $CONDA_CMD create -n "$ENV_NAME" -y

    # Configure environment for x86
    conda activate "$ENV_NAME" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "$ENV_NAME"
    conda config --env --set subdir osx-64

    # Install packages (excluding perl-bioperl and perl-math-cdf which aren't available on macOS)
    print_info "Installing conda packages..."
    $CONDA_CMD install -n "$ENV_NAME" -y \
        -c conda-forge -c bioconda \
        wget \
        "python>=3.10,<3.13" \  # umi_tools incompatible with Python 3.13
        pysam \
        umi_tools \
        perl \
        perl-threaded \
        perl-yaml \
        bedtools \
        samtools \
        htslib \
        bowtie2 \
        bwa \
        "star=2.7.10b" \   # 2.7.11b breaks subprocess spawning on macOS Tahoe via Rosetta
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
        ca-certificates \
        openssl \
        certifi \
        clang_osx-64 \
        clangxx_osx-64

    if [[ $? -eq 0 ]]; then
        print_success "Conda packages installed successfully"
    else
        print_error "Failed to install conda packages"
        exit 1
    fi
fi

#-------------------------------------------------------------------------------
# Step 4: Install Perl dependencies via CPAN
#-------------------------------------------------------------------------------
print_step "Installing Perl dependencies via CPAN..."

# Activate environment to use its perl
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Create compiler wrapper scripts
# The conda perl package references non-existent compilers (x86_64-apple-darwin13.4.0-clang)
# These wrapper scripts redirect to system clang/ar/ranlib/ld
print_info "Creating compiler wrappers for CPAN..."
CONDA_BIN="$(conda info --base)/envs/$ENV_NAME/bin"

# Remove existing symlinks first (conda creates symlinks to clang-21 which interfere)
rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang"
rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++"
rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-ar"
rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-ranlib"
rm -f "$CONDA_BIN/x86_64-apple-darwin13.4.0-ld"

cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang" << 'CLANG_WRAPPER'
#!/bin/bash
# Add -Wno-incompatible-pointer-types to fix XML::LibXML compilation on newer clang
exec /usr/bin/clang -Wno-incompatible-pointer-types "$@"
CLANG_WRAPPER

cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++" << 'CLANGXX_WRAPPER'
#!/bin/bash
exec /usr/bin/clang++ -Wno-incompatible-pointer-types "$@"
CLANGXX_WRAPPER

cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-ar" << 'WRAPPER'
#!/bin/bash
exec /usr/bin/ar "$@"
WRAPPER

cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-ranlib" << 'WRAPPER'
#!/bin/bash
exec /usr/bin/ranlib "$@"
WRAPPER

cat > "$CONDA_BIN/x86_64-apple-darwin13.4.0-ld" << 'WRAPPER'
#!/bin/bash
exec /usr/bin/ld "$@"
WRAPPER

chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang"
chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-clang++"
chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-ar"
chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-ranlib"
chmod +x "$CONDA_BIN/x86_64-apple-darwin13.4.0-ld"
print_success "Compiler wrappers created"

# FIX: Create xlocale.h compatibility stub
# Conda's perl 5.32 was built expecting xlocale.h, which was removed in
# macOS Catalina (10.15) and later. This stub redirects to locale.h which
# now contains the same definitions. The existence check prevents overwriting
# on reinstall.
if [[ ! -f "$CONDA_PREFIX/include/xlocale.h" ]]; then
    print_info "Creating xlocale.h compatibility stub..."
    echo '#include <locale.h>' > "$CONDA_PREFIX/include/xlocale.h"
    print_success "xlocale.h stub created"
fi

# Set compiler flags for Homebrew libraries
if [[ -d "/opt/homebrew/opt/openssl" ]]; then
    # Apple Silicon
    export LDFLAGS="-L/opt/homebrew/opt/openssl/lib -L/opt/homebrew/opt/libxml2/lib -L/opt/homebrew/opt/expat/lib -L/opt/homebrew/opt/berkeley-db/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/openssl/include -I/opt/homebrew/opt/libxml2/include -I/opt/homebrew/opt/expat/include -I/opt/homebrew/opt/berkeley-db/include"
    export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl/lib/pkgconfig:/opt/homebrew/opt/libxml2/lib/pkgconfig:/opt/homebrew/opt/expat/lib/pkgconfig"
elif [[ -d "/usr/local/opt/openssl" ]]; then
    # Intel Mac
    export LDFLAGS="-L/usr/local/opt/openssl/lib -L/usr/local/opt/libxml2/lib -L/usr/local/opt/expat/lib -L/usr/local/opt/berkeley-db/lib"
    export CPPFLAGS="-I/usr/local/opt/openssl/include -I/usr/local/opt/libxml2/include -I/usr/local/opt/expat/include -I/usr/local/opt/berkeley-db/include"
    export PKG_CONFIG_PATH="/usr/local/opt/openssl/lib/pkgconfig:/usr/local/opt/libxml2/lib/pkgconfig:/usr/local/opt/expat/lib/pkgconfig"
fi

# Check if cpanm is available, if not install it
if ! command -v cpanm &> /dev/null; then
    print_info "Installing cpanminus..."
    curl -L https://cpanmin.us | perl - App::cpanminus 2>/dev/null || true
fi

# Install modules needed for BioPerl compatibility
print_info "Installing XML::Parser (required for Bio::SeqIO)..."
cpanm --notest XML::Parser 2>/dev/null || true

print_info "Installing DB_File (required for Bio::SeqIO)..."
cpanm --notest DB_File 2>/dev/null || true

# Install only the modules CTK actually needs (not full BioPerl)
print_info "Installing Math::CDF (required for CIMS/CITS)..."
cpanm --notest Math::CDF 2>/dev/null || true
if perl -MMath::CDF -e '1' 2>/dev/null; then
    print_success "Math::CDF installed successfully"
else
    print_warning "Math::CDF installation failed - CIMS/CITS analysis will not work"
fi

# XML::LibXML has minor test failures but works - install with --force
print_info "Installing XML::LibXML (dependency for Bio::SeqIO)..."
cpanm --notest --force XML::LibXML 2>/dev/null || true

print_info "Installing Bio::SeqIO (required for sequence operations)..."
cpanm --notest Bio::SeqIO 2>/dev/null || true
if perl -MBio::SeqIO -e '1' 2>/dev/null; then
    print_success "Bio::SeqIO installed successfully"
else
    print_warning "Bio::SeqIO installation failed - some CTK features may not work"
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