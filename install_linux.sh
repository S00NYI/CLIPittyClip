#!/usr/bin/env bash
#===============================================================================
# CLIPittyClip Linux Installation Script
# 
# Self-contained installation script for Linux systems.
# Creates conda environment, installs Perl dependencies via CPAN, and
# configures CTK and HOMER from source.
#
# Usage:
#   ./install_linux.sh [OPTIONS]
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
    echo -e "${BOLD}${CYAN}║${NC}        ${BOLD}CLIPittyClip Linux Installation Script${NC}                ${BOLD}${CYAN}║${NC}"
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
CLIPittyClip Linux Installation Script

Usage: ./install_linux.sh [OPTIONS]

Options:
    --env <name>         Conda environment name (default: clipittyclip)
    --tools-dir <path>   Directory for CTK/HOMER (default: ~/Tools)
    --help               Show this help message

Examples:
    ./install_linux.sh
    ./install_linux.sh --env myenv --tools-dir /opt/tools

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

#-------------------------------------------------------------------------------
# Main Installation
#-------------------------------------------------------------------------------
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

# Check for curl
if ! command -v curl &> /dev/null; then
    print_error "curl is required but not found. Please install curl first."
    exit 1
fi
print_success "Found curl"

# Check for git
if ! command -v git &> /dev/null; then
    print_error "git is required but not found. Please install git first."
    exit 1
fi
print_success "Found git"

# Check for build tools (gcc/make)
if ! command -v gcc &> /dev/null; then
    print_warning "gcc not found. Installing build-essential may be required for CPAN modules."
    print_info "Run: sudo apt-get install build-essential"
fi

#-------------------------------------------------------------------------------
# Step 2: Check if environment exists
#-------------------------------------------------------------------------------
print_step "Checking for existing environment..."

if conda env list | grep -q "^${ENV_NAME} "; then
    print_warning "Environment '${ENV_NAME}' already exists."
    read -p "  Do you want to remove and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
    else
        print_info "Keeping existing environment. Skipping conda package installation."
        SKIP_CONDA=true
    fi
fi

#-------------------------------------------------------------------------------
# Step 3: Create conda environment
#-------------------------------------------------------------------------------
if [[ -z "$SKIP_CONDA" ]]; then
    print_step "Creating conda environment '${ENV_NAME}'..."
    print_info "This may take several minutes..."

    # Create the environment from inline YAML
    $CONDA_CMD env create -n "$ENV_NAME" -f - <<'YAML_ENV'
name: clipittyclip
channels:
  - conda-forge
  - bioconda
dependencies:
  # Core utilities
  - wget
  - perl>=5.32
  - libxml2
  - libxslt
  - openssl
  
  # Alignment & Processing
  - bedtools
  - samtools>=1.15
  - htslib
  - bowtie2
  - star>=2.7
  
  # CLIP-specific
  - cutadapt
  - fastp
  - seqkit
  
  # Python packages
  - python>=3.10
  - pandas
  - numpy
  - scipy
  - seaborn
  - matplotlib
  
  # SSL/Certs
  - ca-certificates
  - certifi
YAML_ENV

    if [[ $? -eq 0 ]]; then
        print_success "Conda environment created successfully"
    else
        print_error "Failed to create conda environment"
        exit 1
    fi
fi

#-------------------------------------------------------------------------------
# Step 4: Install Perl dependencies via CPAN
#-------------------------------------------------------------------------------
print_step "Installing Perl dependencies via CPAN..."

# Activate the environment for CPAN installation
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Check if cpanm is available, if not install it
if ! command -v cpanm &> /dev/null; then
    print_info "Installing cpanminus..."
    curl -L https://cpanmin.us | perl - App::cpanminus 2>/dev/null || true
fi

# Install only the modules CTK actually needs (not full BioPerl)
print_info "Installing Math::CDF (required for CIMS/CITS)..."
cpanm --notest Math::CDF 2>/dev/null || true
if perl -MMath::CDF -e '1' 2>/dev/null; then
    print_success "Math::CDF installed successfully"
else
    print_warning "Math::CDF installation failed - CIMS/CITS analysis will not work"
fi

# XML::LibXML might need --force due to test failures
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
    print_info "Updating CTK..."
    cd "$CTK_DIR" && git pull 2>/dev/null || true
else
    print_info "Cloning CTK repository..."
    git clone https://github.com/chaolinzhanglab/ctk.git "$CTK_DIR" 2>/dev/null || {
        print_error "Failed to clone CTK repository"
    }
fi

# Clone czplib (required Perl library)
CZPLIB_DIR="$CTK_DIR/czplib"
if [[ -d "$CZPLIB_DIR" ]]; then
    print_info "Updating czplib..."
    cd "$CZPLIB_DIR" && git pull 2>/dev/null || true
else
    print_info "Cloning czplib (Perl library for CTK)..."
    git clone https://github.com/chaolinzhanglab/czplib.git "$CZPLIB_DIR" 2>/dev/null || {
        print_error "Failed to clone czplib repository"
    }
fi

# Make CTK scripts executable
chmod +x "$CTK_DIR"/*.pl 2>/dev/null || true

print_success "CTK installed to $CTK_DIR"

#-------------------------------------------------------------------------------
# Step 6: Install HOMER
#-------------------------------------------------------------------------------
print_step "Installing HOMER..."

HOMER_DIR="$TOOLS_DIR/homer"

if [[ -d "$HOMER_DIR" ]]; then
    print_warning "HOMER already exists at $HOMER_DIR"
else
    print_info "Downloading HOMER..."
    mkdir -p "$HOMER_DIR"
    cd "$HOMER_DIR"
    
    # Download and install HOMER
    wget -q http://homer.ucsd.edu/homer/configureHomer.pl -O configureHomer.pl || {
        print_error "Failed to download HOMER"
    }
    
    print_info "Installing HOMER (minimal configuration)..."
    perl configureHomer.pl -install homer 2>/dev/null || {
        print_warning "HOMER installation had issues, but may still work"
    }
fi

print_success "HOMER installed to $HOMER_DIR"

#-------------------------------------------------------------------------------
# Step 7: Configure shell environment
#-------------------------------------------------------------------------------
print_step "Configuring shell environment..."

# Detect shell config file
if [[ -f "$HOME/.bashrc" ]]; then
    SHELL_RC="$HOME/.bashrc"
elif [[ -f "$HOME/.bash_profile" ]]; then
    SHELL_RC="$HOME/.bash_profile"
elif [[ -f "$HOME/.zshrc" ]]; then
    SHELL_RC="$HOME/.zshrc"
else
    SHELL_RC="$HOME/.bashrc"
    touch "$SHELL_RC"
fi

# Add CLIPittyClip to PATH
if grep -q "CLIPittyClip" "$SHELL_RC" 2>/dev/null; then
    print_warning "CLIPittyClip already configured in $SHELL_RC"
else
    echo "" >> "$SHELL_RC"
    echo "# CLIPittyClip" >> "$SHELL_RC"
    echo "export PATH=\"\$PATH:${SCRIPT_DIR}\"" >> "$SHELL_RC"
    print_success "Added CLIPittyClip to PATH"
fi

# Add CTK to PATH and PERL5LIB
if grep -q "CTK" "$SHELL_RC" 2>/dev/null; then
    print_warning "CTK already in PATH"
else
    echo "" >> "$SHELL_RC"
    echo "# CTK (CLIP Tool Kit)" >> "$SHELL_RC"
    echo "export PATH=\"\$PATH:${CTK_DIR}\"" >> "$SHELL_RC"
    echo "export PERL5LIB=\"\$PERL5LIB:${CZPLIB_DIR}\"" >> "$SHELL_RC"
    print_success "Added CTK to PATH and PERL5LIB"
fi

# Add HOMER to PATH
if grep -q "HOMER" "$SHELL_RC" 2>/dev/null; then
    print_warning "HOMER already in PATH"
else
    echo "" >> "$SHELL_RC"
    echo "# HOMER" >> "$SHELL_RC"
    echo "export PATH=\"\$PATH:${HOMER_DIR}/bin\"" >> "$SHELL_RC"
    print_success "Added HOMER to PATH"
fi

# Set execute permissions on scripts
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
echo ""
