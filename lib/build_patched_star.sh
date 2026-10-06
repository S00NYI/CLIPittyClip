#!/usr/bin/env bash
#===============================================================================
# build_patched_star.sh -- build STAR 2.7.11b for macOS with the libc++ fix
#
# Usage:
#   lib/build_patched_star.sh --dest DIR [--omp-prefix PREFIX] [--jobs N] [--keep-build]
#
#   --dest DIR         where the finished STAR binary is installed (e.g. $CONDA_PREFIX/bin)
#   --omp-prefix PFX   prefix holding lib/libomp.dylib and include/omp.h
#                      (default: $CONDA_PREFIX; the conda package is `llvm-openmp`)
#   --jobs N           parallel compile jobs (default: number of CPU cores)
#   --keep-build       keep the temporary build directory (for debugging)
#
# What it does: downloads the OFFICIAL STAR 2.7.11b source, verifies its SHA-256,
# applies the vendored patch (lib/patches/), compiles with Apple's clang, and
# installs the binary. About 30 seconds on a 10-core Mac. Nothing is installed
# system-wide; an existing STAR in --dest is kept as STAR.orig.
#
# Why: every macOS build of STAR 2.7.11b reads 0 reads (libc++ ignores
# pubsetbuf). Full explanation, evidence and removal criteria:
# lib/patches/README.md. Called by install_macos.sh only when
# lib/star_selftest.sh shows the conda STAR is broken.
#===============================================================================
set -euo pipefail

STAR_VERSION="2.7.11b"
STAR_URL="https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz"
STAR_SHA256="3f65305e4112bd154c7e22b333dcdaafc681f4a895048fa30fa7ae56cac408e7"   # same value Homebrew pins
PATCH_NAME="star-${STAR_VERSION}-libcxx-streambuf.patch"
PATCH_SHA256="ed4415b036d4ac9ba15ea67039cacd0c3c2293b531753f15875e831b508091c9"   # upstream PR #2691 diff, unmodified

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH_FILE="$HERE/patches/$PATCH_NAME"

DEST=""
OMP_PREFIX="${CONDA_PREFIX:-}"
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
KEEP_BUILD=false

die() { echo "build_patched_star: ERROR: $*" >&2; exit 1; }
say() { echo "build_patched_star: $*"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dest)       DEST="${2:-}"; shift 2 ;;
        --omp-prefix) OMP_PREFIX="${2:-}"; shift 2 ;;
        --jobs)       JOBS="${2:-}"; shift 2 ;;
        --keep-build) KEEP_BUILD=true; shift ;;
        -h|--help)    sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)            die "unknown option: $1" ;;
    esac
done

#--- preflight ------------------------------------------------------------------
[[ "$(uname)" == "Darwin" ]] || die "this is only needed on macOS (Linux STAR is unaffected)"
[[ -n "$DEST" ]] || die "--dest is required"
[[ -d "$DEST" ]] || die "--dest directory does not exist: $DEST"
[[ -n "$OMP_PREFIX" ]] || die "no OpenMP prefix: pass --omp-prefix or activate the conda env"
[[ -f "$OMP_PREFIX/lib/libomp.dylib" && -f "$OMP_PREFIX/include/omp.h" ]] \
    || die "OpenMP not found in $OMP_PREFIX (need lib/libomp.dylib and include/omp.h). Install: mamba install -n <env> llvm-openmp"
[[ -f "$PATCH_FILE" ]] || die "patch file missing: $PATCH_FILE"
for tool in /usr/bin/clang++ make xxd patch curl shasum tar; do
    command -v "$tool" > /dev/null 2>&1 || die "required tool not found: $tool (install Xcode Command Line Tools: xcode-select --install)"
done

sha256_of() { shasum -a 256 "$1" | awk '{print $1}'; }
[[ "$(sha256_of "$PATCH_FILE")" == "$PATCH_SHA256" ]] \
    || die "patch file checksum mismatch ($PATCH_FILE was modified; expected $PATCH_SHA256)"

BUILD="$(mktemp -d "${TMPDIR:-/tmp}/star_build.XXXXXX")"
cleanup() { if $KEEP_BUILD; then say "build directory kept: $BUILD"; else rm -rf "$BUILD"; fi; }
trap cleanup EXIT

#--- fetch + verify + patch -----------------------------------------------------
say "downloading official STAR $STAR_VERSION source..."
curl -fsSL "$STAR_URL" -o "$BUILD/star.tar.gz" || die "download failed: $STAR_URL"
[[ "$(sha256_of "$BUILD/star.tar.gz")" == "$STAR_SHA256" ]] \
    || die "source checksum mismatch (expected $STAR_SHA256). Refusing to build an unverified download."

tar xzf "$BUILD/star.tar.gz" -C "$BUILD"
SRC="$BUILD/STAR-$STAR_VERSION"
[[ -d "$SRC/source" ]] || die "unexpected source layout in the tarball"
(cd "$SRC" && patch -p1 --quiet < "$PATCH_FILE") || die "patch did not apply cleanly"
say "applied $PATCH_NAME"

#--- compiler shim: Apple clang spells OpenMP differently ------------------------
# STAR's Makefile passes a bare -fopenmp, which Apple clang rejects. The shim
# rewrites it to the supported spelling and points at a directory that holds ONLY
# omp.h (so none of conda's other headers leak into the build).
mkdir -p "$BUILD/omp_include"
ln -s "$OMP_PREFIX/include/omp.h" "$BUILD/omp_include/omp.h"
SHIM="$BUILD/clang-openmp"
cat > "$SHIM" <<EOF
#!/bin/bash
args=()
for a in "\$@"; do
    if [ "\$a" = "-fopenmp" ]; then args+=(-Xpreprocessor -fopenmp "-I$BUILD/omp_include"); else args+=("\$a"); fi
done
exec /usr/bin/clang++ "\${args[@]}"
EOF
chmod +x "$SHIM"

#--- compile ---------------------------------------------------------------------
# CXXFLAGS_SIMD=            STAR's default (-mavx2) is x86-only; its SIMD code uses the
#                           portable SIMDe layer, so no flag is needed.
# -DCOMPILE_FOR_MAC         STAR's own macOS switch.
# LDFLAGS_shared            Apple's linker has no -Bstatic/-Bdynamic; link STAR's bundled
#                           static htslib directly (a newer htslib on the path lacks a
#                           symbol STAR still uses).
# libomp by full path       so only OpenMP comes from the conda env; libc++ and zlib are
#                           the system ones.
say "compiling with $JOBS jobs (about 30 s on 10 cores; a few minutes on fewer)..."
if ! (cd "$SRC/source" && make -j"$JOBS" STAR \
        CXX="$SHIM" \
        CXXFLAGS_SIMD= \
        CXXFLAGSextra="-DCOMPILE_FOR_MAC" \
        LDFLAGS_shared="-pthread htslib/libhts.a -lz" \
        LDFLAGSextra="$OMP_PREFIX/lib/libomp.dylib -Wl,-rpath,$OMP_PREFIX/lib" \
        > "$BUILD/build.log" 2>&1); then
    echo "---- last lines of the build log ----" >&2
    grep -E "error|Error|ld:" "$BUILD/build.log" | cut -c1-200 | tail -15 >&2 || true
    KEEP_BUILD=true
    die "compilation failed (full log: $BUILD/build.log)"
fi

NEW="$SRC/source/STAR"
[[ -x "$NEW" ]] || die "build finished but no STAR binary was produced"
file "$NEW" | grep -q "$(uname -m)" || die "built binary is not native ($(file "$NEW" | cut -d: -f2))"

#--- install ---------------------------------------------------------------------
if [[ -f "$DEST/STAR" && ! -f "$DEST/STAR.orig" ]]; then
    cp -p "$DEST/STAR" "$DEST/STAR.orig"
    say "kept the previous STAR as $DEST/STAR.orig"
fi
cp "$NEW" "$DEST/STAR.tmp.$$" && chmod 755 "$DEST/STAR.tmp.$$" && mv -f "$DEST/STAR.tmp.$$" "$DEST/STAR"

cat > "$DEST/STAR.patch-info.txt" <<EOF
STAR $STAR_VERSION, built from source with a libc++ fix, by CLIPittyClip (lib/build_patched_star.sh)
built:        $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(uname -m), macOS $(sw_vers -productVersion 2>/dev/null || echo ?)
source:       $STAR_URL
source sha256 $STAR_SHA256
patch:        $PATCH_NAME (sha256 $PATCH_SHA256)
why:          stock STAR $STAR_VERSION reads 0 reads on macOS; see lib/patches/README.md
EOF

say "installed patched STAR -> $DEST/STAR"
