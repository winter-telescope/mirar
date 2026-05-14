#!/usr/bin/env bash
# install_mirar.sh
# Creates a conda environment, installs mirar Python dependencies via poetry,
# and installs astrometry.net, SWarp, and PSFex for Mac M1 (arm64).
#
# Usage: bash install_mirar.sh [conda-env-name]
#        Defaults to env name "mirar".
#        Must be run from the root of the mirar repo.

set -euo pipefail

ENV_NAME="${1:-mirar}"

# ── Load dependency versions ──────────────────────────────────────────────────

if [[ ! -f astro-deps.cfg ]]; then
    echo "ERROR: astro-deps.cfg not found. Run this script from the root of the mirar repo."
    exit 1
fi

source astro-deps.cfg

echo "==> Dependency versions:"
echo "    scamp:             $SCAMP_VERSION"
echo "    source-extractor:  $SOURCE_EXTRACTOR_VERSION"
echo "    swarp:             $SWARP_VERSION"
echo "    psfex:             $PSFEX_VERSION"

# ── Check architecture ────────────────────────────────────────────────────────

ARCH=$(uname -m)
if [[ "$ARCH" != "arm64" ]]; then
    echo "WARNING: This script is intended for Apple Silicon (arm64). Detected: $ARCH"
    read -p "Continue anyway? [y/N] " -n 1 -r; echo
    [[ $REPLY =~ ^[Yy]$ ]] || exit 1
fi

# ── Check prerequisites ───────────────────────────────────────────────────────

if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniforge from:"
    echo "       https://github.com/conda-forge/miniforge/releases/latest"
    exit 1
fi

if ! command -v brew &>/dev/null; then
    echo "ERROR: Homebrew not found. Install it from https://brew.sh and re-run."
    exit 1
fi

if ! command -v poetry &>/dev/null; then
    echo "==> poetry not found, installing via pipx..."
    pip install --user pipx
    pipx install poetry
fi

if [[ ! -f pyproject.toml ]]; then
    echo "ERROR: pyproject.toml not found. Run this script from the root of the mirar repo."
    exit 1
fi

# Make sure conda shell functions are available
source "$(conda info --base)/etc/profile.d/conda.sh"

# ── Detect Python version from pyproject.toml ────────────────────────────────

PYTHON_VERSION=$(
    grep -E '^\s*requires-python\s*=' pyproject.toml \
    | grep -oE '[0-9]+\.[0-9]+' \
    | head -1
)

if [[ -z "$PYTHON_VERSION" ]]; then
    echo "ERROR: Could not detect Python version from pyproject.toml."
    exit 1
fi

echo ""
echo "==> Detected Python version: $PYTHON_VERSION"

# ── Create conda environment ──────────────────────────────────────────────────

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo ""
    echo "==> Removing existing conda environment: $ENV_NAME"
    conda env remove -y -n "$ENV_NAME"
fi

echo ""
echo "==> Creating conda environment: $ENV_NAME (python=$PYTHON_VERSION)"
conda create -y -n "$ENV_NAME" python="$PYTHON_VERSION"
conda activate "$ENV_NAME"

PREFIX="$CONDA_PREFIX"
echo "==> Environment prefix: $PREFIX"

# ── Install conda packages ────────────────────────────────────────────────────

echo ""
echo "==> Installing conda packages (native arm64)..."
conda install -y -c conda-forge \
    "astromatic-source-extractor=$SOURCE_EXTRACTOR_VERSION" \
    gsl \
    cfitsio fftw openblas \
    autoconf automake libtool pkg-config wcstools

echo "    source-extractor: $(sex -v 2>&1 | head -1)"

# ── astrometry.net via Homebrew ───────────────────────────────────────────────

echo ""
echo "==> Installing astrometry.net via Homebrew..."
brew install astrometry-net
echo "    $(solve-field --version 2>&1 | head -1)"

# ── Helper: build an astromatic tool from source ──────────────────────────────

build_astromatic() {
    local NAME="$1"
    local REPO="$2"
    local VERSION="$3"
    shift 3
    local EXTRA_CONFIGURE=("$@")

    echo ""
    echo "==> Building $NAME $VERSION from source..."

    local SRCDIR="$BUILDDIR/$NAME"
    git clone --depth 1 --branch "$VERSION" \
        "https://github.com/astromatic/$REPO.git" "$SRCDIR"
    pushd "$SRCDIR" > /dev/null

    sh autogen.sh
    ./configure \
        --prefix="$PREFIX" \
        "${EXTRA_CONFIGURE[@]}" \
        LDFLAGS="-Wl,-rpath,$PREFIX/lib"

    make -j"$(sysctl -n hw.logicalcpu)"
    make install
    popd > /dev/null

    # Patch rpath if configure didn't embed it
    local BIN
    BIN=$(which "$NAME" 2>/dev/null || true)
    if [[ -n "$BIN" ]] && ! otool -l "$BIN" | grep -q LC_RPATH; then
        echo "    Patching rpath in $NAME binary..."
        install_name_tool -add_rpath "$PREFIX/lib" "$BIN"
    fi
}

BUILDDIR=$(mktemp -d)
trap 'rm -rf "$BUILDDIR"' EXIT

# ── SCAMP ─────────────────────────────────────────────────────────────────────
# Scamp tag is vX.Y.Z on GitHub, but the binary reports just X.Y.Z, for... reasons
build_astromatic "scamp" "scamp" "v$SCAMP_VERSION" \
    --enable-openblas \
    --with-openblas-libdir="$PREFIX/lib" \
    --with-openblas-incdir="$PREFIX/include" \
    --with-fftw-incdir="$PREFIX/include" \
    --with-fftw-libdir="$PREFIX/lib"

echo "    scamp:            $(scamp -v 2>&1 | head -1)"

# ── SWarp ─────────────────────────────────────────────────────────────────────
# Swarp tag is just X.Y.Z without "v", for... reasons
build_astromatic "swarp" "swarp" "$SWARP_VERSION" \
    --with-cfitsio-incdir="$PREFIX/include" \
    --with-cfitsio-libdir="$PREFIX/lib"

echo "    $(swarp -v 2>&1 | head -1)"

# ── PSFex ─────────────────────────────────────────────────────────────────────
# PSFex tag is just X.Y.Z without "v", for... reasons
build_astromatic "psfex" "psfex" "$PSFEX_VERSION" \
    --enable-openblas \
    --with-openblas-libdir="$PREFIX/lib" \
    --with-openblas-incdir="$PREFIX/include" \
    --with-fftw-incdir="$PREFIX/include" \
    --with-fftw-libdir="$PREFIX/lib"

echo "    $(psfex -v 2>&1 | head -1)"

# ── Install mirar Python dependencies via poetry ──────────────────────────────

echo ""
echo "==> Installing mirar Python dependencies via poetry..."
poetry config virtualenvs.create false --local
poetry install -E dev

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "✅ mirar environment ready!"
echo ""
echo "   conda env:   $ENV_NAME"
echo "   python:      $(python --version)"
echo ""
echo "   solve-field  $(solve-field --version 2>&1 | head -1)"
echo "   swarp        $(swarp -v 2>&1 | head -1)"
echo "   psfex        $(psfex -v 2>&1 | head -1)"
echo "   sex          $(sex -v 2>&1 | head -1)"
echo "   scamp        $(scamp -v 2>&1 | head -1)"
echo ""
echo "   To activate: conda activate $ENV_NAME"
