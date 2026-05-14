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

# ── Check prerequisites ───────────────────────────────────────────────────────

if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniforge from:"
    echo "       https://github.com/conda-forge/miniforge/releases/latest"
    exit 1
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

# After PREFIX is set, near the top of the build section:
NCORES=$(python -c "import os; print(os.cpu_count() or 1)")

# -- Work out what versions of packages to install based on pyproject.toml dependencies --

get_deps() {
python - <<'EOF'
import tomllib, re

with open("pyproject.toml", "rb") as f:
    deps = tomllib.load(f)["project"]["dependencies"]

# packages you want from conda-forge
conda_set = {
    "numpy", "scipy", "astropy", "pandas",
    "pyarrow", "matplotlib"
}

conda = []
pip = []

for d in deps:
    if not isinstance(d, str):
        continue

    name = re.match(r"^[a-zA-Z0-9_.\-]+", d).group(0)

    if name in conda_set:
        conda.append(d)
    else:
        pip.append(d)

# output as two shell-eval friendly lines
print("CONDA_DEPS=\"" + " ".join(conda) + "\"")
print("PIP_DEPS=\"" + " ".join(pip) + "\"")
EOF
}

# Get deps split between pip and conda
eval "$(get_deps)"

echo "==> Conda deps:"
echo "$CONDA_DEPS"

echo "==> Pip deps:"
echo "$PIP_DEPS"

# ── Install conda packages ────────────────────────────────────────────────────

echo ""
echo "==> Installing conda packages (native arm64)..."
conda install -y -c conda-forge \
    "astromatic-source-extractor=$SOURCE_EXTRACTOR_VERSION" \
    gsl \
    cfitsio fftw openblas curl \
    autoconf automake libtool pkg-config wcstools \
    libjpeg-turbo cairo swig netpbm zlib freetype expat "libblas=*=*openblas" \
    "liblapack=*=*openblas" \
    $CONDA_DEPS

echo "    source-extractor: $(sex -v 2>&1 | head -1)"

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

    make -j"$NCORES"
    make install
    popd > /dev/null

    # Patch rpath if configure didn't embed it
    local BIN
    BIN=$(which "$NAME" 2>/dev/null || true)
    if [[ -n "$BIN" ]] && command -v otool &>/dev/null && ! otool -l "$BIN" | grep -q LC_RPATH; then
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
    --with-fftw-libdir="$PREFIX/lib" \
    --with-curl-incdir="$PREFIX/include" \
    --with-curl-libdir="$PREFIX/lib"

echo "    scamp:            $(scamp -v 2>&1 | head -1)"
#
## ── SWarp ─────────────────────────────────────────────────────────────────────
## Swarp tag is just X.Y.Z without "v", for... reasons
build_astromatic "swarp" "swarp" "$SWARP_VERSION" \
    --with-cfitsio-incdir="$PREFIX/include" \
    --with-cfitsio-libdir="$PREFIX/lib"

echo "    $(swarp -v 2>&1 | head -1)"
#
## ── PSFex ─────────────────────────────────────────────────────────────────────
## PSFex tag is just X.Y.Z without "v", for... reasons
build_astromatic "psfex" "psfex" "$PSFEX_VERSION" \
    --enable-openblas \
    --with-openblas-libdir="$PREFIX/lib" \
    --with-openblas-incdir="$PREFIX/include" \
    --with-fftw-incdir="$PREFIX/include" \
    --with-fftw-libdir="$PREFIX/lib" \
    --with-curl-incdir="$PREFIX/include" \
    --with-curl-libdir="$PREFIX/lib"

echo "    $(psfex -v 2>&1 | head -1)"

# ── Install astrometry.net via clone and make  ─────────────────────────

_install_astrometry_activation() {
    local PREFIX="$1"

    mkdir -p "$PREFIX/etc/conda/activate.d"
    mkdir -p "$PREFIX/etc/conda/deactivate.d"

    cat > "$PREFIX/etc/conda/activate.d/astrometry.sh" <<EOF
export PATH="\$CONDA_PREFIX/astrometry/bin:\$PATH"
export ASTROMETRY_DATA="\$CONDA_PREFIX/astrometry/data"
export DYLD_LIBRARY_PATH="\$CONDA_PREFIX/lib:\$DYLD_LIBRARY_PATH"
EOF

    cat > "$PREFIX/etc/conda/deactivate.d/astrometry.sh" <<EOF
unset ASTROMETRY_DATA
EOF
}

build_astrometry_net() {
    local NAME="astrometry.net"
    local REPO="astrometry.net"
    local VERSION="$1"
    local PREFIX="$CONDA_PREFIX"
    local SRCDIR="$BUILDDIR/$NAME"

    echo "==> Building $NAME $VERSION"

    rm -rf "$SRCDIR"
    git clone --depth 1 --branch "$VERSION" \
        "https://github.com/dstndstn/$REPO.git" \
        "$SRCDIR"

    pushd "$SRCDIR" > /dev/null

    # ---- environment for compilation ----
    export LDFLAGS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib"
    export CPPFLAGS="-I$PREFIX/include"
    export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/share/pkgconfig"

    # cores (portable mac + linux)
    local NCORES
    NCORES=$(python -c "import os; print(os.cpu_count() or 1)")

    # ---- build ----
    make -j"$NCORES"

    # IMPORTANT: astrometry uses INSTALL_DIR, not just prefix
    make install INSTALL_DIR="$PREFIX/astrometry"

    popd > /dev/null

    _install_astrometry_activation "$PREFIX"
}

build_astrometry_net "$ASTROMETRY_NET_VERSION"

# ── Install mirar Python dependencies via pip ──────────────────────────────

echo "==> Installing mirar Python dependencies via pip..."
conda run -n "$ENV_NAME" pip install $PIP_DEPS --upgrade-strategy only-if-needed
conda run -n "$ENV_NAME" pip install -e ".[dev]" --no-deps

echo "==> Installing pre-commit hooks..."
git config --global --add safe.directory "$(pwd)"
conda activate "$ENV_NAME"
pre-commit install

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
