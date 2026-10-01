#!/bin/bash
# Build expat + udunits2 from source into stage1/build/local.
#
# PISM's bundled "calcalcs" calendar library (src/external/calcalcs) is a
# hard dependency that requires udunits2.h and libudunits2, which are not
# available as system packages on this cluster. udunits2 itself requires
# expat. Both are small; we build them from source into a local prefix.
#
# Usage: bash stage1/build/build-local-deps.sh [jobs]

set -euo pipefail

JOBS="${1:-4}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PREFIX="$REPO_ROOT/stage1/build/local"
SRC_DIR="$REPO_ROOT/stage1/build/local-src"
mkdir -p "$SRC_DIR"

export PATH="$PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:${LD_LIBRARY_PATH:-}"

# --- expat ---
EXPAT_VER="2.6.4"
EXPAT_TGZ="expat-$EXPAT_VER.tar.gz"
if [[ ! -f "$PREFIX/include/expat.h" ]]; then
  echo "==> Building expat $EXPAT_VER..."
  cd "$SRC_DIR"
  if [[ ! -d "expat-$EXPAT_VER" ]]; then
    curl -sL -o "$EXPAT_TGZ" \
      "https://github.com/libexpat/libexpat/releases/download/R_2_6_4/$EXPAT_TGZ"
    tar xzf "$EXPAT_TGZ"
  fi
  cd "expat-$EXPAT_VER"
  ./configure --prefix="$PREFIX" >/dev/null
  make -j "$JOBS" >/dev/null
  make install >/dev/null
fi

# --- udunits2 ---
UDUNITS_VER="2.2.28"
UDUNITS_TGZ="udunits-$UDUNITS_VER.tar.gz"
if [[ ! -f "$PREFIX/include/udunits2.h" ]]; then
  echo "==> Building udunits2 $UDUNITS_VER..."
  cd "$SRC_DIR"
  if [[ ! -d "udunits-$UDUNITS_VER" ]]; then
    curl -sL -o "$UDUNITS_TGZ" \
      "https://downloads.unidata.ucar.edu/udunits/$UDUNITS_VER/$UDUNITS_TGZ"
    tar xzf "$UDUNITS_TGZ"
  fi
  cd "udunits-$UDUNITS_VER"
  ./configure --prefix="$PREFIX" \
    CPPFLAGS="-I$PREFIX/include" \
    LDFLAGS="-L$PREFIX/lib" >/dev/null
  make -j "$JOBS" >/dev/null
  make install >/dev/null
fi

echo "==> Local deps complete:"
echo "    $PREFIX/include/{expat.h,udunits2.h}"
echo "    $PREFIX/lib/libudunits2.a"
