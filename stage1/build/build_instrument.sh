#!/usr/bin/env bash
# Build the Stage-1 instrumentation drivers against the installed PISM library.
#
# Requires a completed PISM install (see build-pism.sh -> stage1/build/install).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../runs/env.sh"   # working toolchain + PETSc pkg-config

INSTALL="${EQICE_ROOT:-$REPO_ROOT}/stage1/build/install"
LOCAL="${EQICE_ROOT:-$REPO_ROOT}/stage1/build/local"
export PKG_CONFIG_PATH="${INSTALL}/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
LDIR="-L${LOCAL}/lib"

for driver in "${EQICE_ROOT:-$REPO_ROOT}"/stage1/instrument/instrument_*.cc; do
  [ -e "$driver" ] || continue
  bin="${driver%.cc}"
  if [ -n "$(command -v mpicxx)" ]; then
    CXX=mpicxx
  else
    CXX="${CXX:-g++}"
  fi
  echo "Compiling ${driver} -> ${bin}"
  ${CXX} -std=c++17 -O2 \
    $(pkg-config --cflags pism) $(pkg-config --cflags petsc) \
    "${driver}" \
    -o "${bin}" \
    $(pkg-config --libs pism) \
    $(pkg-config --libs petsc) \
    $(pkg-config --libs gsl) \
    $(pkg-config --libs netcdf) \
    $(pkg-config --libs fftw3) \
    ${LDIR} -ludunits2
done
echo "All instrumentation drivers built."
