#!/bin/bash
# Build PISM v2.3.2 against the Stage-1 PETSc build.
#
# Reproducible configuration:
#   - PISM v2.3.2 (stage1/pism, git tag v2.3.2)
#   - PETSc v3.26.0 built by build-petsc.sh (PETSC_DIR/PETSC_ARCH below)
#   - MPI: openmpi 5.0.1 (gcc 13.3.0) from the cluster module system
#   - NetCDF/HDF5: from PETSc's downloaded builds (via pkg-config)
#   - GSL 2.8, FFTW 3.3.10: cluster module prefixes (pkg-config)
#   - Optional deps (PROJ, UDUNITS2, PNETCDF) disabled
#   - Release build, 8 parallel jobs (shared node)
#
# Usage: bash stage1/build/build-pism.sh [jobs]

set -euo pipefail

JOBS="${1:-8}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD_DIR="$REPO_ROOT/stage1/build/pism"
INSTALL_PREFIX="$REPO_ROOT/stage1/build/install"

export PETSC_DIR="$REPO_ROOT/stage1/petsc"
export PETSC_ARCH="arch-linux-c-opt"

MPI_PREFIX="/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0"
GCC_PREFIX="/sw/apps/gcc/13.3.0"
FFTW_PREFIX="/sw/apps/fftw/3.3.10/gcc-13.3.0/openmp-5.0.1"
GSL_PREFIX="/sw/apps/gsl/2.8"

LOCAL_PREFIX="$REPO_ROOT/stage1/build/local"   # expat + udunits2 (build-local-deps.sh)
export UDUNITS2_ROOT="$LOCAL_PREFIX"

export PATH="$MPI_PREFIX/bin:$GCC_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_PREFIX/lib:$GCC_PREFIX/lib64:$PETSC_DIR/$PETSC_ARCH/lib:$LOCAL_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:$FFTW_PREFIX/lib/pkgconfig:$GSL_PREFIX/lib/pkgconfig:$LOCAL_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

if [[ ! -f CMakeCache.txt ]]; then
  echo "==> Configuring PISM..."
  cmake "$REPO_ROOT/stage1/pism" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DPism_USE_PROJ=OFF \
    -DPism_USE_UDUNITS=OFF \
    -DPism_USE_PNETCDF=OFF \
    -DPism_USE_PARALLEL_NETCDF4=OFF \
    -DPism_BUILD_PYTHON_BINDINGS=OFF \
    -DPism_USE_YAC=OFF \
    -DPism_BUILD_EXTRAS=OFF \
    -DPism_DEBUG=OFF
fi

echo "==> Building PISM (-j $JOBS)..."
make -j "$JOBS"

echo "==> Installing PISM to $INSTALL_PREFIX..."
make install -j "$JOBS"

echo "==> PISM build complete."
echo "    Build directory: $BUILD_DIR"
echo "    Install prefix:  $INSTALL_PREFIX"
