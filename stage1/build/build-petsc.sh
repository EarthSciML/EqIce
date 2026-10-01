#!/bin/bash
# Build PETSc (self-contained MPI/HDF5/NetCDF) for the eqice Stage-1 work.
#
# Reproducible configuration:
#   - PETSc v3.26.0 (stage1/petsc, git tag v3.26.0)
#   - MPI: openmpi 5.0.1 (gcc 13.3.0) from the cluster module system
#   - BLAS/LAPACK: system flexiblas (/usr/lib64/libflexiblas.so)
#   - NetCDF + HDF5: downloaded and built by PETSc configure
#   - Debugging off, -O3
#
# Usage: bash stage1/build/build-petsc.sh [jobs]
#   jobs: number of parallel make jobs (default 8; keep modest -- shared node)

set -euo pipefail

JOBS="${1:-8}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PETSC_SRC="$REPO_ROOT/stage1/petsc"
export PETSC_DIR="$PETSC_SRC"
export PETSC_ARCH="arch-linux-c-opt"

MPI_PREFIX="/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0"
GCC_PREFIX="/sw/apps/gcc/13.3.0"
FFTW_PREFIX="/sw/apps/fftw/3.3.10/gcc-13.3.0/openmp-5.0.1"
GSL_PREFIX="/sw/apps/gsl/2.8"

export PATH="$MPI_PREFIX/bin:$GCC_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_PREFIX/lib:$GCC_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$FFTW_PREFIX/lib/pkgconfig:$GSL_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

cd "$PETSC_SRC"

if [[ ! -f "$PETSC_ARCH/lib/petsc/conf/petscvariables" ]]; then
  echo "==> Configuring PETSc ($PETSC_ARCH)..."
  ./configure \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --with-blaslapack-lib="[/usr/lib64/libflexiblas.so]" \
    --with-zlib-dir=/usr \
    --download-hdf5 \
    --download-netcdf \
    --with-debugging=0 \
    --with-shared-libraries=1 \
    COPTFLAGS="-O3" \
    CXXOPTFLAGS="-O3" \
    FOPTFLAGS="-O3"
fi

echo "==> Building PETSc (-j $JOBS)..."
MAKE_NP="$JOBS" make -j "$JOBS" all

echo "==> PETSc build complete."
echo "    PETSC_DIR=$PETSC_DIR"
echo "    PETSC_ARCH=$PETSC_ARCH"
