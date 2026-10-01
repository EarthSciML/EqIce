# Shared runtime environment for PISM Stage-1 runs.
# Source this from run scripts:  source stage1/runs/env.sh
#
# Provides:
#   PISM_BIN   - path to the PISM build directory
#   PATH       - MPI + gcc toolchain
#   LD_LIBRARY_PATH - PETSc, HDF5/NetCDF (from PETSc), FFTW, GSL, expat, udunits2

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

export PETSC_DIR="$REPO_ROOT/stage1/petsc"
export PETSC_ARCH="arch-linux-c-opt"

MPI_PREFIX="/sw/apps/mpi/openmpi/5.0.1/gcc/13.3.0"
GCC_PREFIX="/sw/apps/gcc/13.3.0"
FFTW_PREFIX="/sw/apps/fftw/3.3.10/gcc-13.3.0/openmp-5.0.1"
GSL_PREFIX="/sw/apps/gsl/2.8"
LOCAL_PREFIX="$REPO_ROOT/stage1/build/local"

export PATH="$MPI_PREFIX/bin:$GCC_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_PREFIX/lib:$GCC_PREFIX/lib64:$PETSC_DIR/$PETSC_ARCH/lib:$LOCAL_PREFIX/lib:$FFTW_PREFIX/lib:$GSL_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:$FFTW_PREFIX/lib/pkgconfig:$GSL_PREFIX/lib/pkgconfig:$LOCAL_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

export PISM_BIN="$REPO_ROOT/stage1/build/install/bin"
export PISM_PREFIX="$REPO_ROOT/stage1/build/install"
export LD_LIBRARY_PATH="$PISM_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PISM_PREFIX/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
