#!/usr/bin/env bash
# Environment setup for building and running the instrumented PISM (Stage 1).
#
# Usage: source stage1/env.sh
#
# Records the exact toolchain used to build PETSc/PISM so runs are reproducible
# (see stage1/runs/ for run configurations).
set -euo pipefail

EQICE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Toolchain: GCC 13.3.0 + OpenMPI 5.0.1 (the OpenMPI module is built against
# GCC 13.3.0, so load both to keep the ABI consistent).
module purge
module load gcc/13.3.0
module load openmpi/5.0.1-gcc-13.3.0
module load gsl/2.8
module load fftw/3.3.10

export EQICE_ROOT
export PETSC_DIR="${EQICE_ROOT}/stage1/petsc"
export PETSC_ARCH="arch-linux-c-opt"
export PISM_SRC="${EQICE_ROOT}/stage1/pism"
export PISM_BUILD="${EQICE_ROOT}/stage1/pism-build"
export PISM_INSTALL="${EQICE_ROOT}/stage1/pism-install"
export UDUNITS_ROOT="${EQICE_ROOT}/stage1/udunits/install"
export LD_LIBRARY_PATH="${UDUNITS_ROOT}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${UDUNITS_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

echo "Environment ready:"
echo "  gcc:      $(gcc --version | head -1)"
echo "  mpicc:    $(mpicc --version | head -1)"
echo "  mpirun:   $(mpirun --version | head -1)"
echo "  PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH}"
