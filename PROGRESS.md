# PROGRESS — eqice Stage 1 implementation

Tracking progress on PLAN.md (Stage 1 — instrument the C++ codes). Last
updated: 2026-09-30.

## Milestone 1 — Confirm PISM source access + build; pick reference run configs ✅

- **Source access** ✅ — PISM v2.3.2 at `stage1/pism` (git tag `v2.3.2`,
  commit `fa11747`); PETSc v3.26.0 at `stage1/petsc` (in-place build,
  `arch-linux-c-opt`).
- **Build** ✅ — PETSc (self-contained MPI/HDF5/NetCDF via downloads; system
  flexiblas; openmpi 5.0.1 + gcc 13.3.0 from cluster modules) and PISM
  (Release, installed to `stage1/build/install`). Build scripts:
  `stage1/build/build-petsc.sh`, `stage1/build/build-pism.sh`,
  `stage1/build/build-local-deps.sh` (expat 2.6.4 + udunits2 2.2.28, built
  from source into `stage1/build/local` because the cluster has neither).
  Runtime environment: `stage1/runs/env.sh`.
- **Reference run configs** ✅ — `stage1/runs/README.md`: PISM verification
  tests A, B, C, D, F, G, H, K, L (smoke-tested ✅), SSA via
  `pism_ssa_test_*` executables, BTU via `pism_btutest`. Test V (full-model
  SSA) fails in v2.3.2 — see bugs.md #3.
- **Boundary catalog** ✅ — `stage1/boundaries.md`: the component inventory
  (PLAN.md §3) mapped to concrete PISM 2.3.2 classes/methods with signatures,
  inputs→outputs, and derivative-vs-integrated classification. This is the
  deliverable of Stage-1 step 1 and the target list for instrumentation.

## Stage 1 — Instrumentation (in progress)

- [x] Step 1: identify subassembly boundaries → `stage1/boundaries.md`
- [ ] Step 2: instrument each boundary to dump inputs/outputs
      (in progress — `stage1/instrument/instrument_rheology.cc` started;
      rheology is the first component per PLAN.md §10)
- [ ] Step 3: dump subassembly I/O (e.g. SSA linear solve vs. SIA stencil)
- [ ] Step 4: record run configurations → `stage1/runs/README.md` ✅ (configs);
      curated dumps pending step 2
- [ ] Step 5: extract select dumps into numeric test tuples

## Stage 2 — Stubs → physics → EarthSciModels PRs (not started)

Per PLAN.md §5: hand-author stub `.esm` files with `tests` blocks from
Stage-1 tuples, fill in physics, review + merge into EarthSciModels.

## Stage 3 — Top-level `eqice.esm` (not started)

## Key facts

- PISM: v2.3.2 (`stage1/pism`), built to `stage1/build/install`
- PETSc: v3.26.0 (`stage1/petsc`, `arch-linux-c-opt`)
- Toolchain: gcc 13.3.0 + openmpi 5.0.1 (cluster modules), GSL 2.8,
  FFTW 3.3.10 (module prefixes), expat 2.6.4 + udunits2 2.2.28 (local build)
- EarthSciAST `esm` CLI: built from origin/main (17c57b5bc) at
  `../EarthSciAST/pkg/earthsci-ast-rs/target/release/esm`
- Local PISM patch: `stage1/pism/CMake/FindUDUNITS2.cmake` (UDUNITS-2 optional
  when expat is absent — see bugs.md #1)
- Known issue: full-model `-test V` fails (bugs.md #3); use `pism_ssa_test_*`

## Repository layout (see PLAN.md §7)

```
eqice/
  AGENTS.md, PLAN.md, PROGRESS.md, bugs.md
  stage1/
    boundaries.md      # subassembly boundary catalog (step 1 deliverable)
    dumps/             # instrumented C++ I/O traces (pending)
    runs/              # run configurations + env.sh
    instrument/        # instrumentation drivers (in progress)
    build/             # build scripts + artifacts (gitignored)
    pism/, petsc/      # vendored sources (gitignored)
```
