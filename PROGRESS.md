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

## Milestone 2 — Stage-1 instrumentation scaffolding + first component traces (rheology, SIA) ✅

- **Instrumentation scaffolding** ✅ — `stage1/build/build_instrument.sh`
  builds all `stage1/instrument/instrument_*.cc` drivers against the installed
  PISM library (fixed to use the working env + correct pkg-config paths).
- **Rheology traces** ✅ — `instrument_rheology.cc` (completed; fixes: correct
  `secondInvariant_2D` namespace, typed config dump with units) → 36 dump files
  in `stage1/dumps/rheology/`: softness/hardness/flow (+vectorized `flow_n`,
  +effective viscosity, +averaged hardness) for all 7 flow laws, enthalpy
  converter, second invariant, parameters. Validated: Paterson-Budd
  A(223.15 K) = 3.26e-27 Pa⁻³ s⁻¹, isothermal-Glen B̄ = 6.81e7 Pa s^(1/n).
- **SIA traces** ✅ — new `instrument_sia.cc`: one full SIA update from the
  test-F exact-solution state (flat bed, `arr` flow law, no sliding) → 11 dump
  files in `stage1/dumps/sia/`, including the complete
  `(alpha, pressure, E, stress, flow, delta) -> D` chain (verified: trapezoid of
  delta = D to full precision) and per-column computed-vs-exact u3/v3/w3/Σ.
- **Run configs recorded** ✅ — `stage1/runs/rheology.json`, `stage1/runs/sia.json`
  (exact commands, grids, sample distributions, and check values).
- **New C++ bug found** — `StressBalance::Inputs::dump()` writes an undefined
  NetCDF variable (bugs.md #4); worked around in the SIA driver.

## Stage 1 — Remaining instrumentation

- [x] Step 1: identify subassembly boundaries → `stage1/boundaries.md`
- [x] Step 2: instrument boundaries to dump inputs/outputs — **rheology + SIA +
      SSA done** (the three stress-balance boundaries; next: energy/age, basal,
      hydrology, bed, surface/ocean/calving, geometry)
- [x] Step 3 (partial): dump subassembly I/O — SSA `nuH`/`taud` (FD
      subassemblies) and SIA `delta`/`D`/`q` chains are dumped; deeper
      subassemblies (e.g. the assembled KSP matrix) are deferred
- [x] Step 4: record run configurations → `stage1/runs/README.md` +
      `stage1/runs/{rheology,sia,ssa}.json`; dumps stored in `stage1/dumps/`
- [ ] Step 5: extract select dumps into numeric test tuples (feeds Stage 2)

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
- Known issue: `StressBalance::Inputs::dump()` broken (bugs.md #4); drivers
  write NetCDF directly

## Dumps

- `stage1/dumps/rheology/` — 36 files, ~3 MB (regenerate:
  `instrument/instrument_rheology`)
- `stage1/dumps/sia/` — 11 files, ~28 MB (regenerate:
  `instrument/instrument_sia -Mx 31 -My 31 -Mz 61`)
- `stage1/dumps/ssa/` — 4 files, ~0.3 MB (regenerate:
  `instrument/instrument_ssa -Mx 61 -My 3`)
- Full dumps are gitignored (`stage1/dumps/*`); `stage1/dumps/README.md`
  documents the layout and regeneration. Numeric test tuples extracted from
  these go into Stage-2 `.esm` `tests` blocks.

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
