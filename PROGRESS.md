# PROGRESS — eqice Stage 1 implementation

Tracking progress on PLAN.md (Stage 1 — instrument the C++ codes). Last
updated: 2026-10-01.

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

## Milestone 3 (partial) — energy + age boundaries instrumented ✅

- **Energy traces** ✅ — new `instrument_energy.cc`: one `dt = 10 yr` step of
  `EnthalpyModel` forced entirely by the test-F exact solution (u3/v3/w3, Σ,
  T_s(r), Ggeo = 0.042). Initial state = exact steady-state temperature
  (converted with the model's own standard `EnthalpyConverter`), so the
  interior response is the discrete-vs-continuous gap (~2 J/kg/step); the top
  ice level is pulled ~1 K toward the surface Dirichlet BC; where the exact T
  exceeds T_pmp (deep dome interior) the model goes temperate and melts up to
  0.0037 m/yr (57 temperate cells). Dumps in `stage1/dumps/energy/`.
- **Age traces** ✅ — same driver: one `dt` step of `AgeModel` from age = 0
  with the exact test-F velocity field. Deep-interior age = dt exactly (max
  = dt), surface transition carried by the top coarse levels. Dumps in
  `stage1/dumps/age/`.
- **Run config** ✅ — `stage1/runs/energy_age.json` (grid, dt, forcing,
  check values).
- Note: the first *integrated* (time-advancing) traces; the driver exposes the
  protected enthalpy state via a subclass to set the exact-T initial condition.

## Milestone 3 (continued) — basal strength boundary instrumented ✅

- **Basal strength traces** ✅ — new `instrument_basal_strength.cc`:
  - pointwise `MohrCoulombPointwise` grids: `effective_pressure`
    `(delta, P_overburden, W_till) -> N_till` (with the un-clamped interior
    value), `yield_stress -> tauc = c0 + tan(phi) N_till`, and the
    `till_friction_angle` inverse (round trip `phi -> tauc -> phi` exact to
    7e-15 deg). Verified: `N_till(W=0) = P_overburden`, `N_till(W=Wmax) =
    delta*P`, monotone non-increasing in `W`.
  - basal resistance laws (`plastic`, `pseudo_plastic`, `regularized`):
    `drag(tauc, vx, vy)` + `drag_with_derivative` over zero/slow/fast speeds;
    `drag_with_derivative` agrees with finite differences to ~1e-10; plastic
    limit `|tau_b| -> tauc` at high speed.
  - a `MohrCoulombYieldStress` tauc map on the test-F exact geometry (flat
    bed, Mx=My=31) with smooth prescribed `W_till(r)` and `phi(r)` fields:
    ice-free `tauc = 1e6` (ice_free_bedrock), grounded range
    [1.19e4, 1.55e7] Pa, dome-center `tauc = tan(30°)·P = 1.546e7`.
- **Run config** ✅ — `stage1/runs/basal_strength.json` (sample grids, forcing
  choice, check values).

## Milestone 3 (continued) — hydrology boundary instrumented ✅

- **Hydrology traces** ✅ — new `instrument_hydrology.cc`: one `dt = 0.25 yr`
  step of `hydrology::Routing` (PISM's default `hydrology.model = "routing"`,
  Shreve `q = -K grad psi`) on the test-F exact geometry (flat bed, Mx=My=31).
  Prescribed forcing: constant surface input rate 0.2 m/yr water-equivalent,
  dome-centered basal melt blob 0.1 m/yr `exp(-(r/200km)^2)`, constant sliding
  speed 100 m/yr (part of the Hydrology Inputs contract but **ignored by
  Routing** — only `Distributed` uses it; documented). State-in: dry till
  `W_till = 0` and a dome-centered transportable-water blob
  `W = 0.5 m exp(-(r/150km)^2)`.
  - The model takes 19 internal substeps (CFL-limited to ~5 days by
    `max_timestep_W_cfl()`, dominated by the `eps = 1e-6` regularization at
    `V/dx ~ 5e-8`; `dt_diff ~ 21 yr` not binding).
  - Validation: water-mass conservation closes to round-off both in the model's
    own accounting (residual 0 kg) and in an independent state-based
    recomputation (~8e-15 relative); `W_till` non-negative (max 0.0725 m =
    3 months of input absorbed by the till); `overburden = rho_ice g H` exactly
    (max error 0 Pa over 517 grounded cells); advective flux satisfies
    `q . grad R <= 0` (max `q.gradR = 0`); flux centrally anti-symmetric to
    4e-16 relative; max |q| = 8.15e-4 m²/s.
  - Pointwise flux-law dump (`flux_law.csv`): `K = k W^(alpha-1) (G²+eps²)^((beta-2)/2)`
    (with the beta<2 regularization) and `|q| = K W G` over `(W, G)` grids — an
    analytic .esm test target. Staggered subassembly dump (`staggered.csv`):
    `Wstag, Kstag, Vstag, Qstag` from the last internal substep.
- **Run config** ✅ — `stage1/runs/hydrology.json` (grid, dt, forcing, checks).
- `Distributed`, `SteadyState`/`EmptyingProblem`, `NullTransport` hydrology
  models are **not** instrumented yet (deferred to a later pass; noted in
  `stage1/boundaries.md` §7 and the hydrology README section).

## Stage 1 — Remaining instrumentation

- [x] Step 1: identify subassembly boundaries → `stage1/boundaries.md`
- [x] Step 2: instrument boundaries to dump inputs/outputs — **rheology + SIA +
      SSA + energy + age + basal strength + hydrology (Routing) done**;
      remaining: Distributed/SteadyState hydrology, bed, surface/ocean/calving,
      geometry
- [x] Step 3 (partial): dump subassembly I/O — SSA `nuH`/`taud` (FD
      subassemblies), SIA `delta`/`D`/`q` chains, the
      `MohrCoulombPointwise`/basal-resistance-law pointwise functions, and the
      Routing staggered `Wstag/Kstag/Vstag/Qstag` substep fields are dumped;
      deeper subassemblies (e.g. the assembled KSP matrix) are deferred
- [x] Step 4: record run configurations → `stage1/runs/README.md` +
      `stage1/runs/{rheology,sia,ssa,energy_age,basal_strength,hydrology}.json`;
      dumps stored in `stage1/dumps/`
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
- `stage1/dumps/energy/` + `stage1/dumps/age/` — 5 + 4 files, ~32 MB total
  (regenerate: `instrument/instrument_energy -Mx 31 -My 31 -Mz 61
  -dumps_dir stage1/dumps -dt_years 10`)
- `stage1/dumps/basal_strength/` — 10 files, ~0.8 MB (regenerate:
  `instrument/instrument_basal_strength -Mx 31 -My 31 -dumps_dir
  stage1/dumps/basal_strength`)
- `stage1/dumps/hydrology/` — 6 files, ~1.2 MB (regenerate:
  `instrument/instrument_hydrology -Mx 31 -My 31 -dt_years 0.25 -dumps_dir
  stage1/dumps/hydrology`)
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
