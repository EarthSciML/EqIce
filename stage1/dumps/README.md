# Stage 1 — Instrumented C++ I/O traces

Reference traces produced by the instrumentation drivers in `stage1/instrument/`.
These are the authoritative inputs/outputs of the instrumented PISM components
that Stage-2 `.esm` tests are built from (PLAN.md §4 step 5: "extract a select
subset of these dumps into numeric test tuples").

The full dump directories are **not committed** (see `.gitignore`); they are
regenerated with one command each and kept on disk for the medium horizon. The
numeric tuples that Stage-2 `.esm` tests actually use are committed inside the
`.esm` `tests` blocks.

## Layout

| Directory | Component | Driver | Run config | Regenerate |
|---|---|---|---|---|
| `rheology/` | Flow laws, enthalpy converter, effective viscosity | `instrument/instrument_rheology.cc` | `runs/rheology.json` | `instrument/instrument_rheology` |
| `sia/` | SIA stress balance (test-F exact state) | `instrument/instrument_sia.cc` | `runs/sia.json` | `instrument/instrument_sia -Mx 31 -My 31 -Mz 61` |
| `ssa/` | SSA stress balance (test V / van der Veen shelf) | `instrument/instrument_ssa.cc` | `runs/ssa.json` | `instrument/instrument_ssa -Mx 61 -My 3` |
| `energy/` | Ice energy (EnthalpyModel, one dt step, test-F forcing) | `instrument/instrument_energy.cc` | `runs/energy_age.json` | `instrument/instrument_energy -Mx 31 -My 31 -Mz 61` |
| `age/` | Ice age (AgeModel, one dt step, test-F velocity) | `instrument/instrument_energy.cc` | `runs/energy_age.json` | (same driver; see `energy/`) |

The energy and age traces come from the same driver (they share the test-F
state). All drivers accept `-dumps_dir <path>` to redirect the output; the
energy/age driver writes `energy/` and `age/` under it and takes `-dt_years`.

## Regenerating

```bash
source stage1/runs/env.sh
bash stage1/build/build_instrument.sh          # builds all drivers
stage1/instrument/instrument_rheology -dumps_dir stage1/dumps/rheology
stage1/instrument/instrument_sia -Mx 31 -My 31 -Mz 61 -dumps_dir stage1/dumps/sia
stage1/instrument/instrument_ssa -Mx 61 -My 3 -dumps_dir stage1/dumps/ssa
stage1/instrument/instrument_energy -Mx 31 -My 31 -Mz 61 -dumps_dir stage1/dumps -dt_years 10
```

Each driver writes `meta.txt` recording the PISM revision, config file, and
command line, plus `parameters.csv` with every config key that affects the
component (names, values, and units) — the constants a Stage-2 `.esm`
component must reproduce.

## Using the dumps for tests

- **Rheology** (`rheology/flowlaw_<law>.csv`): `(stress, E, p, gs) -> (softness,
  hardness, flow)` tuples for all seven flow laws; the `_batched` files exercise
  the vectorized `flow_n` path used by SIA/SSA; `effective_viscosity_<law>.csv`
  and `averaged_hardness_<law>.csv` cover the viscosity and column-integral
  boundaries.
- **SIA** (`sia/columns.csv`): per column `(geometry, E, u3, v3, w3, sigma3)`
  plus the exact test-F reference. `sia/delta.csv` contains the full
  `(alpha, pressure, E, stress, flow, delta)` chain that produces the
  diffusivity `D` (verified: the trapezoidal integral of `delta` equals `D` to
  full precision). `sia/inputs.nc` is the complete input+output NetCDF.
- **SSA** (`ssa/columns.csv`): per cell `(geometry, tauc, u, v, F_b, taud_x/y,
  nuH, Bbar)` plus the exact van der Veen reference `(H, u, v)`. The solve
  matches the exact solution to 0.077% average error at Mx=61 (Bbar = 1.9e8
  exactly, constant-flux relation u*H = V0*H0 verified). `ssa/inputs.nc` has
  all inputs + outputs including the FD subassemblies `nuH` and `taud`.
- **Energy** (`energy/columns.csv`): per `(i,j,k)` the full forcing
  `(surface_temp, u3, v3, w3, sigma)` and `(E_in, E_out)` for one dt step,
  plus the exact reference temperature. The initial state is the exact
  test-F steady-state temperature, so the interior response is tiny
  (~2 J/kg/step); the top ice level is pulled ~1 K toward the surface
  Dirichlet BC, and where the exact T exceeds the pressure-melting
  temperature (deep dome interior) the model goes temperate and melts up to
  0.0037 m/yr (`energy/basal_melt_rate.csv`).
- **Age** (`age/columns.csv`): per `(i,j,k)` `(u3, v3, w3, age_in, age_out)`
  for one dt step from age = 0. Deep-interior age = dt exactly (max = dt);
  the top coarse levels carry the age-0 surface transition advected down by
  w*dt.
- The test-F/V states mean every dumped quantity has an exact reference value
  (or a clean null state, as with age = 0); the discrete-vs-exact differences
  are the expected discretization error of the C++ model, which the `.esm`
  reimplementation must reproduce.

## Known C++ issue hit while dumping

`StressBalance::Inputs::dump()` (the natural input-dump hook) is broken in
v2.3.2: it calls `io::write_config()` on a `pism_config` variable that was
never defined, aborting with "NetCDF: Variable not found". Recorded in
`bugs.md` #4. The SIA driver writes `inputs.nc` directly instead.
