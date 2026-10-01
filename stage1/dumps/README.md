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

Both drivers accept `-dumps_dir <path>` to redirect the output.

## Regenerating

```bash
source stage1/runs/env.sh
bash stage1/build/build_instrument.sh          # builds both drivers
stage1/instrument/instrument_rheology -dumps_dir stage1/dumps/rheology
stage1/instrument/instrument_sia -Mx 31 -My 31 -Mz 61 -dumps_dir stage1/dumps/sia
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
- The test-F state means every dumped quantity has an exact reference value
  (see `sia/exact_solution.csv`); the discrete-vs-exact differences (surface
  velocity ~0.2 m/yr at 31x31) are the expected discretization error of the
  C++ model, which the `.esm` reimplementation must reproduce.

## Known C++ issue hit while dumping

`StressBalance::Inputs::dump()` (the natural input-dump hook) is broken in
v2.3.2: it calls `io::write_config()` on a `pism_config` variable that was
never defined, aborting with "NetCDF: Variable not found". Recorded in
`bugs.md` #4. The SIA driver writes `inputs.nc` directly instead.
