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
| `basal_strength/` | Basal yield stress (Mohr-Coulomb) + basal resistance laws | `instrument/instrument_basal_strength.cc` | `runs/basal_strength.json` | `instrument/instrument_basal_strength -Mx 31 -My 31` |
| `hydrology/` | Subglacial hydrology (Routing/Shreve, one dt step, test-F geometry) | `instrument/instrument_hydrology.cc` | `runs/hydrology.json` | `instrument/instrument_hydrology -Mx 31 -My 31 -dt_years 0.25` |

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
stage1/instrument/instrument_basal_strength -Mx 31 -My 31 -dumps_dir stage1/dumps/basal_strength
stage1/instrument/instrument_hydrology -Mx 31 -My 31 -dt_years 0.25 -dumps_dir stage1/dumps/hydrology
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
- **Basal strength** (`basal_strength/`): the pointwise Mohr-Coulomb functions
  are pure `.esm` test targets:
  - `effective_pressure.csv`: `(delta, P_overburden, water_thickness) ->
    N_till = min(P, N0 (delta P/N0)^s 10^((e0/Cc)(1-s)))` (with the un-clamped
    interior value dumped too). At `W=0`, `N_till = P_overburden`; at
    `W = W_till_max`, `N_till = delta*P`; monotone non-increasing in `W`.
  - `yield_stress.csv`: `tauc = c0 + tan(phi) N_till`; `till_friction_angle.csv`
    is the exact inverse (round trip `phi -> tauc -> phi` to 1e-15 deg).
  - `drag_{plastic,pseudo_plastic,regularized}.csv`: `(tauc, vx, vy) -> beta`
    with `tau_b = -beta v`, `dbeta` (w.r.t. `alpha = 0.5|v|^2`), and a
    finite-difference `dbeta_fd` that agrees to ~1e-10. Plastic limit
    `|tau_b| -> tauc` at high speed verified.
  - `inputs.nc` + `tauc_map.csv`: the full `MohrCoulombYieldStress::update`
    on the test-F exact geometry (flat bed, thickness = exact H, Mx=My=31)
    with smooth prescribed `W_till(r)` and `phi(r)` fields — the real
    `update_impl` path including the ice-free `tauc = 1e6` branch. The
    effective pressure field `N_till` is recomputed pointwise from the same
    inputs (the component does not store it).
- **Hydrology** (`hydrology/`): one `dt = 0.25 yr` step of `hydrology::Routing`
  (PISM's default hydrology model, Shreve `q = -K grad psi`) on the test-F
  exact geometry, forced by a constant surface input rate (0.2 m/yr
  water-equivalent), a dome-centered basal melt blob (0.1 m/yr peak), and a
  constant sliding speed (ignored by Routing). State-in: dry till `W_till = 0`
  and a dome-centered transportable-water blob `W = 0.5 m * exp(-(r/150km)^2)`.
  - `hydrology/columns.csv`: per cell `(geometry, raw + converted inputs,
    W_till_in/out, W_in/out, overburden_pressure, flux_u, flux_v)` — the full
    I-boundary contract trace. The dome cell shows `W_till` growing 0 → 0.0725 m
    (3 months of input absorbed by the till) and `W` draining 0.5 → 0.1488 m;
    the flux is radially outward (e.g. `(15,13)` at `y = -116 km` has
    `q = (0, -8.15e-4) m^2/s`).
  - `hydrology/flux_law.csv`: the pointwise flux law
    `K = k W^(alpha-1) (G^2 + eps^2)^((beta-2)/2)` (with the `eps = 1` for
    `beta < 2` regularization from `compute_conductivity`) and
    `|q| = K W G` over `(W, G = |grad R|)` grids — a pure `.esm` test target.
  - `hydrology/staggered.csv`: the internal substep subassembly
    (`Wstag, Kstag, Vstag, Qstag`, plus the accumulated `Qstag_average`) on the
    edge-centered grid, for later factoring of the substep update.
  - `hydrology/inputs.nc` has all inputs, state in/out, overburden pressure and
    the flux components.
  - Validation (see `runs/hydrology.json`): water-mass conservation closes to
    round-off both in the model's own accounting (0 kg residual) and in an
    independent state-based recomputation (~8e-15 relative); `W_till >= 0`;
    `overburden = rho_ice g H` exactly; the advective flux satisfies
    `q . grad R <= 0`; and the flux is centrally anti-symmetric to 1e-15
    (radially symmetric forcing). Note the dumped flux is the *advective* flux
    `V W_upwind`; the diffusive term `rho_w g K W grad W` is applied separately
    inside the `W` update and is not part of `flux()`.
- The test-F/V states mean every dumped quantity has an exact reference value
  (or a clean null state, as with age = 0); the discrete-vs-exact differences
  are the expected discretization error of the C++ model, which the `.esm`
  reimplementation must reproduce.

## Known C++ issue hit while dumping

`StressBalance::Inputs::dump()` (the natural input-dump hook) is broken in
v2.3.2: it calls `io::write_config()` on a `pism_config` variable that was
never defined, aborting with "NetCDF: Variable not found". Recorded in
`bugs.md` #4. The SIA driver writes `inputs.nc` directly instead.
