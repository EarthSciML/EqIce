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
| `bed/` | Bed deformation (PointwiseIsostasy, two update intervals, test-F load history) | `instrument/instrument_bed.cc` | `runs/bed.json` | `instrument/instrument_bed -Mx 31 -My 31` |
| `calving/` | Calving (Eigen / von Mises / Hayhurst / thickness / float-kill) + frontal melt (physics kernels + Constant) on the van der Veen CFBC shelf | `instrument/instrument_calving.cc` | `runs/calving.json` | `instrument/instrument_calving -Mx 65 -My 65` |
| `surface_ocean/` | Atmosphere + ocean + surface forcing boundaries (PIK atmosphere, PIK / Beckmann–Goosse ocean + Constant cross-check, PIK surface) | `instrument/instrument_surface_ocean.cc` | `runs/surface_ocean.json` | `instrument/instrument_surface_ocean -Mx 31 -My 31` |
| `geometry/` | Ice geometry update (GeometryEvolution `flow_step` + `source_term_step`, one dt step each on the test-F exact geometry) + pointwise `part_grid_threshold_thickness` | `instrument/instrument_geometry.cc` | `runs/geometry.json` | `instrument/instrument_geometry -Mx 31 -My 31 -dt_years 1` |

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
stage1/instrument/instrument_bed -Mx 31 -My 31 -dumps_dir stage1/dumps/bed
stage1/instrument/instrument_calving -Mx 65 -My 65 -dumps_dir stage1/dumps/calving
stage1/instrument/instrument_surface_ocean -Mx 31 -My 31 -dumps_dir stage1/dumps/surface_ocean
stage1/instrument/instrument_geometry -Mx 31 -My 31 -dt_years 1 -dumps_dir stage1/dumps/geometry
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
- **Bed deformation** (`bed/`): two full update intervals of
  `bed::PointwiseIsostasy` (`bed_deformation.model = "iso"`) on the test-F
  exact geometry, with a deterministic load history (uniform 300 m
  thickening over interval 1, constant in interval 2; 5 `update()` calls of
  dt = 20 "365day yr" per 100-yr365 interval). **Important:** `PointwiseIsostasy`
  is the *instantaneous* local isostasy law `topg_out = topg_last -
  f*(load - load_last)`, `f = rho_ice/rho_mantle`, with `load_last <- load`
  and **no relaxation time** — the exponential viscous half-space relaxation
  belongs to the Lingle-Clark model (`bed_deformation.model = "lc"`), which is
  a deferred subassembly (not instrumented). The `load` seen by `update_impl`
  is the *time-averaged* load over the update interval (the base class
  accumulates `load*dt` on every `update()` call and only triggers the model
  when `t_final` hits `m_t_last + update_interval`).
  - `bed/columns.csv`: per cell the full I-boundary trace — `(H0, H1,
    load_avg_1, load_avg_2, bed_in_1, bed_out_1, uplift_1, bed_in_2,
    bed_out_2, uplift_2, bed_analytic_1, bed_analytic_2)`. The dome cell
    (15,15) shows `bed_out_1 = -41.363636 m`, `bed_out_2 = -82.727273 m`
    (exactly `-f*150` and `-f*300`); ice-free far-field cells stay at 0.
  - `bed/pointwise.csv`: the two pointwise laws as pure `.esm` test targets —
    `compute_load(bed, H, sea_level)` (the ice-equivalent load, 0 where the
    ocean load exceeds the ice load, i.e. floating ice excluded) and the
    update law `topg_out = bed_in - f*(load - load_last)` with
    `f = rho_ice/rho_mantle`.
  - `bed/inputs.nc` has all inputs + state in/out for both intervals;
    `bed/load_history.csv` records the exact `(call, t, dt, dH)` load
    sequence.
  - Validation (see `runs/bed.json`): the pointwise contract holds to 0 m
    (machine precision), `uplift = (bed_out - bed_in)/dt_beddef` exactly
    (max err 1.65e-24 m/s), and the load accumulator reproduces the
    time-average of the call loads exactly (the load is `H0 + dH/2` over
    interval 1 and `H0 + dH` over interval 2).
- **Calving + frontal melt** (`calving/`): the calving and frontal-melt
  components run on the van der Veen CFBC shelf state (the same state as
  `ssa/`, but on a *square* grid `Mx = My = 65` because `EigenCalving`,
  `vonMisesCalving` and `HayhurstCalving` reject non-square grid cells in
  `init()`, and with 5 ocean columns beyond the calving front at `i = 60` so
  the `±2`-cell calving-rate sampling stays in the interior). The velocity the
  calving models act on is the actual SSA solution (0.062% average error vs the
  exact van der Veen solution).
  - `calving/eigen_calving.csv`: per cell `(H, mask, u, v, eigen1, eigen2,
    du_dx_analytic = C H^3, eigen1_avg, eigen2_avg, N_avg, rate)` for the
    exact-velocity run. `compute_2D_principal_strain_rates` gives `eigen1 =
    u_x` (matched to the analytic `C H^3` to 0.15% on interior cells) and
    `eigen2 = 0` *exactly* (1-D shelf), so the eigen-calving condition
    `eigen2 > 0 and eigen1 > 0` fails and the rate is 0 everywhere — a clean
    test of the compressive branch. Note the strain rates near the west
    Dirichlet boundary (i = 0, 1) are contaminated by the stale velocity ghost
    at `i = -1` (non-periodic boundary, as in the model).
  - `calving/eigen_calving_divergent.csv`: the same quantities for a
    manufactured divergent velocity `v = 0.05 u cos(2 pi y / L)` (periodic in
    y), exercising the nonzero branch: `rate = K * eigen1_avg * eigen2_avg`
    (with `K = calving.eigen_calving.K = 3e16 m s`, set from the config
    default 0.0) at 32 of the 65 front cells; max rate 297.8 m/yr. The dumped
    `eigen1_avg`/`eigen2_avg` reproduce the component's rate exactly.
  - `calving/vonmises_calving.csv`: per cell the eigen rates and the rate; at
    the 65 front cells also the reconstructed averaging `(|v|_avg, hardness_avg,
    e_s, sigma_tilde)` so the `.esm` target `rate = |v|_avg * sigma_tilde /
    sigma_max` with `sigma_tilde = sqrt(3) B e_s^(1/n)` is reproduced exactly
    (front cell `(60,32)`: `e_s = eigen1/sqrt(2) = 1.1785e-11 s^-1`,
    `sigma_tilde = 74891 Pa`, `rate = 71.13 m/yr`).
  - `calving/hayhurst_calving.csv`: per cell `(H, water_depth, omega, sigma_0,
    threshold, rate)` for the pointwise Hayhurst law (with the floating-shelf
    adjustment `omega = water_depth / (water_depth + (1 - rho_i/rho_w) H)`);
    rate range 1077–7906 m/yr over the 3900 icy cells. The component then
    (intentionally) propagates the mean of icy neighbors' rates to ice-free
    cells next to ice — verified to the last bit in `runs/calving.json`.
  - `calving/calving_at_thickness.csv`: `CalvingAtThickness::update(t=0,
    dt=1 yr)` with `calving.thickness_calving.threshold = 300 m` removes
    exactly the 65-cell front column (i = 59, `H = 188.77 m < 300`,
    floating, next to ice-free ocean): mask 3 → 4, thickness → 0.
  - `calving/float_kill.csv`: `FloatKill` with default config removes all 3900
    floating cells (mask 3 → 4); the analytic flotation thickness
    `h_f = (sea_level - bed) rho_w/rho_i = 1129.67 m` confirms every removed
    cell is floating.
  - `calving/frontal_melt_undercutting.csv` / `frontal_melt_ismip6.csv`: the
    two pointwise `FrontalMeltPhysics` kernels `(h, q_sg, TF) -> q_m` over
    deterministic grids (the undercutting kernel clamps negative/zero inputs to
    `q_m = 0`; the ISMIP6 kernel applies the raw power law). Pure `.esm` test
    targets with `q_m = (A h q_sg^alpha + B) TF^beta`.
  - `calving/frontal_melt_constant.csv`: the `Constant` frontal-melt model on
    the fully floating shelf. Default config (`include_floating_ice = no`)
    gives a zero melt rate (no grounded ice); with `include_floating_ice =
    yes` the melt rate is 1 m/day on icy cells and the front-cell retreat rate
    is `(rho_i/rho_w) * 1 = 0.885214 m/day` (`retreat_rate =
    H_submerged/H_threshold * melt`, floating front: `H_submerged = (rho_i/
    rho_w) H_threshold` with the part-grid threshold thickness).
  - `calving/inputs.nc` has the geometry, the SSA velocity, the enthalpy, and
    all three calving-rate outputs. `calving/thickness_threshold_input.nc` is
    the minimal NetCDF-3 input file the `CalvingAtThickness` constructor
    requires (it contains no `thickness_calving_threshold` variable, so the
    constant config threshold is used).
  - The front cell is 2 cells from the y-periodic edges; the `j ± 2` sampling
    there reads non-periodic ghost memory that happens to be ice-free, so the
    rates are identical to the interior front cells (verified: the dumps are
    byte-for-byte deterministic across reruns). This is a latent model edge
    case, not something the `.esm` reimplementation needs to reproduce.
- **Atmosphere / ocean / surface forcing** (`surface_ocean/`): the three
  forcing-boundary model families on one shared 31×31 grid — the test-F exact
  radial dome (extent 750 km) on a Gaussian ocean-trough bed
  `bed(r) = -1500 - 1500 exp(-((r-450km)/200km)^2)` m, with latitude
  `-70 + 10 y/Ly` deg north (a −80°..−60° latitudinal temperature gradient;
  69 grounded / 448 floating / 444 ice-free cells).
  - `atmosphere_temperature.csv`: the pointwise martin-2011 kernel
    `(latitude, usurf) -> T_ma = 303.15 + 0.68775·lat - 0.0075·usurf`
    (`T_ms = T_ma`, no seasonal cycle). The `air_temperature()` map matches
    this exactly (max err 0 K); `atmosphere_precipitation.csv` is the
    driver-prescribed field `P(r) = (200 + 100 cos(pi r/Lx))` kg m^-2 yr^-1
    (the real model reads it from a file), and
    `atmosphere_temperature_timeseries.csv` verifies the yearly-cycle time
    series is flat for martin (zero amplitude cosine).
  - `shelf_base_temperature.csv`: the pointwise pressure-melting-point kernel
    `T_shelf = T0 - beta_CC rho_i g H` (max err 0 K). **This is the linear
    pressure melting point and *decreases* with depth** — `T_shelf(2000) =
    271.7395 K` vs `T_shelf(0) = 273.15 K`; the naive "deep cavity → T_max"
    expectation is wrong for this model. The same value is used by the
    Constant ocean (cross-check, max err 0 K).
  - `beckmann_goosse_mass_flux.csv`: the pointwise kernel
    `(H, z_b, T_f, T_ocean, ocean_heat_flux) -> mass_flux = Q/L` with
    `Q = melt_factor rho_w c_p gamma_T (T_ocean - T_f)` and `z_b =
    -(rho_i/rho_w) H` (max err 1.7e-21 kg m^-2 s^-1, round-off). **Positive =
    melting** (the geometry update does `dH_BMB = -dt·mf/rho_i`); with the
    default `T_ocean = -1.7 °C` the flux is positive everywhere:
    `mf(0) = 1.229863e-06`, `mf(2000) = 9.501970e-06` kg m^-2 s^-1
    (= 0.3295 m/yr), range [1.23e-06, 1.36e-05]. The `ConstantPIK.cc` comment
    claiming the flux is "always negative" / "positive if freezing on" is a
    documentation bug (bugs.md #5).
  - `water_column_pressure.csv`: the pointwise kernel
    `P = 0.5 rho_w g h_w^2 / H` — the **depth-averaged** water column pressure
    at a margin (`h_w` = ocean depth below the ice base, 0 for `H = 0`),
    *not* `rho_w g depth` (max err 0 Pa; `P(2000 m, bed -3000) =
    7.902394e6 Pa`).
  - `surface_partition.csv`: `smb -> accumulation = max(smb,0)`,
    `melt = runoff = max(-smb,0)` (the `dummy_*` partition), so
    `mass_flux = accumulation - runoff` always. The driver wires the coupling
    that `surface::PIK` itself ignores (it holds SMB constant in real runs):
    SMB = atmosphere precipitation minus a Gaussian ablation blob centered on
    the warm-north floating shelf (peak 300 kg m^-2 yr^-1, sigma 250 km), so a
    32-cell band of the floating shelf has melt > 0. The surface temperature
    uses the *same* martin formula as the atmosphere, so the two maps are
    field-for-field identical (max err 0 K).
  - `columns.csv`: the per-cell trace of all maps plus the analytic references
    (`T_ma`, `T_shelf`, `mf`, `P`) the component outputs must match;
    `inputs.nc` has the geometry and all atmosphere/ocean/surface maps.
  - Validation (see `runs/surface_ocean.json`): air temperature, shelf base
    temperature, mass flux and water column pressure all match their analytic
    references to round-off; the Constant-ocean cross-check matches the PIK
    ocean on temperature and water column pressure and gives
    `mf = melt_rate·rho_i = 1.497038e-06` kg m^-2 s^-1 everywhere; the
    surface SMB partition identity and the surface==atmosphere temperature
    identity hold to 0; all 32 melting cells are on the floating shelf.
- **Geometry / mass continuity** (`geometry/`): the two integrated steps of
  `GeometryEvolution` (`flow_step` — advective + diffusive mass transport —
  and `source_term_step` — surface + basal mass balance), each over `dt = 1 yr`
  from the SAME test-F exact geometry (flat bed, `H` = exact FG, sea level 0:
  517 grounded + 444 ice-free bedrock cells), plus the pointwise
  `part_grid_threshold_thickness`.
  - The forcing is exactly specified: the advective velocity is the exact
    test-F **surface** velocity `U_surf(r)·(x/r, y/r)` (a cold no-slide SIA
    solution, so these speeds are small, max ~2.6 m/yr — the advective part is
    a small correction); the diffusive flux is a prescribed Gaussian radial
    field `q = q0·(r/Rq)·exp(-(r/Rq)^2)·e_r` with `q0 = 2e-2 m^2/s`,
    `Rq = 300 km` on the staggered edge midpoints (a documented, calibrated
    surrogate for the real test-F SIA flux); the SMB is the exact test-F
    `M(r)·rho_ice`; the basal melt rate is `0.05 m/yr·exp(-(r/200 km)^2)`; the
    thickness BC mask is 0 everywhere.
  - `geometry/columns.csv`: per cell the full I-boundary trace for both steps
    — `(H, surface, bed, mask)`, the forcing `(u, v, qx, qy, smb_rate,
    basal_melt_rate)`, and the geometry out + thickness change after each step.
  - `geometry/flow_inputs.nc` / `geometry/source_inputs.nc`: the complete
    in/out NetCDF traces (geometry in, forcing, the model's internal fields —
    `flux_staggered`, `flux_divergence`, `thickness_change`,
    `conservation_error`, `effective_SMB`, `effective_BMB` — and the geometry
    out with `ensure_consistency` applied, named `thk_flow_out`/`usurf_flow_out`/
    `mask_flow_out` and `thk_source_out`/`usurf_source_out`/`mask_source_out`).
  - `geometry/pointwise.csv`: 560 samples of
    `part_grid_threshold_thickness(cell_type, thickness, surface, bed)` over
    7 mask patterns × `(H, h, bed)` grids — a pure `.esm` test target. The
    function averages `H` and `h` over the icy N/S/E/W neighbors and returns
    `max(min(h_avg - bed, H_avg), 0)` (0 if no icy neighbors); all 560 rows
    match this closed form to < 1e-14.
  - Validation (see `runs/geometry.json`): the flow_step domain-sum mass
    balance closes to round-off (`sum(H_change) = 6.59e-12 m`; the flux
    limiters zero the fluxes on the ice-free outer ring, so `sum(div Q)`
    telescopes to 0); the contract residual
    `sum(dH + dHref - (-dt divQ) - conservation_error) = 6.6e-12 m`; the dome
    cell `(15,15)` thins by `dH = -4.25123 m = -(4.1684 diffusive + 0.0829
    advective)`; the flow advances the ice front one cell (76 cells go mask
    0 → 2 with `H_out` in [0.033, 0.093] m). The source step sums
    `dH_SMB = 0.0233458 m` and `dH_BMB = -1.86362 m`, each exactly equal to
    `dt·sum(smb/rho)` and `-dt·sum(bmr)` (the only non-negativity truncation
    is 1.7e-6 m at the 444 ice-free cells); the dome gains
    `dH = 0.0353269 m = dt·M(0) - dt·0.05 m/yr` exactly.

- The test-F/V states mean every dumped quantity has an exact reference value
  (or a clean null state, as with age = 0); the discrete-vs-exact differences
  are the expected discretization error of the C++ model, which the `.esm`
  reimplementation must reproduce.

## Known C++ issue hit while dumping

`StressBalance::Inputs::dump()` (the natural input-dump hook) is broken in
v2.3.2: it calls `io::write_config()` on a `pism_config` variable that was
never defined, aborting with "NetCDF: Variable not found". Recorded in
`bugs.md` #4. The SIA driver writes `inputs.nc` directly instead.
