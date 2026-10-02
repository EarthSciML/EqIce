# PROGRESS — eqice implementation

Tracking progress on PLAN.md (Stage 1 — instrument the C++ codes; Stage 2 —
stubs → physics → EarthSciModels PRs). Last updated: 2026-10-01.

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

## Milestone 3 (continued) — bed deformation boundary instrumented ✅

- **Bed deformation traces** ✅ — new `instrument_bed.cc`: two full update
  intervals of `bed::PointwiseIsostasy` (PISM config `bed_deformation.model =
  "iso"`) on the test-F exact geometry (flat bed, Mx=My=31). Deterministic
  load history: interval 1 (t = 0..100 "365day yr", 5 calls of dt = 20 yr365)
  grows the ice load uniformly from H0(r) to H0(r)+300 m in 5 equal steps, so
  the time-averaged load is exactly H0 + 150 m; interval 2 (t = 100..200
  yr365) holds the load constant at H0 + 300 m.
  - **Key finding:** `PointwiseIsostasy` is the *instantaneous* local isostasy
    law `topg_out = topg_last - f*(load - load_last)`,
    `f = rho_ice/rho_mantle = 910/3300 = 0.275758`, with `load_last <- load`
    and **no relaxation time** — a load change is fully compensated in the
    single update that follows it. The exponential viscous half-space
    relaxation belongs to the Lingle-Clark model (`bed_deformation.model =
    "lc"`), which is **not** instrumented (deferred subassembly; the
    boundaries.md §8 catalog and the dumps README note this).
  - The load fed to `update_impl` is the *time-averaged* load over the update
    interval: `BedDef::update()` accumulates `load*dt` on every call and only
    triggers the derived model when `t_final` hits `m_t_last +
    update_interval` (within `time_stepping.resolution`); it then uses
    `load = accumulator/dt_beddef`, resets the accumulator, and sets
    `uplift = (topg - topg_last)/dt_beddef`, `topg_last = topg`.
  - Validation: `bed_out_1 = -f*150 = -41.3636 m` and
    `bed_out_2 = -f*300 = -82.7273 m` exactly (max |bed - analytic| = 0 m over
    the 517 ice cells); 444 ice-free far-field cells stay at 0; the uplift
    identity `uplift = (bed_out-bed_in)/dt_beddef` holds to 1.65e-24 m/s;
    `uplift = -1.3116e-8 m/s = -0.4136 m/365day-yr` per interval; the load
    accumulator reproduces the time average of the call loads exactly
    (H0 + 150 m / H0 + 300 m).
  - Dumps in `stage1/dumps/bed/`: `meta.txt`, `parameters.csv`, `inputs.nc`
    (loads, bed in/out, uplift, analytic references for both intervals),
    `columns.csv` (per-cell I-boundary trace), `pointwise.csv` (the
    `compute_load` and update-law pointwise functions — analytic .esm test
    targets), `load_history.csv` (the exact (call, t, dt, dH) sequence).
- **Run config** ✅ — `stage1/runs/bed.json` (grid, dt, load history, checks).
- The `Null` (`bed_deformation.model = "none"`) cross-check is trivially
  verified from source (empty `update_impl`: bed unchanged, uplift 0), not
  dumped. The `Given` and `LingleClark`/`LingleClarkSerial` models are not
  instrumented (deferred; noted in `stage1/boundaries.md` §8).

## Milestone 3 (continued) — calving / front retreat + frontal melt instrumented ✅

- **Calving traces** ✅ — new `instrument_calving.cc` on the van der Veen CFBC
  shelf state (same state as `instrument_ssa.cc`, but on a *square* grid
  Mx=My=65: the rate-based calving components reject non-square cells in
  `init()`, so the SSA run's 61x3 grid cannot be used — see
  `stage1/runs/calving.json`; 5 ocean columns beyond the front keep the ±2-cell
  calving sampling in the interior). The velocity the calving models act on is
  the actual SSA solution (0.062% avg error vs exact).
  - **EigenCalving** (D): dumps the strain-rate invariants `eigen1`/`eigen2`
    from `compute_2D_principal_strain_rates` plus the rate. With the exact
    shelf velocity `eigen2 = 0` *exactly* (1-D shelf), so the rate is 0
    everywhere (compressive branch); `eigen1` matches the analytic
    `du/dx = C H^3` to 0.15% on interior cells (the i=0,1 boundary strain
    rates are contaminated by the stale velocity ghost at i=-1, as in the
    model). A manufactured divergent velocity `v = 0.05 u cos(2 pi y / L)`
    exercises the nonzero branch `rate = K * eigen1_avg * eigen2_avg` (K =
    3e16 m s, set from the config default 0.0): max rate 297.8 m/yr at 32 of
    65 front cells; the reconstructed formula matches the component to 0.
  - **vonMisesCalving** (D, clean pointwise rate — not deferred): it uses the
    2-D strain-rate invariants + the Glen flow-law averaged hardness
    (`sigma_tilde = sqrt(3) B e_s^(1/n)`), no iterative 3-D stress solve.
    Front-cell `(60,32)`: `e_s = 1.1785e-11 s^-1`, `sigma_tilde = 74891 Pa`,
    `rate = 71.13 m/yr`; reconstruction matches the component to 0.
  - **HayhurstCalving** (D, clean pointwise): rate range 1077–7906 m/yr over
    the 3900 icy cells; the floating-shelf omega adjustment and the
    `sigma_0 -> max(sigma_0, threshold)` clamp verified. The component's
    intended propagation of the mean icy-neighbor rate to ice-free cells next
    to ice verified to 0.
  - **CalvingAtThickness** (I, in-place): with threshold 300 m removes exactly
    the 65-cell front column (H = 188.8 m, floating, next to ice-free ocean);
    mask 3 -> 4, thickness -> 0 dumped per cell.
  - **FloatKill** (I, in-place): removes all 3900 floating cells (default
    config); flotation thickness h_f = 1129.67 m confirms the removal set.
- **Frontal melt traces** ✅ — same driver:
  - `FrontalMeltPhysics::frontal_melt_from_undercutting` / `_from_ismip6`
    pointwise kernels `(h, q_sg, TF) -> q_m` over deterministic grids (pure
    `.esm` targets; `q_m = (A h q_sg^alpha + B) TF^beta`); undercutting clamps
    negative/zero inputs to 0, ISMIP6 does not.
  - `frontalmelt::Constant` on the fully floating shelf: default config
    (`include_floating_ice = no`) -> 0 melt (no grounded ice); with
    `include_floating_ice = yes` -> 1 m/day on ice and a front-cell retreat
    rate `(rho_i/rho_w) * 1 = 0.885214 m/day`.
- **Dumps** ✅ — `stage1/dumps/calving/` (13 files, ~4.8 MB): `meta.txt`,
  `parameters.csv`, `inputs.nc`, `eigen_calving.csv`, `eigen_calving_divergent.csv`,
  `vonmises_calving.csv`, `hayhurst_calving.csv`, `calving_at_thickness.csv`,
  `float_kill.csv`, `frontal_melt_undercutting.csv`, `frontal_melt_ismip6.csv`,
  `frontal_melt_constant.csv`, plus `thickness_threshold_input.nc` (the minimal
  input file `CalvingAtThickness`'s constructor requires; no threshold variable,
  so the constant config threshold is used).
- **Run config** ✅ — `stage1/runs/calving.json` (grid, config overrides,
  dump files, check values).
- Deterministic: dumps byte-for-byte identical across reruns.
- Deferred (documented in boundaries.md §9/§6.4 and the dumps README):
  `FrontRetreat::update_geometry` (retreat-rate application over dt, **I**),
  `PrescribedRetreat` (**I**, ISMIP6 parameterized retreat mask), the
  `Given`/`DischargeGiven`/`DischargeRouting` frontal-melt models, and the
  `calving.rate_scaling` modifier.
- **Grid note (important for Stage 2):** the flow-line 61x3 grid of the SSA
  dump cannot be used for the calving components — `EigenCalving`,
  `vonMisesCalving` and `HayhurstCalving` all throw on non-square cells
  (`|dx-dy|/min(dx,dy) > 1e-2`). The calving dumps use a square 65x65 grid with
  the same uniform-in-y shelf state.

## Milestone 3 (continued) — atmosphere / ocean / surface forcing boundaries instrumented ✅

- **Surface/ocean forcing traces** ✅ — new `instrument_surface_ocean.cc`: the
  three forcing-boundary model families (PIK atmosphere, PIK / Beckmann–Goosse
  ocean + Constant cross-check, PIK surface) on one shared 31×31 grid. The
  geometry is the test-F exact radial dome (extent 750 km) on a Gaussian ocean
  trough/cavity bed `bed(r) = -1500 - 1500 exp(-((r-450km)/200km)^2)` m
  (r = 450 km → −3000 m trough under the mid-shelf) with sea level 0 and a
  latitude gradient `-70 + 10 y/Ly` deg (spanning −80..−60, so the south edge
  is cold and the north edge warm): 69 grounded / 448 floating / 444
  ice-free-ocean cells. Models constructed directly (not via the factory); the
  PIK atmosphere and surface are subclassed in-driver to prescribe the
  precipitation / SMB fields instead of reading files.
  - **Atmosphere** (atmosphere::PIK, default `martin` parameterization):
    `air_temperature() = T_ma = 303.15 + 0.68775·lat - 0.0075·usurf`
    (note the source's `-0.68775*lat*(-1.0)` = `+0.68775·lat`), verified to the
    analytic formula with max err 0 K over the map (`T_ma(-80,0) = 248.13 K`,
    `T_ma(-60,0) = 261.885 K`; dome cell (15,15): lat −70, usurf 1490.51 m →
    243.829 K). `T_ms = T_ma` for martin, so the yearly-cycle time series is a
    flat zero-amplitude cosine (verified through the component's own
    `temp_time_series`). Precipitation is driver-prescribed
    `P(r) = (200 + 100 cos(pi r/Lx)) kg m^-2 yr^-1` (the real model reads it
    from a file; `atmosphere.pik.file` left empty).
  - **PIK ocean** (ocean::PIK, Beckmann–Goosse-style): `shelf_base_temperature`
    is the **linear** pressure melting point `T0 - beta_CC rho_i g H` — it
    *decreases* with depth (`T_shelf(0) = 273.15 K`, `T_shelf(2000) =
    271.7395 K`), not the naive "deep cavity → T_max" (task guess corrected;
    max err 0 K). `shelf_base_mass_flux` is `Q/L` with
    `Q = melt_factor rho_w c_p gamma_T (T_ocean - T_f)`, `T_f = 273.15 + 0.0939
    - 0.057 S + 7.64e-4 z_b`, `z_b = -(rho_i/rho_w) H`; **positive = melting**
    (the geometry update path is `combine_basal_melt_rate` → `dH_BMB =
    -dt·mf/rho_i`). With the default `T_ocean = -1.7 °C` the flux is positive
    everywhere: `mf(0) = 1.229863e-06`, `mf(2000) = 9.501970e-06` kg m^-2 s^-1
    (= 0.3295 m/yr), range [1.23e-06, 1.36e-05]; max err vs analytic
    1.7e-21 kg m^-2 s^-1 (round-off). `average_water_column_pressure` is the
    **depth-averaged** form `0.5 rho_w g h_w^2/H` (NOT `rho_w g depth`), max err
    0 Pa (`P(2000 m ice, bed −3000) = 7.902394e6 Pa`; range [0, 9.979e6] Pa).
  - **Constant ocean cross-check** (ocean::Constant): `shelf_base_temperature`
    identical to the PIK pressure melting point (max err 0 K);
    `shelf_base_mass_flux = melt_rate·rho_i = 1.497038e-06` kg m^-2 s^-1
    everywhere (the default `0.051914 m/yr`); water column pressure identical
    to PIK (max err 0 Pa).
  - **Surface** (surface::PIK, wired by the driver): the SMB is the atmosphere
    precipitation minus a deterministic Gaussian ablation blob centered on the
    warm-north floating shelf (x=0, y=+450 km, peak 300 kg m^-2 yr^-1, sigma
    250 km), so a 32-cell band of the north shelf has melt > 0 while the rest
    accumulates. The `dummy_*` partition is `accumulation = max(smb,0)`,
    `melt = runoff = max(-smb,0)`, so `mass_flux = accumulation - runoff`
    exactly (max err 0). The surface temperature uses the same martin formula
    as the atmosphere, so `surface temperature == atmosphere air_temperature`
    field-for-field (max err 0 K).
- **Dumps** ✅ — `stage1/dumps/surface_ocean/` (11 files): `meta.txt`,
  `parameters.csv`, `inputs.nc` (geometry + all atmosphere/ocean/surface maps),
  `columns.csv` (per-cell trace with analytic references),
  `atmosphere_temperature.csv`, `atmosphere_precipitation.csv`,
  `atmosphere_temperature_timeseries.csv`, `shelf_base_temperature.csv`,
  `beckmann_goosse_mass_flux.csv`, `water_column_pressure.csv`,
  `surface_partition.csv`.
- **Run config** ✅ — `stage1/runs/surface_ocean.json` (grid, geometry, model
  construction, config keys that matter, check values). Key config:
  `atmosphere.pik.parameterization = "martin"` (default),
  `atmosphere.pik.file = ""` (driver prescribes precipitation),
  `ocean.pik_melt_factor = 5e-3`, `ocean.constant.melt_rate = 0.0519 m/yr`
  (default); the ocean-side hard-coded constants are in `ConstantPIK.cc`:
  `c_p_ocean = 3974`, `gamma_T = 1e-4`, `salinity = 35`, `T_ocean = -1.7 °C`.
- **New C++ documentation bug found** — the comment in `ConstantPIK.cc` claims
  the shelf base mass flux is "always negative" and "positive if ice is
  freezing on"; the actual sign convention is **positive = melting** (bugs.md
  #5).
- Deferred (documented in boundaries.md §10): `atmosphere::Given`/
  `atmosphere::FifoSIA`/`atmosphere::Forcing`, `ocean::Given`, the
  `surface::Given` and `surface::Delta_T` SMB models, and the coupler-level
  `SurfaceForcing` time interpolation.
- Deterministic: dumps byte-for-byte identical across reruns.

## Milestone 3 (continued) — geometry (mass continuity) boundary instrumented ✅

- **Geometry traces** ✅ — new `instrument_geometry.cc`: the two integrated
  steps of `GeometryEvolution` on the test-F exact geometry (flat bed,
  H = exact FG, sea level 0, Mx=My=31; 517 grounded + 444 ice-free bedrock
  cells), each over `dt = 1 yr` from the same input state, plus a pointwise
  sample of `part_grid_threshold_thickness`.
  - **flow_step** (I-boundary trace): advective velocity = the exact test-F
    **surface** velocity `U_surf(r)·(x/r, y/r)` (analytic from `exactFG`;
    test F is a cold no-slide SIA solution, so these speeds are small — max
    2.58 m/yr — and the advective part is a small correction); diffusive flux
    = a prescribed Gaussian radial field `q = q0·(r/Rq)·exp(-(r/Rq)^2)·e_r`
    (`q0 = 2e-2 m^2/s`, `Rq = 300 km`) on the staggered edge midpoints — a
    documented, exactly-specified surrogate for the real test-F SIA flux
    (max |q_SIA| ~ 2e-2 m^2/s, so the magnitude is calibrated), with the
    continuum divergence `div q = (2 q0/Rq) exp(-s)(1-s)` as an analytic
    target. `thickness_bc_mask = 0` everywhere.
    - Validation: the domain-sum mass balance closes to round-off —
      `sum(H_change) = 6.59e-12 m`, `sum(conservation_error) = 0 m`, and
      `sum(-dt·div Q) = -2.24e-14 m` (= −7.5e-5 m³): because the outermost
      ring of the domain is ice-free, the flux limiters zero every flux on
      the outer boundary, so `sum(div Q)` telescopes to 0. The contract
      residual `sum(dH + dHref − (−dt divQ) − conservation_error) = 6.6e-12 m`.
    - The dome cell (15,15) thins by `dH_flow = -4.25123 m` =
      `−(4.1684 diffusive + 0.0829 advective)` m (reconstructed from the
      dumped edge fluxes exactly; the continuum diffusive-only estimate is
      `-dt·2q0/Rq = -4.20759 m`).
    - The flow advances the ice front one grid cell: **76 cells change mask
      0 (ice-free bedrock) → 2 (grounded ice)** in the margin ring at
      r ∈ [754.8, 789.8] km, each gaining 0.033–0.093 m of ice over the year
      (the upwind advective + diffusive inflow crosses the 0.01 m ice-free
      threshold). No cell thins below zero (conservation_error = 0).
  - **source_term_step** (I-boundary trace): smb_rate = the exact test-F
    surface mass balance `M(r)·rho_ice` (analytic); basal_melt_rate =
    `0.05 m/yr·exp(-(r/200 km)^2)`; same `thickness_bc_mask = 0`.
    - `sum(dH_SMB) = 0.0233458 m` and `sum(dH_BMB) = -1.86362 m`, each
      **exactly** equal to `dt·sum(smb/rho_ice)` and `-dt·sum(bmr)`
      (0 SMB-clamped cells); the only non-negativity truncation is 1.7e-6 m at
      the 444 ice-free cells (BMB-clamped: nothing to melt there).
    - Dome cell: `dH_SMB = 0.0853269 m = dt·M(0)` exactly (M(0) = 0.0853 m/yr,
      the exact test-F SMB near the dome), `dH_BMB = -0.05 m = -dt·0.05 m/yr`
      exactly; `dH_source = 0.0353269 m`. No mask changes.
  - **part_grid_threshold_thickness** (D pointwise): 560 samples over
    7 mask patterns (all grounded / all floating / grounded+floating+ocean /
    3 icy+1 ocean / 1 icy+3 ice-free / all ice-free ocean / all ice-free land)
    × (H ∈ {0.1, 100, 1000, 3000} m, h ∈ {−500, 0, 500, 2000, 4000} m,
    bed ∈ {−1000, 0, 500, 2000} m). **All 560 rows match the closed form
    `max(min(h_avg − bed, H_avg), 0)` (0 if no icy neighbors) to < 1e-14** —
    a pure analytic .esm test target.
  - The part-grid scheme itself is inert in the steps (no ice-free-ocean cells
    in the flat-bed test-F state, so `geometry.part_grid.enabled = no` as
    default); the pointwise function is sampled separately, and the flow/source
    steps use the flat-bed fully-grounded state as the trace geometry.
- **Dumps** ✅ — `stage1/dumps/geometry/` (6 files, ~1.0 MB): `meta.txt`,
  `parameters.csv`, `flow_inputs.nc`, `source_inputs.nc` (complete in/out
  traces with distinct `*_flow_out`/`*_source_out` variable names),
  `columns.csv` (per-cell trace for both steps), `pointwise.csv`.
- **Run config** ✅ — `stage1/runs/geometry.json` (grid, dt, forcing choice +
  formulas, check values).
- Deterministic: dumps byte-for-byte identical across reruns.
- This completes the last Stage-1 boundary in PLAN.md §10 milestone 3.

## Stage 1 — Remaining instrumentation

- [x] Step 1: identify subassembly boundaries → `stage1/boundaries.md`
- [x] Step 2: instrument boundaries to dump inputs/outputs — **rheology + SIA +
      SSA + energy + age + basal strength + hydrology (Routing) + bed
      (PointwiseIsostasy) + calving/frontal melt + surface/ocean (PIK
      atmosphere, PIK/Beckmann–Goosse ocean + Constant cross-check, PIK
      surface) + geometry (GeometryEvolution flow_step/source_term_step +
      part_grid_threshold_thickness) done**; remaining (deferred):
      Distributed/SteadyState hydrology, LingleClark/Given bed models
- [x] Step 3 (partial): dump subassembly I/O — SSA `nuH`/`taud` (FD
      subassemblies), SIA `delta`/`D`/`q` chains, the
      `MohrCoulombPointwise`/basal-resistance-law pointwise functions, the
      Routing staggered `Wstag/Kstag/Vstag/Qstag` substep fields, the bed
      `compute_load`/update-law pointwise functions, the
      atmosphere/ocean/surface pointwise kernels
      (`(lat,usurf)->T_ma`, `H->T_shelf`, `(H,z_b,T_f,T_ocean)->mf`,
      `(H,bed)->P_awcp`, `smb->(accumulation,melt,runoff)`), and the geometry
      subassemblies (`flux_staggered`, `flux_divergence`, `thickness_change`,
      `conservation_error`, `effective_SMB`/`effective_BMB`, and the
      `part_grid_threshold_thickness` pointwise function) are dumped; deeper
      subassemblies (e.g. the assembled KSP matrix) are deferred
- [x] Step 4: record run configurations → `stage1/runs/README.md` +
      `stage1/runs/{rheology,sia,ssa,energy_age,basal_strength,hydrology,bed,calving,surface_ocean,geometry}.json`;
      dumps stored in `stage1/dumps/`
- [ ] Step 5: extract select dumps into numeric test tuples (feeds Stage 2)

## Milestone 4 — Stage 2: ice-rheology component authored and passing ✅

The first Stage-2 deliverable: the PISM ice-rheology boundary as a
hand-authored, compositional `.esm` component built from the Stage-1 rheology
dumps (`stage1/dumps/rheology/`).

- **Deliverable** ✅ — `stage2/ice_rheology.esm` (10 models: EnthalpyConverter,
  EffectiveViscosity, SecondInvariant, and the 7 FlowLawFactory flow laws —
  isothermal Glen, Paterson–Budd, Arrhenius cold/warm, GPBLD, Hooke,
  Goldsby–Kohlstedt) plus `stage2/ice_rheology_templates.esm` (the shared
  expression-template library, imported by reference per model).
- **Physics** ✅ — matches the instrumented C++ v2.3.2 (`stage1/dumps/rheology/`):
  - EnthalpyConverter: `E -> (T, T_pa, omega, is_temperate, E_cts, E_l, T_m)`
    with exact config constants (`T_melting = 273.15`, `L = 3.335e5`,
    `c_p_ice`, `beta_CC_grad = 7.9e-4 K/m`, `Tm_0 = 273.15`).
  - Flow laws: `softness = A` (Arrhenius `A(T_pa)` for pb/arr/arrwarm/gk;
    Hooke `A(T)`; GPBLD temperate branch `A(T_melting)·(1 + 181.25·min(ω, 0.01))`),
    `hardness = B = A^(-1/n)` (gk `B = A_pb(T_pa)^(-1/n)`), and the strain rate
    `flow = A·sigma^(n-1)` (gk via the harmonic combination of dislocation,
    diffusional, and grain-boundary-sliding creep; grain-size sweep tested).
  - EffectiveViscosity: `nu = (0.5 A)^(-1/n) eps_2D^((1-n)/(2n))` + `dnu`;
    SecondInvariant: `eps^2 = 0.5 eps_ij eps_ij` with `w_z = -(u_x + v_y)`.
- **Composition** ✅ — every shared calculation (L(T), Arrhenius A, Hooke A,
  hardness, flow, effective viscosity, second invariant, …) is an expression
  template in `ice_rheology_templates.esm`, imported by reference via
  `expression_template_imports`; match-based op-call rules keep the equations
  readable as math (e.g. `hardness(T, p)` with the model's own unknowns). No
  scripts generated equations; only the numeric test tuples were extracted
  from the C++ dumps.
- **Tests** ✅ — 62 inline test groups, 313 assertions, passing under both the
  native and interpreter compilers (`./esm test stage2/ice_rheology.esm` and
  `--compiler interpreter`): 8 EnthalpyConverter tuples (all 7 outputs), 12
  effective-viscosity (6 configs × `eps = 0` / Schoof `eps_Schoof`), 9
  second-invariant, 4–5 per flow law (softness/hardness/flow/T/T_pa/ω/
  is_temperate), 7 Goldsby–Kohlstedt (incl. grain-size sweep). Relative
  tolerances 1e-12 (1e-10 Hooke, 1e-9 GK).
- **Verification** ✅ — `./esm validate` clean on both files; `./esm info`
  lists the 10 models; `./esm units --check`: 43 consistent / 0 mismatched /
  20 not checked (the not-checked entries are the intentional non-literal `^`
  exponents with parameter `n` — no unit mismatches).
- **Deferred** — `AveragedHardness` (the flow-law factory's depth-averaged
  hardness): the Stage-1 dumps are not self-contained (they depend on the full
  column solution), so it waits for the SIA/SSA components' column traces.
- **C++ bugs** — none new in the rheology boundary. Two by-design quirks
  (not bugs): the arr cold flow uses unadjusted `T` while softness uses
  `T_pa`; gk `softness()` throws in C++ (dump NaN), so the gk model omits it
  and exposes `eps_disl/eps_diff/eps_gbs/eps_basal` instead. `bugs.md`
  unchanged.

## Stage 2 — Stubs → physics → EarthSciModels PRs (in progress)

Per PLAN.md §5: hand-author stub `.esm` files with `tests` blocks from
Stage-1 tuples, fill in physics, review + merge into EarthSciModels.
Rheology (milestone 4 above) is done; the remaining boundaries follow the
PLAN.md §5 sequencing (SIA/SSA → energy → basal/hydrology → bed →
surface/ocean/calving → subassemblies).

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
- Known issue: `ConstantPIK.cc` shelf-base mass flux sign comment wrong
  (bugs.md #5; the flux is positive = melting)

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
- `stage1/dumps/bed/` — 6 files, ~0.8 MB (regenerate:
  `instrument/instrument_bed -Mx 31 -My 31 -dumps_dir
  stage1/dumps/bed`)
- `stage1/dumps/calving/` — 13 files, ~4.8 MB (regenerate:
  `instrument/instrument_calving -Mx 65 -My 65 -dumps_dir
  stage1/dumps/calving`)
- `stage1/dumps/surface_ocean/` — 11 files, ~0.8 MB (regenerate:
  `instrument/instrument_surface_ocean -Mx 31 -My 31 -dumps_dir
  stage1/dumps/surface_ocean`)
- `stage1/dumps/geometry/` — 6 files, ~1.0 MB (regenerate:
  `instrument/instrument_geometry -Mx 31 -My 31 -dt_years 1 -dumps_dir
  stage1/dumps/geometry`)
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
