# PROGRESS — eqice implementation

Tracking progress on PLAN.md (Stage 1 — instrument the C++ codes; Stage 2 —
stubs → physics → EarthSciModels PRs). Last updated: 2026-10-02.

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

## Milestone 4 (continued) — Stage 2: SIA stress-balance component authored and passing ✅

The second Stage-2 deliverable: the PISM shallow-ice-approximation (SIA)
stress-balance boundary as a hand-authored, compositional `.esm` component
built from the Stage-1 SIA dumps (`stage1/dumps/sia/`).

- **Deliverable** ✅ — `stage2/sia.esm` (4 models: SIADiffusivityColumn,
  SIADiffusiveFlux, SIA3DHorizontalVelocity, SIAVerticalVelocity) plus
  `stage2/sia_templates.esm` (SIA expression-template library: `sia_pressure`,
  `sia_effective_stress`, `sia_flow_arr`, `sia_delta`, `trapezoid_step`),
  imported by reference per model together with
  `./ice_rheology_templates.esm`.
- **Physics** ✅ — matches the instrumented C++ v2.3.2 (`stage1/dumps/sia/`):
  - vertical chain per level: `pressure = rho_ice*g*depth`,
    `stress = alpha*p`, `flow = A_cold*exp(-Q_cold/(R*T))*stress^(n-1)` with
    `T = enthalpy_to_temperature(E,p)` (the standard converter), and
    `delta = e_factor*theta*2*p*flow` — verified bit-exact against the dumped
    `delta` column (the authoritative one, see bugs.md #6);
  - `I(z) = ∫_0^z delta ds` (cumulative trapezoid, prefix-reduction faq) and
    `D = ∫_0^H depth*delta dz` (trapezoid with the partial top interval)
    — verified to 1e-12 (dump bit-exact);
  - diffusive flux `q = -D grad h` (o=0 → q_x, o=1 → q_y; the instrument's
    diffusive_flux.csv pairs the two staggered components on one row, so the
    `.esm` tests assert the per-point vector `-D*h_x, -D*h_y`);
  - 3D horizontal velocity `U = sliding_velocity - 0.25 Σ_{e,w,n,s} I·∇h`
    over the four staggered neighbors — reproduced bit-exact (max err 0.0);
  - vertical velocity `w3` from incompressibility
    `w[k] = w[k-1] - 0.5*dz*(div[k]+div[k-1])`, `w[0]=0` — reproduced to
    4e-25 relative.
- **Composition** ✅ — every shared calculation is an expression template in
  `sia_templates.esm` / `ice_rheology_templates.esm`, imported by reference;
  `sia_flow_arr` composes rheology templates (`flow_power_law`,
  `arrhenius_softness`, `enthalpy_to_temperature`) inside its body. The
  vertical integrals use the shared `trapezoid_step` template and faq
  prefix-reduction / sum-with-filter forms. `theta` (bed smoothness) is an
  input parameter; the BedSmoother is deferred to stage 3. No scripts
  generated equations — only numeric arrays/tuples were extracted from the
  dumps into the hand-authored `tests` blocks.
- **Tests** ✅ — 48 assertions, passing under both the native and interpreter
  compilers (`./esm test stage2/sia.esm` and `--compiler interpreter`):
  3 diffusivity-column tests ((20,15,0), dome (15,15,0), margin (10,3,0))
  covering D, delta, I, stress, pressure, T; 4 diffusive-flux tests; 1
  3D-horizontal-velocity test at cell (20,16) (u3/v3 at 5 levels); 1
  vertical-velocity test at (20,15) (w3 at 6 levels). Relative tolerances
  1e-12.
- **Verification** ✅ — `./esm validate` clean on both files; `./esm info`
  lists the 4 models; `./esm units --check`: 8 consistent / 0 mismatched /
  5 not checked (non-literal `^` exponent with parameter `n`, same
  informational pattern as ice_rheology.esm); ice_rheology.esm regression
  still 313/313.
- **C++ bugs** — new entry in `bugs.md` #6: the instrumented delta.csv `flow`
  column is computed with a **Cold** EnthalpyConverter (`T_melting` ~1e6 K),
  inconsistent with the authoritative `delta`/`D` columns at temperate
  columns (e.g. (15,15,0,k=0): dumped flow 6.950e-16 vs true 6.069e-16
  Pa⁻¹ s⁻¹). The `.esm` reproduces the authoritative chain with the standard
  converter.
- **CLI learning** — matched expression templates auto-expand inside faq
  bodies (proved via the exp6/7/8 and exp43/44 construction chain); a
  hand-authored column model with `enthalpy_to_temperature` in a faq body
  fails to compile unless the template library is actually imported via
  `expression_template_imports` (the missing-import failure mode:
  `unlowered_operator: rewrite-target operator 'enthalpy_to_temperature'
  reached evaluation`).

## Milestone 4 (continued) — Stage 2: SSA stress-balance component authored and passing ✅

The third Stage-2 deliverable: the PISM shallow-shelf-approximation (SSA)
stress-balance boundary as a hand-authored, compositional `.esm` component
built from the Stage-1 SSA dumps (`stage1/dumps/ssa/`), flow-line form of
verification test V (van der Veen / CFBC floating shelf).

- **Deliverable** ✅ — `stage2/ssa.esm` (model SSAFlowline: the 1-D membrane
  balance on the 60-cell test-V shelf) plus `stage2/ssa_templates.esm` (the
  SSAFD expression-template library: `ssa_face_stiffness`,
  `ssa_surface_gradient_x`, and the membrane-divergence rule
  `ssa_flowline_membrane_divergence`), importing `effective_viscosity` from
  `./ice_rheology_templates.esm` by reference. The templates declare the
  flow-line geometry contract themselves (metaparameter N=60, index sets x
  and face) rather than importing the cartesian_uniform_1d grid library, so
  the file is self-contained for the eventual EarthSciDiscretizations PR.
- **PDE as written** ✅ — the membrane momentum balance
  div(4·ν·H·grad u) = ρ·g·H·h_x with PISM's SSAFD discretization:
  `fluxdiv = D(k*D(u,x),x)` is rewritten (rule priority 10) into the
  conservative face-flux-difference stencil — interior cells
  (kf[i+½]·(u[i+1]−u[i]) − kf[i−½]·(u[i]−u[i−1]))/dx² with the **face**
  stiffness kf = 4·(ν·H_face + ε), ν from the face strain rate and
  H_face = (H_i+H_j)/2 (PISM's `compute_nuH_cfbc`, *not* the mean of cell
  values — the operator-form `k` is declared for the div(k·grad u) spelling
  but not read by the rule); the Dirichlet inflow row (u = V0 at the left
  wall, PISM's bc_mask identity row) and the calving-front row with the
  prescribed east-face flux δ_p (PISM's CFBC, `average_water_column_pressure`
  for P_water; reduces to ½·ρ_i·g·(1−ρ_i/ρ_w)·H² at the front). Mapping to
  SSAFDBase.cc is documented in the rule description for the
  EarthSciDiscretizations PR. Basal drag is structured for later addition
  (τ_b = 0 floating shelf; the β-row would enter the rule body).
- **Pseudo-transient** ✅ — the harness cannot solve the implicit/elliptic
  steady problem, so the SSAFD operator is relaxed:
  du/dt = U_char·(fluxdiv − ρ·g·H·h_x)/(ρ·g·H) with a penalty row
  du₁/dt = (V0 − u₁)/τ_dir (u₁ = V0 exactly at the fixed point) and the
  front row carrying the prescribed flux. Time is dimensionless; the
  slowest mode decays on τ̃ ≈ 3 pseudo-time units (fastest ~6e5), so
  t_end = 100 ≈ 35 decay times reaches the fixed point to machine precision
  (t ≈ 30). Verified sign: the dump's u closes this stencil to ~0.1% of the
  local driving stress (Picard residual), and the relaxation converges to the
  same fixed point.
- **Tests** ✅ — 26 assertions, passing under the default (native, BDF),
  `--compiler interpreter`, and `--solver sdirk`:
  `testV_warm_start_ic` pins the ic (the van der Veen exact solution
  u = Q0·(4C/Q0·X + H0⁻⁴)^¼ spelled from literals + N; rel 1e-9);
  `testV_steady_state_ssafd` asserts at t_end=100: u vs the dump at cells
  1/15/30/45/60 (rel 2e-3; cell 1 = V0 to rel 1e-6, the penalty row is
  exact), u vs the exact solution (rel 4e-3; SSAFD discretization error
  2.1–2.8e-3), taud_x vs the dump (rel 1e-10, bit-exact geometry),
  nuH_face vs the dump at faces 15/30/45/60 (rel 2e-3; the front face is a
  Schoof-regularized constant, rel 1e-3), flux_uH = Q0 (rel 4e-3), and the
  membrane residual fluxdiv − ρgH·h_x < 1e-3 Pa at cells 15/30/45/60
  (measured 2e-9…6e-7 Pa, solver-tolerance-limited).
- **Dump caveat (documented in the test description)** — PISM's Picard
  iteration stops on the relative nuH change < 1e-4
  (`stress_balance.ssa.fd.relative_convergence`), which leaves an equation
  residual of ~0.1% of the local driving stress; the conditioning of the
  SSAFD operator amplifies that to a ~1.5e-3 relative spread in u between
  the dump and the true fixed point. The assertions are calibrated to the
  fixed point (the dump's u closes the discrete equations to ~0.1%).
- **Harness findings (documented in the file comments)** —
  1. the native compiler refuses a `faq` in a build-time (initial-condition)
     scope that reads an observed produced by a rewrite rule (e.g. `fluxdiv`):
     the state-free analysis collects only the names in the body, classifies
     the faq as static, and refuses the per-cell walk it cannot take.
     `residual` and `taud_x` are therefore written as array-level elementwise
     expressions (skipped at build time, evaluated per step with the other
     observed rules);
  2. the native compiler refuses faqs that read an x-shaped array by a
     `face`-set symbol (cross-axis index map); `nuH_face` is a positional
     makearray with plain integer ranges;
  3. diffsol (the harness solver) integrates the pseudo-transient to
     t ≈ 108 (BDF, the default) / 290 (SDIRK) pseudo-time units and then
     fails in the quasi-static regime: its finite-difference-Jacobian Newton
     (step √ε·max(‖y‖,1)) cannot converge once the solution sits at the
     fixed point (the weighted Newton norm is dominated by the huge
     Schoof-regularized front-cell stiffness against the tight abstol).
     t_end = 100 stays short of the wall with deterministic margin;
  4. `./esm test` defaults to the BDF solver regardless of the document's
     `solver.stiffness` hint (the hint only selects tolerances), so the
     shipped test must pass under BDF — it does.
- **Verification** ✅ — `./esm validate` clean on both files; `./esm info`
  lists SSAFlowline (27 vars, 14 eqs); `./esm units --check`: 9 consistent /
  0 mismatched / 5 not checked (the same non-literal `^`-exponent
  informational pattern as ice_rheology/sia); regressions: ice_rheology.esm
  313/313 and sia.esm 48/48, native and interpreter.
- **C++ bugs** — none new; the SSA analysis found no additional issues in
  the instrumented code (bugs.md #1–#6 unchanged).

## Milestone 4 (continued) — Stage 2: PISM ice-energy (EnthalpyModel) component authored and passing ✅

The fourth Stage-2 deliverable: the PISM ice-energy boundary as a
hand-authored, compositional `.esm` component built from the Stage-1 energy
dumps (`stage1/dumps/energy/`, one 10-yr `EnthalpyModel` step of the test-F
exact-solution field). The scope is the **column solve** of
`enthSystemCtx`/`EnthalpyModel` (diffusion + vertical advection + explicit
horizontal advection + strain heating + surface Dirichlet + basal heat flux +
fine/coarse interpolation + basal melt rate), plus the DrainageCalculator
rate curve via the templates. TemperatureModel, BedThermalUnit, CHSystem and
the coupled 31×31 3-D field are **deferred** (documented below).

- **Deliverable** ✅ — `stage2/ice_energy.esm` (model EnthalpyColumn: 56
  vars, 19 eqs) importing the shared math by reference from
  `stage2/energy_templates.esm` (column-grid contract: metaparameters Mz=61 /
  Mz_fine=62, index sets `lev`/`lev_f`; templates `energy_thermal_diffusivity`,
  `energy_vertical_diffusion`, `energy_vertical_advection`,
  `energy_horizontal_advection`, `energy_coarse_to_fine`,
  `energy_fine_to_coarse`, `energy_drainage_rate`) and
  `stage2/ice_rheology_templates.esm` (enthalpy-converter free names). Both
  template libraries are first-time commits alongside the component.
- **Equations as math** ✅ — the state is the fine-grid enthalpy E_f on
  `lev_f`; the semi-discrete model is
  `dE_f/dt = D(K·D(E_f,z),z) − w·dE_f/dz + horizontal_advection + strain`,
  with the BOMBPROOF column stencils of enthSystemCtx::solve as rewrite
  rules (`energy_vertical_diffusion` rewrites the compound `D(K*D(E_f,z),z)`
  at priority 10; the basal row is the Neumann row of set_basal_neumann_bc
  with the geothermal gradient G, interior rows 2..ks the mean-face stencil,
  rows ks+1..0). The surface row and the atmosphere are pinned to the
  surface Dirichlet value Enth_ks in the initial condition; all four terms
  vanish there, so the pinned rows stay put under the harness integration.
  The cold/temperate K switch (assemble_R, K = k/(rho·c) or ratio·k/(rho·c))
  is assembled **once from the initial enthalpy** with the C++ basal special
  case K[0]←Enth[1] (the model writes it as a frozen faq over `E_f_init`).
- **Sign finding (vertical advection)** — the C++ tridiagonal is
  `D = I − dt·M`, so the vertical advection's contribution to the ODE
  `dE/dt = M·E + b` is the **negative** of the `energy_vertical_advection`
  stencil body (physical `−w·∂E/∂z`). The model writes
  `... − energy_vertical_advection(E_f, w_eff) + ...`; the wrong-sign ODE
  diverges from the dump by ~960 J/kg, the corrected one by ≤35 J/kg (dome
  base) / ≤1.5 J/kg elsewhere. The basal advection carries w[1]·A_b·G, which
  is zero because w[1] = 0 at the bed for the test columns (no basal melt).
- **C++-sign bug found in the head-start template (fixed)** —
  `energy_horizontal_advection` negated the v-term (through its `-` op) but
  **not** the u-term: it produced `+u·(E−E_w)/dx − v·(E−E_s)/dx` where the
  C++ RHS subtracts both upwind fluxes. The error was invisible for the
  zero-velocity columns and only surfaced in the shorter non-marginal columns
  with nonzero horizontal flow (c21_15, c26_15; c15_21 passed because its
  velocity is along v). Fixed by wrapping the sum of the two upwind fluxes in
  the leading negation: `hadv = −(u·(E−E_w)/dx + v·(E−E_s)/dx)`. Verified: the
  wrong-sign model reproduces the harness's failing values exactly
  (67550.9555 at c21_15 lev 10), the corrected one reproduces the dump to
  the ODE-integration gap.
- **Tests** ✅ — 104 assertions across 10 tests, passing under
  `--compiler interpreter` (BDF default and `--solver sdirk`): `ic_*` at
  t=0 (6 exact assertions each on the pinned fine-grid enthalpy and the
  storage-grid E_out) and `step_*` at t = dt for all 5 columns (E_out vs the
  dump at selected storage levels, the surface/atmosphere pinning, and the
  basal melt rate). Because the harness cannot force PISM's single
  implicit-Euler step, it integrates the semi-discrete operator and the
  assertions use abs tolerances = 2.5× the measured integration-scheme gap
  (dome base ≤35 J/kg, ≤1.5 J/kg elsewhere, c27_15 ≤9.25 J/kg near the
  surface; bmelt dump 1.0891e-10 m/s at the dome, 0 elsewhere).
- **Harness findings (documented in the file comments)** —
  1. the native compiler refuses the `K` faq ("const-cadence … unresolved
     symbol i"); tests run with `--compiler interpreter` (same documented
     harness limitation as ssa.esm);
  2. `coords` in assertions are **1-based** (like the `index` op): `lev` and
     `lev_f` coordinates are written +1 relative to the 0-based dump arrays,
     and the C++ 0-based m_c2f/m_f2c maps are baked in as 1-based parameters;
  3. booleans are not allowed in the schema, so `marginal` and the three
     `exclude_*` gates are 0/1 numbers;
  4. no document-level `index_sets` (would conflict with the template
     imports); `lev`/`lev_f` come from the templates' metaparameters;
  5. a match-less template (advection, f2c) is invoked via
     `apply_expression_template {name, bindings}`; the compound
     `energy_vertical_diffusion` fires through its match rule;
  6. the initial-condition equation `rhs` must be the bare parameter
     `"E_f_init"` (a faq RHS is rejected).
- **Verification** ✅ — `./esm validate` clean on both files; `./esm info`
  lists EnthalpyColumn (56 vars, 19 eqs); `./esm units --check`: 19
  consistent / 0 mismatched; regressions: ice_rheology.esm 313/313, sia.esm
  48/48, ssa.esm 26/26 (all `--compiler interpreter`).
- **C++ bugs** — none new; the horizontal-advection sign discrepancy was a
  template-authoring error on this side, not a C++ defect (bugs.md #1–#6
  unchanged).
- **Deferred (documented)** — TemperatureModel, BedThermalUnit, CHSystem
  (chemistry), the warm/floating basal branches of basal_melt_rate, the
  post-solve drainage/bulge-limit post-processing of DrainageCalculator
  (only the rate curve is tested), and the coupled 31×31 3-D field with
  horizontal exchange; these are follow-up components/subassemblies.

## Milestone 4 (continued) — Stage 2: PISM ice-age (AgeModel column) component authored and passing ✅

The fifth Stage-2 deliverable: the PISM ice-age boundary as a
hand-authored, compositional `.esm` component built from the Stage-1 age
dumps (`stage1/dumps/age/`, one dt = 315569259.747 s step from age = 0 of the
test-F exact-solution velocity field). The scope is the **column solve** of
`AgeColumnSystem`/`AgeModel` (age advection PDE, surface and basal-accumulation
BCs, fine-grid storage interpolation). Coupling to real velocities
(SIA/SSA), the full 3-D multi-column solve, the `AgeModel::update`
negative-age clamp, and the ks-from-thickness computation are **deferred**
(documentated in the component).

- **Deliverable** ✅ — `stage2/ice_age.esm` (model AgeColumn: 62-level fine
  grid `lev` (dz = Lz/61 = 65.5738 m — the C++ `ceil(Lz/min_dz)+1`
  FP-jittered construction), 61-level storage grid `lev_s`, 15 vars, 4 eqs,
  4 tests / 24 assertions) importing the shared math by reference from
  `stage2/age_templates.esm` (`age_vertical_advection`,
  `age_horizontal_advection`, `age_storage_interpolate`).
- **Equations as math** ✅ — the age PDE written as
  `D(tau, wrt:t) = 1 − hadv − w·dτ/dz`: first-order upwind advection,
  implicit in the vertical (an ODE term integrated by the harness's stiff
  solver, mirroring the C++ backward-Euler column step) and horizontal at old
  time from neighbor-column age parameters (explicit in the C++); BCs:
  above-surface and surface cells pinned to age 0, basal-accumulation pin
  when `w3[1] > 0`, ice-free columns stay 0.
- **Tests** ✅ — 24/24 passing under native BDF, SDIRK, and the interpreter
  compiler. Measured agreement vs the dump (`stage1/dumps/age/columns.csv`):
  deep/mid levels of the dome (15,15), (11,3) and (6,11) columns within
  ≤9e-5–1.7e-4 relative (asserted rel 1e-3); the two near-surface dome levels
  4.4e-3/6.2e-3 (asserted rel 1e-2); ice-free (0,0) exactly 0. The
  exact-ODE-vs-dump gap is the difference between continuous ODE integration
  and the C++ single backward-Euler step; the agent verified its transcription
  of the C++ tridiagonal reproduces the dump to 1e-14–1e-16.
- **Verification** ✅ — `validate`/`info`/`units --check` clean; full suite
  571/571 native (rheology 313, sia 48, ssa 26, basal 160, age 24).
- **C++ bugs** — none new (bugs.md #1–#6 unchanged); the 62-level fine grid
  from storage-spacing FP jitter is a documented construction detail, not a
  defect.

## Milestone 4 (continued) — Stage 2: PISM basal-strength component authored and passing ✅

The sixth Stage-2 deliverable: the PISM basal-strength boundary as a
hand-authored, compositional `.esm` component built from the Stage-1
basal-strength dumps (`stage1/dumps/basal_strength/`: Mohr-Coulomb pointwise
sample grids, three drag-law grids, and the tauc map on the test-F exact
geometry). All functions are **D** (pointwise, instantaneous).
`OptTillphiYieldStress` (inverse friction-angle optimization, **I**) and the
SSA consumption of β into the SSA matrix (coupling) are **deferred**.

- **Deliverable** ✅ — `stage2/basal_resistance.esm` (6 models:
  `MohrCoulombPointwise` (effective_pressure / yield_stress /
  till_friction_angle), `MohrCoulombYieldStress` (tauc map with the ice-free
  high_tauc branch, P = ρ·g·H), `ConstantYieldStress`, and the three
  `IceBasalResistance` drag laws — plastic (PISM's default SSA law),
  pseudo-plastic, regularized — each its own model so the SSA component can
  import the configured law by reference; 62 tests / 160 assertions) importing
  the shared math by reference from `stage2/basal_templates.esm` (13
  templates: effective pressure, yield stress / inverse friction angle,
  overburden, basal speed machinery shared by all three drag laws).
- **Equations as math** ✅ — Tulaczyk et al. 2000 effective pressure
  `N_interior = N0·(δP/N0)^s·10^((e0/Cc)(1−s))`, `s = W/W_max`,
  `N_till = min(P, N_interior)` (W_max = 2 m, N0 = 1000 Pa, e0 = 0.69,
  Cc = 0.12, δ = 0.02); yield stress `tauc = c0 + tan(φ)·N_till` (c0 = 0) and
  its exact inverse; the three drag laws with analytic `dbeta` (w.r.t.
  α = 0.5|v|²), q = 0.25, u_th = 100 m/yr (converted with the UDUNITS
  tropical year 3.15569259747e7 s — the key constant needed to reproduce the
  dumps to full double precision), ε = 0.01 m/yr plastic regularization.
- **Tests** ✅ — 160/160 passing at rel 1e-12 (bit-exact to dump precision):
  φ→tauc→φ round-trips exact (20 → 20.000000000000004), effective-pressure
  clamping knee, all three drag laws with `dbeta` cross-checked against the
  dump's central-difference `dbeta_fd` column (rel err ~1e-10 in the dump),
  and the tauc map at selected (i,j) cells. Non-finite dump rows
  (`N_interior = inf` at δP = 0, W < 0; `phi_inverse = NaN` at N_till = 0)
  reproduced and handled in the tests.
- **Verification** ✅ — `validate`/`info`/`units --check` clean; full suite
  571/571 native.
- **C++ bugs** — none; no discrepancies found between any dumped value and
  the C++ formulas (bugs.md #1–#6 unchanged).

## Milestone 4 (continued) — Stage 2: PISM bed-deformation (PointwiseIsostasy) component authored and passing ✅

The seventh Stage-2 deliverable: the PISM bed-deformation boundary as a
hand-authored, compositional `.esm` component built from the Stage-1 bed
dumps (`stage1/dumps/bed/`, two full update intervals of
`bed::PointwiseIsostasy`, `bed_deformation.model = "iso"`, on the test-F
exact geometry). `Given`, `Null`, and `LingleClark`/`LingleClarkSerial`
(Fourier spectral collocation) are **deferred**.

- **Deliverable** ✅ — `stage2/bed_deformation.esm` (4 models:
  `ComputeLoad`, `PointwiseIsostasyUpdate`, `BedDefLoadAccumulator`,
  `PointwiseIsostasyTwoIntervals`; 131 tests / 346 assertions, all at
  rel 1e-12) importing the shared math by reference from
  `stage2/bed_templates.esm` (5 templates: isostatic compensation factor
  f = ρ_ice/ρ_mantle, bed ocean depth, ocean load (ρ_w/ρ_i)·max(sea_level−bed,0),
  ice load (floating ice excluded), isostatic bed update).
- **Equations as math** ✅ — `bed_out = bed_in − f·(load − load_last)`,
  `load_last ← load`, with `load_avg = (1/dt_beddef)·Σ load·dt` over the
  update interval (the C++ accumulator divides by dt_beddef at the update
  time); `init_impl` sets load_last to the initial load (H0), so interval 1
  responds to a +150 m load change and interval 2 to +300 m.
- **Tests** ✅ — all 120 `pointwise.csv` rows (compute_load + update_law)
  asserted at rel 1e-12; the two-interval `columns.csv` trace at 7 cells
  (dome, near-dome, mid-radius, near-margin, 3 ice-free) asserting
  bed_out_1/2, uplift_1/2, load_last_1/2; the accumulator divide-by-dt_beddef
  subtlety including the (5,7) cell where FP round-off shifts load_avg_1 by
  ~6e-14 m (reproduced). A Python replication of the C++ FP chain reproduces
  every dump value bit-exactly for all 517 ice cells; the .esm chain matches
  to ~1e-15 rel (uplift within 1–2 ulp — C++ computes
  (topg−topg_last)·(1/dt_beddef), the .esm divides).
- **Harness note** — the event-based state-machine form (discrete events at
  update times) is not usable in the inline-test pathway
  (`UnsupportedConstruct` for `discrete_events`, esm-spec §9.6.6); the
  accumulator + two-interval trace models cover the same semantics
  (documented in the component).
- **Verification** ✅ — `validate`/`info`/`units --check` clean (14/14
  consistent); full native suite 1390/1390.
- **C++ bugs** — none; `bed_analytic_1/2` equal `bed_out_1/2` in the dump
  (self-consistent) (bugs.md #1–#6 unchanged).

## Milestone 4 (continued) — Stage 2: PISM surface/ocean forcing components authored and passing ✅

The eighth Stage-2 deliverable: the PISM coupler boundaries as hand-authored,
compositional `.esm` components built from the Stage-1 surface/ocean traces
(`stage1/dumps/surface_ocean/`, instantaneous-derivative type: PIK
atmosphere, PIK + Constant ocean, PIK surface partition, water-column
pressure, yearly-cycle time series). All required pieces are **D**
(pointwise, instantaneous). The PDD surface model (**I**), Given/Uniform
atmosphere and ocean models, Anomaly/Delta_* modifiers, ISMIP6/DEBM, PICO,
GivenTH (Holland–Jenkins), ForceThickness, sea-level coupling, and
frontal-melt physics are **deferred** (documented in the components).

- **Deliverable** ✅ — `stage2/surface_forcing.esm` (AtmospherePIK,
  SurfacePIK, SurfaceOceanCoupler per-cell wiring), `stage2/ocean_coupling.esm`
  (OceanPIK, OceanConstant), and `stage2/surface_forcing_templates.esm`
  (16 expression templates: martin mean-annual temperature, yearly cycle,
  prescribed precipitation, shelf base temperature (pressure melting point,
  shared PIK+Constant), Beckmann–Goosse freezing temperature, ocean-to-ice
  heat flux, shelf base mass flux, water-column height/pressure, constant
  ocean mass flux, SMB partition + mass-flux identity). 473/473 assertions
  (OceanPIK 221, OceanConstant 6, AtmospherePIK 128, SurfacePIK 36,
  SurfaceOceanCoupler 114), all at rel 1e-12.
- **Equations as math** ✅ — `T_ma = 303.15 + 0.68775·lat − 0.0075·usurf`
  (K), `T_ms = T_ma`, flat yearly cycle; `P(r) = 100·(2 + cos(πr/Lx))`
  kg m⁻² yr⁻¹; `T_shelf = 273.15 − β_CC·ρ_i·g·H` (linear pressure melting
  point, decreasing with depth); `z_b = −(ρ_i/ρ_w)·H`,
  `T_f = 273.15 + 0.0939 − 0.057·S + 7.64e-4·z_b`,
  `Q = melt_factor·ρ_w·c_p·γ_T·(T_ocean − T_f)`, `mf = Q/L`
  (melt_factor 0.005, ρ_w 1028, c_p 3974, γ_T 1e-4, S 35, T_ocean −1.7 °C,
  L 334000); `awcp = ½ρ_w·g·h_w²/H`; Constant ocean `mf = melt_rate·ρ_i`
  (0.051914193590840291 m/yr). SMB partition: accumulation = max(smb,0),
  melt = runoff = max(−smb,0), mass_flux = accumulation − runoff.
- **Tests** ✅ — all pointwise CSVs (shelf_base_temperature 31 rows,
  beckmann_goosse mass-flux chain 31×4, water_column_pressure 20×3,
  atmosphere temperature/precipitation curves, surface partition) at
  rel 1e-12; the per-cell `columns.csv` maps (incl. the `_const` Constant
  ocean cross-check) with the analytic reference columns (T_ma_analytic,
  T_shelf_analytic, mf_analytic, P_analytic) asserted to coincide (manager
  re-verified: 961/961 cells each for T_ma and T_shelf).
- **Verification** ✅ — `validate`/`info`/`units --check` clean (22+10
  consistent); full native suite 1390/1390.
- **C++ notes** — `ocean/ConstantPIK.cc` mass-flux comment is wrong
  (already tracked as bugs.md #5; reproduced as positive = melting); the
  mf map vs mf_analytic differ by ≤2 ULP at a few cells (expression order);
  both PIK and Constant oceans produce nonzero mass flux in ice-free cells
  (real C++ behavior, reproduced). No new bugs (bugs.md #1–#6 unchanged).

## Milestone 4 (continued) — Stage 2: PISM hydrology (Routing/Shreve) component authored and passing ✅

The ninth Stage-2 deliverable: the PISM subglacial-hydrology boundary as a
hand-authored, compositional `.esm` component built from the Stage-1
hydrology traces (`stage1/dumps/hydrology/`, ONE 0.25-yr update step of
`hydrology::Routing`, Shreve q = −K∇ψ, on the test-F exact geometry with the
dome-centered water blob). `Distributed`, `NullTransport`, and
`SteadyState`/`EmptyingProblem` are **deferred**.

- **Deliverable** ✅ — `stage2/hydrology.esm` (the Routing/Shreve model as a
  continuous PDE over the 31×31 field; 356 assertions) importing the shared
  math by reference from `stage2/hydrology_templates.esm` (180 lines), with
  test-generation tooling committed under `stage2/tools/`
  (`build_hydrology_tests.py`, `fill_hydrology_data.py`, and
  `instrument_hydrology_substeps.cc` — a Stage-2 analysis driver that
  replicates `Routing::update_impl`'s internal substep loop with dump hooks,
  self-validated bit-for-bit against the Stage-1 dump, exposing the
  per-substep trajectory the one-step Stage-1 dump cannot show).
- **Equations as math** ✅ — the Shreve flux law q = −K∇ψ with the
  conductivity/flux law K(W) and the W / W_till storage split, forced and
  initialized exactly as in the dump (W = 0.5·exp(−(r/150 km)²) blob,
  dry till, surface input 0.2 m/yr we, basal melt blob 0.1 m/yr). No `ic`
  equations: the native compiler's build-time initial-condition scope cannot
  evaluate observed definitions that read the state (the same harness quirk
  documented in ssa.esm); the initial state is declared as unknown defaults.
- **Tests** ✅ — 356/356 on the native compiler: the staggered edge operators
  (Wstag/Kstag/Vstag/Qstag) at t = 0 and at the last-substep state vs
  `staggered.csv`; the initial state (W blob, dry till, overburden) exact;
  the flux law at all 99 `flux_law.csv` samples (rel 1e-12); the one-step
  evolution (W_out, W_till_out, time-averaged flux, P_over) at documented
  Euler-vs-ODE scheme-gap tolerances: W/flux rel 0.05 (measured max 4.4% —
  the C++ takes one discrete backward-Euler-ish substep, the .esm integrates
  the continuous PDE), W_till rel 1e-6, overburden rel 1e-12.
- **Verification** ✅ — `validate`/`info`/`units --check` clean; full native
  suite 3167/3167.
- **C++ bugs** — none new (bugs.md #1–#6 unchanged).
- **Suite note** — the worker's follow-up commit (7210cba, made after the
  merge) trimmed the suite to 316 assertions by removing the
  `staggered_operators_last_substep` test; main keeps the 356-assertion
  version, verified green under native and SDIRK (the interpreter compiler
  times out on this large suite — a harness slowness, not a failure). The
  removed test pins the staggered operators at the last-substep state via
  the `initial_conditions` seeding machinery.

## Milestone 4 (continued) — Stage 2: PISM calving/front-retreat component authored and passing ✅

The tenth Stage-2 deliverable: the PISM calving / front-retreat boundary as a
hand-authored, compositional `.esm` component built from the Stage-1 calving
traces (`stage1/dumps/calving/`, on the van der Veen CFBC shelf state, square
65×65 grid). PrescribedRetreat (**I**), DischargeGiven/DischargeRouting
frontal melt, FrontRetreat::update_geometry part-grid, and the
`icemodel/frontretreat.cc` orchestration are **deferred** (documented).

- **Deliverable** ✅ — `stage2/calving.esm` (8 models: PrincipalStrainRates,
  EigenCalving, vonMisesCalving, HayhurstCalving, CalvingAtThickness,
  FloatKill, FrontalMeltUndercutting, FrontalMeltISMIP6, ConstantFrontalMelt;
  812 tests / 1421 assertions — 72 hand-authored + 1349 dump-derived)
  importing the shared math by reference from `stage2/calving_templates.esm`
  (17 expression templates).
- **Equations as math** ✅ — verified against the dumps at rel 0
  (bit-exact): eigen-calving rate = K·eigen1_avg·eigen2_avg at ice-free-ocean
  cells next to floating ice (K = 3e16 m·s), with the principal strain rates
  from the centered-difference reduction of
  `compute_2D_principal_strain_rates`; von Mises: e_s, σ̃ = √3·B·e_s^(1/n),
  rate = |v|_avg·σ̃/σ_max (B = 1.9e8 Pa s^(1/3), n = 3, σ_max = 1e6 Pa);
  Hayhurst: B̃·(1e-6)^r/yr·(1−ω^2.8)·(σ₀−σ_th)^r·H_eff with the floating
  H_eff adjustment (B̃ = 65, r = 0.43, σ_th = 170 kPa); CalvingAtThickness
  (H < 300 m at the front) and FloatKill in-place removal; the
  FrontalMeltPhysics kernels q_m = (A·h·q_sg^α + B)·TF^β (full 294/210-row
  grids); ConstantFrontalMelt (1 m/day melt, retreat 0.8852 m/day with
  include_floating_ice, H_submerged = (ρ_i/ρ_w)·H_threshold for floating
  fronts).
- **Tests** ✅ — 1421/1421 on the native compiler; every law asserted on
  both its branches (zero and nonzero rates; both include_floating_ice
  settings). Unit-checker m/yr and m/day conversions handled by declaring
  `seconds_per_year` (s/yr) and `seconds_per_day` (s/day).
- **Dump findings** (documented in the component) — the Hayhurst
  `sigma_0_Pa` dump column is the instrument's reconstruction with the
  *unadjusted* H (it does not equal the component's internal σ₀; the rate,
  which matches exactly, uses the H_eff-based σ₀); the (64,·) ghost column is
  a periodic copy of column 0 (excluded from tests). No C++ bugs
  (bugs.md #1–#6 unchanged).
- **Verification** ✅ — `validate`/`info`/`units --check` clean; full native
  suite 3167/3167.
- **Deferred (documented)** — front-cell sampling averages (+/−2-cell
  stencil), the part-grid threshold-thickness stencil, next-to-* mask
  properties, Hayhurst ice-free-cell neighbor propagation, margin one-sided
  strain-rate stencils.

## Milestone 4 (continued) — Stage 2: PISM geometry (mass continuity) component authored and passing ✅

The eleventh Stage-2 deliverable: the PISM geometry boundary as a
hand-authored, compositional `.esm` component built from the Stage-1
geometry traces (`stage1/dumps/geometry/`, ONE dt = 1 yr step each of
`GeometryEvolution::flow_step` and `source_term_step` from the same test-F
exact geometry, plus the `part_grid_threshold_thickness` pointwise sample).
With this commit, **all 12 Stage-2 leaf components are authored and passing**
(see PLAN.md §10.4); the next phase is the subassemblies.

- **Deliverable** ✅ — `stage2/geometry.esm` (3 models: flow_step,
  source_term_step, part-grid threshold; 81 assertions — 27 flow, 32
  source, 22 part-grid) importing the shared math by reference from
  `stage2/geometry_templates.esm` (347 lines: face-flux construction, mask
  limiter, divergence rules, effective_change, cell mask/surface,
  part_grid_threshold_thickness).
- **Equations as math** ✅ — `flow_step` as the mass-continuity PDE
  dH/dt = −div(Q) with mask-limited diffusive + upwind advective face fluxes
  and the staggered finite-volume divergence on the zero-flux closed domain
  (the outermost ring is ice-free so sum(div Q) telescopes to 0 — reproduced:
  conservation integral residual −3.3e-12 m); `source_term_step` as the
  pointwise effective-change map with non-negativity clamping
  (dH_SMB = dt·smb/ρ_ice, dH_BMB = −dt·bmr); `part_grid_threshold_thickness`
  pointwise (22 cases, machine-exact rel 1e-12).
- **Tests** ✅ — 81/81 on the native compiler: flow assertions at
  dome/margin/ice-free cells carry the measured forward-Euler-vs-exact-ODE
  scheme gap (H/surface rel ≤5e-6, dH_flow rel ≤3e-4, documented); the
  closed-domain mass-balance integral and non-negativity min(H) also
  asserted; the algebraic source step and the part-grid thresholds are
  machine-exact (rel 1e-12).
- **Harness finding (worth an upstream EarthSciAST issue)** — the `ic`
  equation had to be dropped: for ANY multi-model document the inline-test
  build path (`ArrayCompiled::from_flattened`) cannot resolve a
  scoped-reference `ic` whose RHS is a document-inline array (the flattener
  routes `ic` equations into `flat.field_ics`, so the resolver sees no `ic`
  in the synthetic model). Single-model docs and provider-loaded fields
  work. Fix: the initial profile is declared as an inline-array `default`
  (spec-sanctioned, `StateDefault::Field`). Recommended compiler fix (small):
  in `from_flattened`, append `flat.field_ics` as `ic`-LHS equations before
  `from_model`. This matters for the subassembly phase (coupled models).
- **C++ bugs** — none found in `GeometryEvolution.cc` (bugs.md #1–#6
  unchanged).
- **Verification** ✅ — `validate`/`info`/`units --check` clean; full
  native suite 3248/3248.

## Milestone 4 (continued) — Stage 2: stress-balance subassembly authored and passing ✅

The first TRUE cross-component assembly — `stage2/stress_balance.esm` (merged
`<merge>`, from wt/sb commit `36d177e` (merged as part of the subassembly milestone)): the PISM `StressBalance` assembler with
the SIA and SSA leaf components **mounted by reference** and **coupled**, plus
the three shared D-methods, with scenario tests vs the sia/ssa dumps.
This proves the composition machinery end-to-end.

- **Assembly wiring** — 8 top-level `models.<k>` `{ref}` mounts with `"model"`
  selectors (`SIADiffusivityColumn` ×4 (I columns), `SIADiffusiveFlux`,
  `SIA3DHorizontalVelocity`, `SIAVerticalVelocity`, `SSAFlowline`) and 14
  `variable_map` couplings (incl. `{"op":"index",...}` literal-index
  transforms and a `param_to_var` hybrid-sliding coupling
  `StressBalance.ubar/vbar` → `SIAVel.sliding_velocity_*`).
- **Harness capability (empirically established)** — **`esm test` RUNS a
  coupled document**: mounted tests are dropped (§6.6 normative), the
  subassembly's own 47 assertions run under the native compiler (and
  interpreter). Two native constraints shaped the equations: differing-branch
  `ifelse` boxes are refused → ConstantInColumn selection via lev-shaped
  `ubar_col/vbar_col` boxes; index maps must be affine/constant → the
  u_z/v_z bed/interior switch is a `makearray` region split (algebraically
  identical).
- **Assembler semantics (verified against C++)** — PISM 2.3.2 has NO per-cell
  mask-based SIA/SSA selection: `stress_balance.model` picks ZeroSliding+SIAFD
  (`sia`), SSAFD+ConstantInColumn (`ssa`), or SSAFD+SIAFD (`ssa+sia`) at the
  config level, and SIAFD applies everywhere (no grounded gate). Encoded as
  `sb_model` (0/1/2). Σ = 2·e^(−1/n)·B·D2^(0.5(1/n+1)) with B from the
  pressure-adjusted temperature (the C++'s T_pa-vs-absolute-T inconsistency,
  already in bugs.md, confirmed against dumps); w = trapezoid of −div via the
  mounted SIAVerticalVelocity (w[0]=−bmr, bmr=0 in traces); F_b = C·|ubar|²,
  plastic C = tauc/√(ε²+…), skipped on mask ocean (3) AND ice-free (4).
- **Tests** — 47/47 native + interpreter: sia grounded (20,15) 20/20
  (u3/v3/w3/sigma vs columns.csv, q_x vs diffusive_flux.csv, machine
  precision), sia margin (27,15) 15/15 (one-sided x-stencil), ssa floating
  cell 29 (ubar = SSA steady state rel 2e-3; hybrid ubar = u_exact(29) =
  2.5762120147478445e-05 bit-exact, verified by the manager vs the dump),
  F_b synthetic (grounded 100.03317883754333 W/m², ocean/floating 0).
- **Leaf fidelity bug fixed (found by the subassembly)** — `stage2/sia.esm`
  I(z): the cumulative trapezoid added a spurious partial top interval above
  each column's surface; `SIAFD::compute_I` fills levels above the ice with an
  exact copy of the top in-ice value. Fix: `j <= ks` cap (8-line diff); leaf
  tests unchanged (48/48). This is an .esm fidelity bug, not C++.
- **C++ bugs** — none new (T_pa inconsistency re-confirmed; bugs.md #1–#6).
- **Deferred to Stage 3** — full-grid neighbor couplings (dump-sampled
  parameters here), model-1 div/w3 + model-2 hybrid full-field assertions,
  the −bmr base term, BedSmoother, all four flux faces,
  `vertical_velocity_approximation=upstream`, full 2-D SSA hybrid with basal
  drag.
- **Verification** — validate/info/units clean; full suite 3295/3295 native +
  104 interpreter.

## Milestone 4 (continued) — Stage 2: 3-D coupled energy and age subassemblies authored and passing ✅

The deferred 3-D pieces of the energy and age boundaries, merged from wt/energy3d
(commit `84cc203` → merge `404a20f`) and wt/age3d (commit `2896e62` → merge
`1e2f2e5`). Both are full 3-D PDE solves with **live horizontal exchange**
between columns — the horizontal-advection terms the leaf column models
deferred.

- **Energy 3-D** — `stage2/ice_energy_3d.esm` (91 assertions) +
  `energy_3d_templates.esm` (3-D template lifts) + `tools/build_energy3d_data.py`:
  the full enthalpy PDE dE/dt = ∂/∂z(K ∂E/∂z) − w·∂E/∂z − (u·∂E/∂x + v·∂E/∂y)
  + σ/ρ over four 2×3-column sub-domains, built on the validated EnthalpyColumn
  (BOMBPROOF vertical stencils, ColumnInterpolation maps, cold/temperate
  diffusivity switch, basal melt) imported/lifted by reference. 91/91 under the
  interpreter (23m44s). The last open item was NOT a model bug: (16,15) is the
  only cell whose second fine layer starts at the cold-temperate transition at
  the last bit (E_f[2] = E_cts[1] + 4e-11), switching the basal diffusivity
  entries R[0]/R[1] to the temperate ratio·K (ratio=0.1) — a ~392 J/kg
  one-step scheme gap vs the dump, measured and documented (tolerances widened
  to the measured values; sign of the gap flips with local dynamics). The
  **no-coupling demonstration** assertion: a hadv-free control copy fails at
  (7,15) lev 21 by 1.320 J/kg (> 0.1 tol), proving the horizontal exchange is
  live.
- **Age 3-D** — `stage2/ice_age_3d.esm` (model AgeBlock3D, 62 assertions) +
  `age3d_templates.esm` + `tools/build_age3d_data.py`: the age PDE
  ∂τ/∂t + u·∂τ/∂x + v·∂τ/∂y + w·∂τ/∂z = 1 on a 7×7-column block (~9k unknowns,
  ~0.7 GB dense Jacobian), 62-level fine grid, **native** compiler, 62/62.
  Block-edge ghosts frozen at the full-field initial age (exactly zero here).
- **Harness capability fact (answers the 3-D scaling question)** — the ODE
  harness's implicit solvers (BDF/SDIRK, the default) form a **dense
  finite-difference Jacobian** over the whole state vector: the full 31×31
  field (~119k unknowns) needs ~113 GB and was measured and killed in an
  early probe; the 7×7 block (~9k unknowns) needs ~0.7 GB and integrates
  fine. The documented sub-domain pattern is therefore the established route
  for implicit 3-D fields (energy uses four 2×3 blocks, age one 7×7 block).
  **Explicit ERK is available (`--solver erk`) and needs NO Jacobian (O(N)
  memory)** — manager probes: geometry 81/81 under ERK; age 3-D passes (only
  two atol=0 assertions on exactly-zero quantities differ at the 1e-21 level);
  the SSA relaxation reaches only ~1e-3 residual under ERK (vs machine-zero
  BDF) — the elliptic equilibrium is the one piece that wants the implicit
  solver, but its state is 2-D (~2k unknowns, ~30 MB dense Jacobian). The
  113 GB figure applies only to ONE monolithic BDF solve over the whole 3-D
  field — the Stage-3 orchestration should use per-substep solvers (ERK for
  transport on the full grid, BDF for the small stiff pieces).
- **C++ bugs** — none found (energy3d worker: "the only bug found was in the
  previous attempt's own build script — a transposed block_slice"; age3d:
  none).
- **Verification** — validate/info/units clean on both; age3d 62/62 native,
  energy3d 91/91 interpreter (both re-run by the manager); full 13-file native
  regression 3443/3443 per the energy3d worker (re-run on main, see below).

## Milestone 4 (continued) — Stage 2: subassembly phase launched (3 workers in flight)

With all 12 leaf components merged (full native suite 3248/3248 + energy 104
under interpreter, HEAD `a8d5318`), the subassembly phase is underway. The
manager probed the harness's composition machinery first: top-level
`models.<k>` `{ref}` mounts with a `model` selector and `variable_map`
`coupling` entries validate and load in the `esm` CLI (probe mounting
`sia.esm`'s `SIA3DHorizontalVelocity` + `{"type": "variable_map", "from":
"SIA.u3", ...}` passed; coupling entry types per spec §10:
`operator_compose` / `couple` / `variable_map`). Three workers dispatched in
parallel worktrees:

1. **stress_balance subassembly** (`wt-sb`) — `stage2/stress_balance.esm`:
   mounts the SIA and SSA components by reference, couples their outputs, and
   implements the shared StressBalance D-methods (volumetric strain heating
   Σ, vertical velocity w(z), basal frictional heating F_b) with
   scenario-based tests vs the sia/ssa dumps. First true cross-component
   assembly; mounted tests are dropped at the mount (normative), so it
   authors its own tests. Open question being resolved empirically: can
   `esm test` run a coupled document?
2. **Energy 3-D coupled field** (`wt-energy3d`) — completes the deferred
   full-field 31×31×61 enthalpy solve with horizontal exchange, built on the
   single-column EnthalpyColumn.
3. **Age 3-D multi-column solve** (`wt-age3d`) — completes the deferred
   full-field age solve with horizontal exchange, built on the single-column
   AgeColumn.

Next after these: top-level `eqice.esm` (Stage 3) needs full-model C++
simulations with component-boundary dumps (new instrumentation; PLAN.md §6).

## Stage 2 — Stubs → physics → EarthSciModels PRs (in progress)

Per PLAN.md §5: hand-author stub `.esm` files with `tests` blocks from
Stage-1 tuples, fill in physics, review + merge into EarthSciModels.
Rheology, SIA/SSA, and the energy column (milestones 4 above) are done; the
remaining boundaries follow the PLAN.md §5 sequencing (basal/hydrology →
bed → surface/ocean/calving → subassemblies).

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
