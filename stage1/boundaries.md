# Stage 1 — PISM Subassembly Boundary Catalog

This is the deliverable of PLAN.md §4 Stage 1, step 1: *identify discrete
subassembly boundaries in the C++ codes matching the component inventory*.
It maps each row of the PLAN.md §3 component inventory to concrete PISM 2.3.2
C++ classes, methods, and input/output boundaries, so that Stage-1
instrumentation (step 2) and Stage-2 `.esm` stubs (§5) have an authoritative
target list.

Reference implementation: PISM **v2.3.2** (`stage1/pism`, git tag `v2.3.2`,
commit `fa11747`). All paths below are relative to `stage1/pism/src/`.

Legend: **D** = instantaneous/derivative computation (pointwise function of
current state — preferred for tests); **I** = integrated trajectory (advances
state over `dt`).

## 0. Top-level call order (`IceModel::step`, `icemodel/IceModel.cc:398`)

1. `m_stress_balance->update(stress_balance_inputs(), updateAtDepth)` —
   computes u,v,w, strain heating Σ, basal frictional heating F_b
2. `max_timestep(...)`
3. `m_basal_yield_stress_model->update(yield_stress_inputs(), t, dt)` → tauc
4. `m_age_model->update(...)` (if enabled)
5. `energy_step(t, dt)` — `bedrock_thermal_model_step` then
   `m_energy_model->update(t, dt, energy_model_inputs())`
6. mass transport (`flow_step`, `source_term_step`), front retreat, hydrology

Couplings: energy/age consume u,v,w,Σ,F_b from step 2; tauc from step 3 is
used by the *next* step's stress balance. Couplers (atmosphere/ocean/surface)
are updated inside `step()` **after** the energy step (see bugs.md #3).

## 1. Rheology — `rheology/` → `ice_rheology.esm`

| Class | File | Config name |
|---|---|---|
| `FlowLaw` (abstract) | `rheology/FlowLaw.hh` | — |
| `PatersonBudd` | `rheology/PatersonBudd.hh` | `"pb"` |
| `IsothermalGlen` | `rheology/IsothermalGlen.hh` | `"isothermal_glen"` |
| `PatersonBuddCold` / `PatersonBuddWarm` | `rheology/PatersonBuddCold.hh` / `PatersonBuddWarm.hh` | `"arr"` / `"arrwarm"` |
| `Hooke` | `rheology/Hooke.hh` | `"hooke"` |
| `GoldsbyKohlstedt` | `rheology/GoldsbyKohlstedt.hh` | `"gk"` (forward-form only; excluded from SSA) |
| `GPBLD` | `rheology/GPBLD.hh` | `"gpbld"` |
| `FlowLawFactory` | `rheology/FlowLawFactory.hh` | registry |

Key boundary methods (namespace `pism::rheology`, all `const`):

- `double softness_impl(double E, double p) const = 0` — rate factor A (softness).
  **D**. Inputs: enthalpy E (J/kg), pressure p (Pa).
- `double hardness_impl(double E, double p) const` — hardness B = A^(−1/n). **D**.
- `double flow_impl(double stress, double E, double pressure, double grainsize) const`
  — forward form ε̇ = A·σ^(n−1). **D**.
- `void effective_viscosity(double hardness, double gamma, double *nu, double *dnu) const`
  — ν = ½·B·(ε+γ)^((1−n)/(2n)). **D**.
- `double averaged_hardness(const FlowLaw&, double H, unsigned kbelowH, const double *z, const double *E)`
  — vertical trapezoid integral of B̄(E,p). **D** (spatial integral).
- `double secondInvariant_2D(const Vector2d &U_x, const Vector2d &U_y)` — γ. **D**.

Enthalpy→temperature conversion: `util/EnthalpyConverter.hh`
(`temperature(E,p)`, `pressure_adjusted_temperature(E,p)`, `water_fraction(E,p)`,
`enthalpy_cts(p)`, `pressure(depth)`). **D**.

Enhancement factor is applied **outside** the flow law, at three call sites:
SIA diffusivity (`SIAFD::compute_diffusivity`), SSA νH (`SSAFDBase::compute_nuH_*`),
strain heating (`StressBalance::compute_volumetric_strain_heating`).

## 2. Stress balance — `stressbalance/` → `sia.esm`, `ssa.esm`, `stress_balance.esm`

### 2.1 Top-level assembler — `stressbalance/StressBalance.{hh,cc}`

```cpp
void StressBalance::update(const Inputs &inputs, bool full_update);
```

`Inputs` (StressBalance.hh:40): `geometry` (bed, thickness, surface, cell type,
grounded fraction, sea level), `basal_melt_rate`, `basal_yield_stress` (tauc),
`water_column_pressure`, `fracture_density`, `enthalpy`, `age`, `bc_mask`,
`bc_values`, regional no-model fields. `Inputs::dump(filename)` already writes
all inputs to NetCDF — a natural input-instrumentation hook.

Outputs (getters): `advective_velocity()` (2D), `diffusive_flux()` (staggered),
`max_diffusivity()`, `velocity_u/v/w()` (3D), `basal_frictional_heating()`,
`volumetric_strain_heating()`, `max_timestep_cfl_2d/3d()`.

Key internal methods (all **D**):
- `compute_volumetric_strain_heating(const Inputs&)` — Σ = 2·e^(−1/n)·B·D^(1/n+1).
- `compute_vertical_velocity(mask, u, v, bmr, result)` — w(z) = −∫(u_x+v_y)dζ + w_b.
- `ShallowStressBalance::compute_basal_frictional_heating(V, tauc, mask, result)` — F_b = C·|V|².

### 2.2 SIA — `stressbalance/sia/SIAFD.{hh,cc}` (class `SIAFD : public SSB_Modifier`)

```cpp
void SIAFD::update(const array::Vector &sliding_velocity, const Inputs &inputs, bool full_update);
```

Sub-methods (each a clean **D** boundary):
- `compute_surface_gradient(inputs, h_x, h_y)` — ∇h (eta/haseloff/mahaffy methods).
- `compute_diffusivity(full_update, geometry, enthalpy, age, h_x, h_y, result)` —
  **the flow-law→SIA boundary**: per level `pressure = EC.pressure(depth)`,
  `stress = |∇h|·p`, `flow_n(stress, E, p, grain_size, ...)`, then
  `D = ∫(H−z)·δ dz` (trapezoid). Output: D field + `m_D_max`.
- `compute_diffusive_flux(h_x, h_y, D, result)` — q = −D·∇h.
- `compute_3d_horizontal_velocity(geometry, h_x, h_y, sliding_velocity, u, v)` —
  U(z) = −I(z)·∇h + U_b (hybridization: sliding velocity added).

No elliptic solve in SIA — purely local explicit computation on a 5-point
staggered stencil. `BedSmoother` (preprocessing): `preprocess_bed(topg)`,
`theta(usurf, result)` (Schoof roughness), `smoothed_thk(...)`.

### 2.3 SSA — `stressbalance/ssa/`

- `SSA` (abstract): `update(inputs, full_update)`, pure virtual `solve(inputs)`.
  Inputs: thickness, surface, bed, sea level, tauc, enthalpy, fracture density,
  BC mask/values. Outputs: `m_velocity` (2D), `m_basal_frictional_heating`. **D**
  (steady nonlinear elliptic solve; iterative but instantaneous).
- `SSAFDBase` (FD discretization): `initialize_iterations(inputs)` (driving
  stress τ_d = −ρgH∇h, average hardness via `rheology::averaged_hardness`),
  `compute_nuH(...)` (νH = effective_viscosity·H·e^(−1/n) + ε),
  `fd_operator(...)` (13-point stencil; basal drag β = `basal_sliding_law->drag(tauc, v)`).
- `SSAFD` — Picard iteration over νH (`solve`, `picard_iteration`,
  `picard_manager`); inner KSP solve.
- `SSAFD_SNES` — Newton via PETSc SNES (`solve`, `compute_jacobian`).
- `SSAFEM` — FEM (Q1) Newton (`solve`, `PointwiseNuHAndBeta`, `compute_local_function/jacobian`).

### 2.4 Sliding laws

- `WeertmanSliding` (`stressbalance/WeertmanSliding.{hh,cc}`):
  `update(inputs, full_update)` — u_s = −(2A_s/(1−k))·(N·|∇h|)^(n−1)·∇h. **D**.
- `ZeroSliding`, `PrescribedSliding` — trivial/file-driven.

### 2.5 Timestepping — `stressbalance/timestepping.{hh,cc}`

Free functions: `max_timestep_cfl_3d(thickness, mask, u3, v3, w3)`,
`max_timestep_cfl_2d(thickness, mask, velocity)`,
`max_timestep_diffusivity(D_max, dx, dy, ratio)`. **D**.

## 3. Energy — `energy/` → `ice_energy.esm`

- `EnergyModel` (base): `update(t, dt, inputs)` → `update_impl`. Inputs
  (`energy::Inputs`): cell_type, basal_frictional_heating, basal_heat_flux,
  thickness, surface_liquid_fraction, shelf_base_temp, surface_temp,
  till_water_thickness, volumetric_heating_rate, u3, v3, w3. Outputs:
  `enthalpy()` (state), `basal_melt_rate()`. **I**.
- `EnthalpyModel` (default): `update_impl` — per-column `enthSystemCtx` solve
  (BOMBPROOF λ scheme; horizontal advection explicit upwind, vertical implicit).
  Basal BC decision (cold–temperate transition at the bed): floating →
  Dirichlet from `shelf_base_temp`; grounded warm+wet → Neumann/Dirichlet;
  cold → Neumann `q·n = q_lith + F_b`. **I**.
- `TemperatureModel` (cold, legacy): `update_impl` — per-column
  `tempSystemCtx` solve; post-processing clamps to Tpmp, drainage, bulge
  limiter; hard low-temperature check (`T < 200 K`, count > 10 → abort).
- `enthSystemCtx` / `tempSystemCtx` (`energy/enthSystem.{hh,cc}`,
  `energy/tempSystem.{hh,cc}`): per-column tridiagonal systems.
  `init(i, j, ismarginal, H)`, `compute_enthalpy_CTS()` (E_s(p) per level),
  `compute_lambda()`, `assemble_R()` (diffusivity switch at CTS), `solve(x)`. **I**.
- `BedThermalUnit` (base): `update(bedrock_top_temperature, t, dt)` →
  `flux_through_top_surface()` (G₀, the `basal_heat_flux` input to the energy
  model). **I**. `BTU_Full` (per-column `BedrockColumn::solve`, backward Euler),
  `BTU_Minimal` (identity: G₀ ≡ bheatflx).
- `DrainageCalculator`: `get_drainage_rate(omega)` — piecewise-linear D(ω). **D**.
- `CHSystem` (cryo-hydrologic warming, regional): `update_impl`; coupling
  function `cryo_hydrologic_warming_flux(k, R, H, E_ice, E_ch, result)`. **D**.

## 4. Age — `age/` → `ice_age.esm`

- `AgeModel`: `update(t, dt, inputs{thickness, u3, v3, w3})` → `m_ice_age`
  (state). PDE: ∂τ/∂t + u·∇τ = 1; BC: age 0 at surface, 0 at base when w>0. **I**.
- `AgeColumnSystem`: `init(i, j, thickness)`, `solve(x)` — tridiagonal,
  explicit upwind horizontal + implicit vertical advection. **I**.

## 5. Basal strength — `basalstrength/` → `basal_resistance.esm`

- `YieldStress` (base): `update(inputs{geometry, till_water_thickness,
  subglacial_water_thickness, no_model_mask}, t, dt)` →
  `basal_material_yield_stress()` (tauc). **D** for the tauc computation.
- `MohrCoulombYieldStress` (default): `update_impl` — τ_c = c₀ + tan(φ)·N_till
  (Tulaczyk et al. 2000 effective pressure via `MohrCoulombPointwise::effective_pressure`).
  State: `m_till_phi` (till friction angle). **D**.
- `MohrCoulombPointwise`: `yield_stress(delta, P_overburden, water_thickness, phi)`,
  `till_friction_angle(...)`, `effective_pressure(...)` — pure pointwise
  functions, ideal `.esm` test targets. **D**.
- `ConstantYieldStress` — tauc ≡ config value. **D** (identity).
- `OptTillphiYieldStress` — inverse optimization of φ at intervals. **I** (φ), **D** (tauc).
- Basal resistance laws (`basalstrength/basal_resistance.{hh,cc}`):
  `IceBasalResistancePlasticLaw` / `IceBasalResistancePseudoPlasticLaw` /
  `IceBasalResistanceRegularizedLaw` — `drag(tauc, vx, vy)` → β (basal shear
  stress τ_b = −β·V) and `drag_with_derivative(tauc, vx, vy, *drag, *ddrag)`. **D**.
  Consumed by SSA at `SSAFDBase::fd_operator` (β added to matrix diagonal).

## 6. Couplers — `coupler/` → `surface_forcing.esm`, `ocean_coupling.esm`

### 6.1 Atmosphere — `coupler/atmosphere/`

Base `AtmosphereModel`: `update(geometry, t, dt)` → `precipitation()`,
`air_temperature()`. **D**. Models: `Given` (`"given"`), `PIK` (`"pik"`),
`SeaRISEGreenland` (`"searise_greenland"`), `CosineYearlyCycle` (`"yearly_cycle"`),
`Uniform`, `WeatherStation`; modifiers `Anomaly`, `Delta_T`, `Delta_P`, `Frac_P`,
`ElevationChange`, `OrographicPrecipitation`, `PrecipitationScaling`.

### 6.2 Ocean — `coupler/ocean/`

Base `OceanModel`: `update(inputs{geometry}, t, dt)` → `shelf_base_temperature()`,
`shelf_base_mass_flux()`, `average_water_column_pressure()`. **D**. Models:
`Given` (`"given"`), `GivenTH` (`"th"`, Holland–Jenkins 3-equation),
`PIK` (`"pik"`, Beckmann–Goosse quadratic), `Constant` (`"constant"`),
`Pico` (`"pico"`, 3-D box model), `Cache`; modifiers `Anomaly`, `Delta_SMB`,
`Frac_SMB`, `Delta_T`, `Runoff_SMB`, `Delta_MBP`, `Frac_MBP`.
Sea level: `coupler/SeaLevel.hh` — `update(geometry, t, dt)` → `elevation()`.

### 6.3 Surface — `coupler/surface/`

Base `SurfaceModel`: `update(geometry, t, dt)` → `accumulation()`, `melt()`,
`runoff()`, `mass_flux()`, `temperature()`, layer state. Models: `Given` (`"given"`),
`Simple` (`"simple"`), `PIK` (`"pik"`), `TemperatureIndex` (`"pdd"`, PDD scheme
with firn/snow state — **I** for totals/state), `ISMIP6` (`"ismip6"`),
`DEBMSimple` (`"debm_simple"`), `Elevation`; modifiers `Anomaly`, `Cache`,
`Delta_T`, `ForceThickness` (**I**), `ElevationChange`, `NoGLRetreat`.

### 6.4 Frontal melt — `coupler/frontalmelt/`

Base `FrontalMelt`: `update(inputs{geometry, subglacial_water_flux}, t, dt)` →
`frontal_melt_rate()`, `retreat_rate()`. **D**. Models: `Given`, `Constant`,
`DischargeGiven`, `DischargeRouting`; physics kernel `FrontalMeltPhysics`
(`frontal_melt_from_undercutting`, `frontal_melt_from_ismip6`).

## 7. Hydrology — `hydrology/` → `hydrology.esm`

Base `Hydrology`: `update(t, dt, inputs{no_model_mask, geometry,
surface_input_rate, basal_melt_rate, ice_sliding_speed})` →
`till_water_thickness()` (W_till), `subglacial_water_thickness()` (W),
`overburden_pressure()`, `flux()`. **I** (internal substeps; one-way coupling).
Models: `Routing` (`"routing"`, Shreve: q = −K∇ψ), `Distributed` (`"distributed"`,
linked-cavity with P as state), `NullTransport` (`"null"`, till-can),
`SteadyState` (`"steady"`, + `EmptyingProblem`). Effective pressure is computed
downstream in `MohrCoulombYieldStress`, not here.

## 8. Bed deformation — `earth/` → `bed_deformation.esm`

Base `BedDef`: `update(ice_thickness, sea_level_elevation, t, dt)` →
`bed_elevation()`, `uplift()`. Load accumulator integrates ice-equivalent
thickness over `dt`; `update_impl(load, t, dt)` called at
`bed_deformation.update_interval`. **I**. Models: `Null`, `PointwiseIsostasy`
(local relaxing half-space), `Given` (file), `LingleClark` (wrapper) /
`LingleClarkSerial` (Fourier spectral collocation, viscous half-space + elastic
lithosphere; `step(dt, H)` → total/viscous/elastic displacement).

## 9. Front retreat / calving — `frontretreat/` → `calving.esm`

- `FrontRetreat::update_geometry(dt, geometry, bc_mask, retreat_rate, href, thickness)`
  — applies retreat rate over `dt` (part-grid). **I**.
- `PrescribedRetreat::update(t, dt, thickness, area_specific_volume)` — ISMIP6
  parameterized retreat mask. **I**.
- Calving (`frontretreat/calving/`): `EigenCalving::update(cell_type, velocity)`
  (m_K·max eigen strain rate), `vonMisesCalving::update(cell_type, thickness,
  velocity, enthalpy)`, `HayhurstCalving::update(cell_type, thickness, sea_level, bed)`,
  `CalvingAtThickness::update(t, dt, cell_type, thickness)` (in-place),
  `FloatKill::update(cell_type, thickness)` (in-place). Rates **D**; in-place
  removal **I**. Orchestration: `icemodel/frontretreat.cc` sums rates into
  `retreat_rate` then calls `FrontRetreat::update_geometry`.

## 10. Geometry — `geometry/`

- `GeometryEvolution`: `flow_step(geometry, dt, advective_velocity,
  diffusive_flux, bc_mask)` and `source_term_step(geometry, dt, bc_mask,
  smb_rate, basal_melt_rate)` — mass continuity over `dt`. **I**.
  Contract: `H_change + Href_change = dt·(SMB + BMB − flux_divergence) + err`.
- `part_grid_threshold_thickness(cell_type, thickness, surface, bed)` —
  pointwise geometric function. **D**.

## 11. Instrumentation notes

- Existing debug hook: `StressBalance::Inputs::dump(filename)` (NetCDF).
- Diagnostic scaffolding: `SSB_diagnostics.hh` (beta, taub, taud),
  `SIAFD_diagnostics.hh`, `StressBalance_diagnostics.hh`,
  `SSAFDBase::integrated_viscosity()` / `driving_stress()`.
- Preferred **D** boundaries for first instrumentation: `MohrCoulombPointwise::*`,
  `IceBasalResistancePlasticLaw::drag*`, `DrainageCalculator::get_drainage_rate`,
  `StressBalance::compute_volumetric_strain_heating`,
  `StressBalance::compute_vertical_velocity`,
  `ShallowStressBalance::compute_basal_frictional_heating`,
  `cryo_hydrologic_warming_flux`, `MohrCoulombYieldStress::update_impl` (tauc map),
  `SIAFD::compute_diffusivity` (flow law → D), `FlowLaw::*` (softness/hardness/flow).
- **I** boundaries needing (t, dt) + state in/out dumps: `EnthalpyModel::update_impl`,
  `enthSystemCtx::solve`, `BTU_Full::update_impl`, `AgeModel::update`,
  `AgeColumnSystem::solve`, `Hydrology::update`, `BedDef::update`,
  `FrontRetreat::update_geometry`, `GeometryEvolution::flow_step`/`source_term_step`.
