# Bugs and issues found in the PISM C++ codes

Centralized list per PLAN.md §8. Each entry records: location, description,
impact, whether the `.esm` reproduces or corrects it, and status.

## 1. PISM 2.3.2 configure hard-fails when UDUNITS-2 is absent

- **Location:** `CMake/FindUDUNITS2.cmake` (lines 63–86), triggered from
  `CMake/PISM_CMake_macros.cmake` line 157.
- **Description:** `pism_find_prerequisites` calls `find_package(UDUNITS2)`
  without `REQUIRED`, but the Find module itself runs a compile+run probe and,
  when the probe fails (no udunits2), falls through to
  `find_package(EXPAT REQUIRED)`. On a system with neither udunits2 nor
  expat development files, configure aborts with
  `Could NOT find EXPAT (missing: EXPAT_LIBRARY EXPAT_INCLUDE_DIR)`.
  UDUNITS-2 is effectively a **required** dependency to build PISM: the
  bundled `src/external/calcalcs/utCalendar2_cal.c` unconditionally includes
  `udunits2.h` and links `libudunits2` (calendar conversions used by the time
  module). The CMake code, however, treats it as optional, so the failure mode
  is confusing.
- **Impact on results:** none (build-time issue only).
- **Workaround / patch applied locally** (`stage1/pism/CMake/FindUDUNITS2.cmake`):
  `find_package(EXPAT QUIET)`; if EXPAT is unavailable, report
  "building PISM without UDUNITS-2 support" and set `UDUNITS2_FOUND FALSE`
  with empty include/library variables. (With no udunits2 at all the build
  still fails later in `calcalcs`; we built udunits2 from source —
  see `stage1/build/build-local-deps.sh`.)
- **Reproduced/corrected in .esm:** n/a (build tooling).
- **Status:** patched locally; candidate upstream bug report.

## 2. PISM_CONFIG_FILE default points at /usr/local when installed to a custom prefix

- **Location:** `CMakeLists.txt` lines 64–66
  (`set (Pism_CONFIG_FILE "${CMAKE_INSTALL_FULL_DATADIR}/pism/pism_config.nc")`).
- **Description:** when building with a non-default `CMAKE_INSTALL_PREFIX`,
  the compiled-in default config path is still derived from the install prefix
  *at configure time*; there is no environment-variable override (e.g.
  `PISM_PREFIX` only affects shell scripts in examples). Running the binary
  from the build tree without installing fails with
  `PISM ERROR: No such file or directory ... pism_config.nc`. This is standard
  CMake behavior (not strictly a bug), but it is a footgun for stage-1
  instrumentation runs from the build tree.
- **Impact on results:** none; resolved by `make install` to
  `stage1/build/install` (see `stage1/build/build-pism.sh`).
- **Status:** documented; no code change needed.

## 3. PISM 2.3.2 full-model verification test V aborts with "too many low temps"

- **Location:** `src/icemodel/IceModel.cc` (step ordering: `energy_step` at
  line 491 runs before `m_ocean->update` at line 569);
  `src/coupler/ocean/Constant.cc` (`init_impl` does not compute
  `shelf_base_temperature`); `src/verification/iceCompModel.cc`
  (`allocate_energy_model` unconditionally creates `TemperatureModel_Verification`).
- **Description:** `pism -test V ...` (the full-model SSA verification test)
  aborts at the first time step with `PISM ERROR: too many low temps: 18`.
  Root cause chain:
  1. Verification mode always uses `ocean::Constant` and the temperature-based
     energy model (`IceCompModel::allocate_energy_model` ignores `energy.model`).
  2. `ocean::Constant::init_impl` computes only the water-column pressure;
     `shelf_base_temperature` is computed only in `update_impl`.
  3. In `IceModel::step`, `energy_step` runs **before** `m_ocean->update`, so at
     the first step the energy model sees `shelf_base_temp = 0` K.
  4. Floating cells (the upstream 600 m columns of test V) get a 0 K basal
     Dirichlet BC → the column solve produces ~0–151 K near the base → the
     temperature model's hard low-temperature check (`T < 200 K`, count > 10)
     aborts the run.
  Instrumented debug output at the first step (column i=0, j=0, H=600 m,
  mask=3): `T_shelf=0.000`, solution `x = [0.00, 151.26, 221.51, ...]`.
- **Impact on results:** test V cannot be run as documented in the PISM manual
  for v2.3.2. The dedicated SSA executables (`pism_ssa_test_cfbc`, used by
  PISM's own CI) are unaffected — they do not use the energy model. The same
  0 K shelf-base glitch affects the *first* energy step of any run with
  floating ice and the default `constant` ocean model, but the enthalpy model
  has no hard low-temperature abort, so normal runs recover after one step.
- **Reproduced/corrected in .esm:** n/a (C++ runtime issue). The `.esm`
  ocean-coupling component should compute `shelf_base_temperature` at
  initialization, not rely on a first update.
- **Status:** documented; workaround for Stage 1 is to instrument the SSA
  boundary via `pism_ssa_test_*` (see `stage1/runs/README.md`). Candidate
  upstream fix: compute `shelf_base_temperature` in `Constant::init_impl`, or
  move `m_ocean->update` before `energy_step` in `IceModel::step`.

## 5. ConstantPIK.cc sign comment contradicts the code: shelf_base_mass_flux is positive for melting

- **Location:** `src/coupler/ocean/ConstantPIK.cc` lines 136–137
  (`PIK::mass_flux`, the Beckmann–Goosse sub-shelf parameterization).
- **Description:** the comment claims "shelfbmassflux is positive if ice is
  freezing on; here it is always negative" and "same sign as ocean_heat_flux
  (positive if massflux FROM ice TO ocean)". Both claims are wrong for the
  code that follows: `result = ocean_heat_flux / L`, with
  `ocean_heat_flux = melt_factor * rho_w * c_p * gamma_T * (T_ocean - T_f)`
  and `T_f = 273.15 + 0.0939 - 0.057*S + 7.64e-4*z_b` (`z_b = -(rho_i/rho_w) H
  <= 0`). With the default `T_ocean = -1.7 C` the flux is **positive
  everywhere** (`T_ocean - T_f > 0`), and the downstream consumers
  (`IceModel::combine_basal_melt_rate` → `GeometryEvolution` with
  `dH_BMB = -dt*mf/rho_i`) interpret a positive shelf-base mass flux as
  **melting** (thickness loss). So the sign convention is "positive =
  melting", the flux is never negative under the shelf, and the comment is
  wrong on both points.
- **Impact on results:** none (comment only) — but it is a trap for a
  reimplementation: an `.esm` ocean-coupling component must reproduce
  `mf = Q/L` positive = melting, not the comment's "always negative /
  positive = freezing". The Stage-1 surface/ocean dumps record the actual
  sign.
- **Reproduced/corrected in .esm:** the surface/ocean Stage-1 traces record
  `shelf_base_mass_flux > 0` (melting) under the shelf; the Stage-2
  ocean-coupling component will reproduce the code (positive = melting).
- **Status:** documented; candidate upstream comment fix.

## 4. StressBalance::Inputs::dump() aborts: writes an undefined config variable

- **Location:** `src/stressbalance/StressBalance.cc`
  (`Inputs::dump(const char *filename)`), calls
  `io::write_config(*config, "pism_config", output)`.
- **Description:** `Inputs::dump()` is dead code in v2.3.2 (never called by the
  model). It creates an output file with `SynchronousOutputWriter` but never
  `define_variable`s the `pism_config` variable before `write_config()` calls
  `file.write_text("pism_config", ...)`. Running it aborts with
  `PISM ERROR: NetCDF: Variable not found while writing variable 'pism_config'`.
  This is the natural input-dump hook for Stage-1 instrumentation of the
  stress-balance boundary.
- **Impact on results:** none for model runs (dead code); a Stage-1
  workaround was needed — the SIA driver (`stage1/instrument/instrument_sia.cc`)
  writes the input+output NetCDF directly with the same writer/define/write
  pattern instead of calling `Inputs::dump()`.
- **Reproduced/corrected in .esm:** n/a (C++ I/O issue).
- **Status:** documented; candidate upstream fix: `file.define_variable(...)`
  for `pism_config` before writing, or drop the `write_config` call from
  `Inputs::dump()`.

---
(no `.esm`-reproducible physics bugs found yet — Stage 2 will populate this file
as component tests are compared against the C++ implementations.)
