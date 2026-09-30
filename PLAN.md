# EqIce — PISM Implementation Plan (EarthSciAST `.esm`)

This document is the working plan for re-implementing the PISM ice-sheet model
(https://github.com/pism/pism/) as a composition of EarthSciAST `.esm`
components, following the staged process described in `AGENTS.md`. It is a
living document: update it as components are instrumented, migrated, and merged.

## 1. Context and objectives

- Goal: a top-level `eqice.esm` model, assembled by reference from individual
  `.esm` components, that reproduces PISM simulations.
- The original PISM implementation (hereafter "the fortran codes") is the
  reference for physics and numerics. Components that already exist in
  [EarthSciModels](https://github.com/EarthSciML/EarthSciModels) must be reused
  (and their tests extended from instrumented fortran outputs), not duplicated.
- All model logic lives in `.esm` files. Tests and examples live in the
  `tests` and `analysis` sections of `.esm` files.
- All local tests run through the EarthSciAST Rust CLI binary, kept untracked
  in this repo's root. Final CI gates run each file's inline `tests` through the
  EarthSciModels `tools/run_esm_inline_tests.py` / `run_esm_tests.jl` runners
  (the single-pathway rule).

## 2. Tooling and conventions

- **Format:** `.esm`, JSON, current version `1.0.0` (esm-spec §2).
  Variables are `unknown` or `parameter`; everything else is derived (§6.3).
- **Authoring stance (esm-spec §1.1):** AST first, closed-function registry
  second, factoring third. Use `expression_templates` liberally; import shared
  templates by reference (esm-spec §9.6/§9.7) instead of repeating logic.
- **PDEs:** write them as PDEs, not einsum/aggregate ops. Reuse (or submit PRs
  to) [EarthSciDiscretizations](https://github.com/EarthSciML/EarthSciDiscretizations)
  grid+stencil+rule templates; import them with metaparam bindings (e.g. `N`).
  Boundary conditions live inside the rule.
- **Reuse:** check EarthSciModels (`components/<domain>/`) and
  EarthSciDiscretizations first. Some PISM physics (e.g. vertical diffusion,
  enthalpy/heat transport) may already exist as components or discretization
  rules; add tests rather than re-deriving.
- **No scripts for equation authoring.** Scripts may only (a) extract numeric
  test tuples from fortran dumps into a hand-authored `tests` block, and (b)
  translate tabular data (lookup tables, `.eqn` reaction lists) into
  `function_tables`/`reaction_systems`/`data_sources`.

## 3. PISM component inventory

Map PISM process components to `.esm` components/subassemblies. Group by the
PISM source tree so each maps to a well-defined input/output boundary for
instrumentation.

| PISM area | PISM component / process | `.esm` component (draft) | Existing impl? |
|---|---|---|---|
| Rheology | Glen–Nye flow law, Arrhenius/temperature-dependent viscosity, enhancement factors | `ice_rheology.esm` (expression templates) | check |
| Stress balance | SIA (shallow-ice approx.), SSA (shallow-shelf approx.), SIA+SSA coupling | `sia.esm`, `ssa.esm`, `stress_balance.esm` | check (SIA/SSA discretizations) |
| Energy | Enthalpy/temperature evolution, cold-temperate transition, vertical advection+diffusion | `ice_energy.esm` | check vertical diffusion |
| Age | Ice-age transport | `ice_age.esm` | check |
| Bed deformation | Elastic lithosphere + relaxing asthenosphere (Lingle–Clark / EarthDeform) | `bed_deformation.esm` | check (geodynamics) |
| Basal resistance | Sliding laws (Weertman, Budd, Coulomb) | `basal_resistance.esm` | check |
| Subglacial hydrology | Water thickness/effective pressure | `hydrology.esm` | check |
| Ocean coupling | Sub-shelf melt (PICO/PIK/thin–film), thermal forcing | `ocean_coupling.esm` | check |
| Atmosphere / surface | Surface mass balance, temperature forcing, precipitation | `surface_forcing.esm` | check |
| Calving / front retreat | Eigen-calving, von Mises, thickness threshold | `calving.esm` | check |

Each row becomes one or more component PRs into EarthSciModels (stage 2). The
top-level `eqice.esm` (stage 3) composes the merged components.

## 4. Stage 1 — Instrument the fortran codes

Purpose: produce authoritative input/output traces to build tests from.

1. Identify discrete subassembly boundaries in the fortran codes matching the
   component inventory above.
2. Instrument each boundary to dump inputs and outputs during simulation.
   Prefer **instantaneous derivatives**; use integrated trajectories only where
   derivatives are unavailable.
3. For each component, also dump inputs/outputs of significant internal
   subassemblies (e.g. SSA linear solve vs. SIA stencil), to support later
   factoring.
4. Record the exact run configuration (grid, time step, domain, forcing) for
   reproducibility, and store a curated set of dumps per component.
5. Extract a select subset of these dumps into numeric test tuples.

Deliverable: a catalog of instrumented boundaries, their dump files, and the
run configurations. Keep these dumps on disk (not in `/tmp` — RAM-backed) and
reference them from the plan.

## 5. Stage 2 — Stubs → physics → EarthSciModels PRs

For each component (bottom-up; leaf physics first, then subassemblies):

1. **Hand-author a stub** `.esm` declaring the component's unknowns,
   parameters, and a `tests` block with the Stage-1 test tuples. The stub
   equations are placeholders (or absent) so tests initially fail/define the
   contract.
2. **Fill in the physics** in the stub so the tests pass:
   - Author equations as math (PDEs where the fortran discretizes PDEs),
     importing discretization rules from EarthSciDiscretizations.
   - Factor repeated math into `expression_templates` and import them by
     reference.
3. **Review + merge** each PR into EarthSciModels (`components/<domain>/`),
   after a human review, following the existing layout and `reference`/`notes`
   documentation style of sibling components.
4. **Cascade:** once leaf components are merged, create subassembly `.esm`
   files that import the leaves **by reference** and couple them, adding
   subassembly tests. Repeat until the component's subassemblies pass.
5. For any PISM physics already present in EarthSciModels, add the Stage-1
   test tuples to the existing `.esm` rather than writing a new one.

## 6. Stage 3 — Top-level `eqice.esm`

1. Import the merged components by reference into `eqice.esm`.
2. Add the couplings (stress balance ↔ energy, basal resistance ↔ hydrology,
   surface/ocean forcing, bed deformation, calving/front retreat).
3. Test against full-fortran-model simulations (integrated trajectories over
   realistic domains), using the same run configurations as Stage 1.
4. Iterate on coupling order/timestep and any subassembly discrepancies.

## 7. Repository layout (this repo)

```
eqice/
  AGENTS.md            # process rules (authoritative)
  PLAN.md              # this document
  <rust CLI binary>    # untracked local EarthSciAST runner
  stage1/
    dumps/             # instrumented fortran I/O traces (referenced by tests)
    runs/              # run configurations
  stage2/
    <component>.esm    # working stubs/physics before merging upstream
  stage3/
    eqice.esm          # top-level coupled model
  bugs.md              # centralized list of bugs/errors found in fortran
```

Components that pass review are removed from `stage2/` once merged into
EarthSciModels and imported by reference from there.

## 8. Centralized fortran bug list

Maintain `bugs.md` (and a matching entry wherever AGENTS.md points). Record for
each bug: location (file:line/function), description, impact on results,
whether the `.esm` reproduces or corrects it, and status.

## 9. Test strategy summary

- Leaf physics: instantaneous-derivative tuples (preferred).
- Subassemblies: composed-component tests with the same tuples.
- Full model: integrated trajectories vs. full fortran simulations.
- All tests are inline `tests` blocks runnable by the Rust CLI locally and by
  the EarthSciModels CI gate on merge.

## 10. Milestones / sequencing

1. Confirm PISM source access + build; pick reference run configs.
2. Stage-1 instrumentation scaffolding + first component traces (rheology, SIA).
3. Stage-2: rheology → stress balance (SIA, SSA) → energy → basal/hydrology →
   bed → surface/ocean/calving → subassemblies.
4. Stage-3: `eqice.esm` coupling + full-model validation.

Each milestone closes with merged PRs and passing inline tests.
