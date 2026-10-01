# Stage 1 — Reference Run Configurations

Authoritative run configurations for Stage-1 instrumentation (PLAN.md §4,
step 4: *record the exact run configuration for reproducibility*). All runs use
the PISM v2.3.2 build in this repo (`stage1/build/install`, see
`stage1/build/build-pism.sh`) and the shared environment in `stage1/runs/env.sh`.

```bash
source stage1/runs/env.sh   # sets PATH, LD_LIBRARY_PATH, PETSC_DIR, ...
```

## Verification tests (self-contained, no input files)

PISM's built-in verification tests exercise the component inventory with exact
solutions — ideal for instrumentation because they are reproducible and cover
each physics area. Canonical invocations (from the PISM manual,
`doc/sphinx/manual/verification/index.rst`):

| Test | Components exercised | Canonical invocation | Status |
|---|---|---|---|
| A | isothermal SIA, rheology (isothermal_glen) | `pism -test A -Mx 61 -My 61 -Mz 11 -y 25000` | ✅ runs |
| B | isothermal SIA, similarity solution | `pism -test B -Mx 61 -My 61 -Mz 11 -ys 422.45 -y 25000` | ✅ runs |
| C | isothermal SIA, growing accumulation | `pism -test C -Mx 61 -My 61 -Mz 11 -y 15208.0` | ✅ runs |
| D | isothermal SIA, oscillating accumulation | `pism -test D -Mx 61 -My 61 -Mz 11 -y 25000` | ✅ runs |
| F | thermomechanical SIA (energy, strain heating) | `pism -test F -Mx 61 -My 61 -Mz 61 -y 25000` | ✅ runs |
| G | thermomechanical SIA, oscillating accumulation | `pism -test G -Mx 61 -My 61 -Mz 61 -y 25000` | ✅ runs |
| H | bed deformation (iso) + isothermal SIA | `pism -test H -Mx 61 -My 61 -Mz 11 -y 40034 -bed_def iso` | ✅ runs |
| K | pure conduction in ice and bedrock (BTU) | `pism -test K -Mx 6 -My 6 -Mz 401 -Mbz 101 -y 130000` | ✅ runs |
| L | isothermal SIA, non-flat bed | `pism -test L -Mx 61 -My 61 -Mz 31 -y 25000` | ✅ runs |
| V | SSA flow line (van der Veen) | `pism -test V -y 1000 -part_grid -ssa_method fd -cfbc -Mx 51` | ❌ fails — see bugs.md #3 |

Smoke-tested (short runs, small grids) on 2026-09-30 with the repo build:
A, F, G, H, K all complete and report error norms vs. the exact solutions.

## SSA verification executables

PISM's own CI uses the dedicated SSA test executables (no energy model), which
is the right tool for SSA boundary instrumentation:

```bash
$PISM_BIN/pism_ssa_test_cfbc -Mx 201 -verbose 1 -o_size none -My 3 \
  -ssafd_ksp_type richardson -ssafd_pc_type lu -stress_balance.ssa.epsilon 0
$PISM_BIN/pism_ssa_test_i -ssa_method fd -Mx 5 -My 500 -ssafd_picard_rtol 1e-6 -ssafd_ksp_rtol 1e-11
$PISM_BIN/pism_ssa_test_j -ssa_method fd -Mx 60 -My 60 -ssafd_ksp_rtol 1e-12
$PISM_BIN/pism_ssa_test_const -Mx 61 -My 61   # constant viscosity SSA
$PISM_BIN/pism_ssa_test_linear -Mx 61 -My 61  # linear (Stokes) SSA
```

## BTU verification

```bash
$PISM_BIN/pism_btutest    # bedrock thermal unit verification (test K physics)
```

## Notes

- Grid sizes above are the canonical ones; smaller grids (`-Mx 31` etc.) are
  fine for smoke tests but error norms are only meaningful at canonical sizes.
- `-test V` (full-model SSA) currently fails in v2.3.2 — see bugs.md #3. The
  SSA boundary is instrumented via `pism_ssa_test_*` instead.
- Run outputs (NetCDF) are temporary: write them to `/scratch.local/ctessum/`
  (may disappear between sessions), not into the repo.
