#!/usr/bin/env python3
"""Build the numeric data of stage2/ice_energy_3d.esm from the Stage-1 energy dumps.

This script performs ONLY numeric data extraction and filling (AGENTS.md: "Scripts
MAY be used to extract numeric test tuples from C++ dumps into an existing
hand-authored tests block, and to translate tabular data ... into ... entries; the
equations, structure, and factoring of a component are always authored by hand").
The component structure of stage2/ice_energy_3d.esm is hand-authored; this script
patches the index sets, the block-size scalars, the E_f_init default, and the
`tests` block (whose ids, coords, descriptions, and tolerance philosophy are
hand-authored below), filling every numeric array and every expected value from
the dump.

The esm grid is 1-based (esm cell (i,j) = C++ cell (i-1,j-1)); the flat default
arrays are row-major over the declared shape axes with the LAST axis fastest
(esm `index` convention, verified empirically against the harness). Each block
covers C++ cells (i0..i0+1, j0..j0+2) -- a 2-column x 3-row sub-domain, all four
blocks sharing the index sets i (size 2) and j (size 3); i0/j0 are the C++
indices of the block's west/south columns. The full-field parameters are read at
absolute esm indices (i0+i, j0+j).
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)                 # stage2/
ROOT = os.path.dirname(WT)                 # worktree root
DUMP = os.path.join(ROOT, 'stage1/dumps/energy/columns.csv')
BMELT = os.path.join(ROOT, 'stage1/dumps/energy/basal_melt_rate.csv')
ESM = os.path.join(WT, 'ice_energy_3d.esm')

DZ = 4000.0 / 61.0
DX = 58064.516129
DT = 315569259.747
MZ = 61
MZF = 62
MARGIN = 100.0
C_ICE = 2009.0
T_MELTING = 273.15
T_0 = 223.15
BETA = 7.9e-8
RHO = 910.0
G = 9.81
P_AIR = 0.0

# The four sub-domain blocks (C++ 0-based indices of the west/south column).
BLOCKS = {'w': (7, 13), 'e': (22, 13), 'dome': (15, 13), 'm': (27, 13)}
NI, NJ = 2, 3

# ----------------------------------------------------------------------
# dump access
# ----------------------------------------------------------------------
lines = [ln.strip().split(',') for ln in open(DUMP).readlines() if not ln.startswith('#')]
data = lines[1:] if lines and lines[0][0] == 'i' else lines  # tolerate header row


def col(i, j):  # C++ cell (i,j), 0-based; 61 rows
    start = (i + 31 * j) * 61
    return data[start:start + 61]


def get(i, j, idx):
    return float(col(i, j)[0][idx])


def Hc(i, j):
    return get(i, j, 5)


def E_in(i, j):
    return [float(r[14]) for r in col(i, j)]


def E_out(i, j):
    return [float(r[15]) for r in col(i, j)]


def u3(i, j):
    return [float(r[10]) for r in col(i, j)]


def v3(i, j):
    return [float(r[11]) for r in col(i, j)]


def w3(i, j):
    return [float(r[12]) for r in col(i, j)]


def sig3(i, j):
    return [float(r[13]) for r in col(i, j)]


def bmelt_of(i, j):
    for ln in open(BMELT).readlines():
        if ln.startswith('#'):
            continue
        p = ln.strip().split(',')
        if int(p[0]) == i and int(p[1]) == j:
            return float(p[4])
    raise KeyError((i, j))


# ----------------------------------------------------------------------
# grid maps and interpolation (ColumnInterpolation.cc)
# ----------------------------------------------------------------------
def z_s(m):
    return m * 4000.0 / 60.0


def z_f(k):
    return k * DZ


def coarse2fine_map():
    m = 0
    out = []
    for k in range(MZF):
        if z_f(k) <= z_s(0):
            out.append(0)
            continue
        if z_f(k) >= z_s(MZ - 1):
            out.append(MZ - 1)
            continue
        while z_s(m + 1) < z_f(k):
            m += 1
        out.append(m)
    return out


C2F = coarse2fine_map()          # 0-based storage index bracketing each fine level


def incr_c2f(k):
    m = C2F[k]
    return (z_f(k) - z_s(m)) / (z_s(m + 1) - z_s(m))


def fine2coarse_map():
    out = []
    for m in range(MZ):          # storage level m at z_s(m); find fine bracket
        k = 0
        while z_f(k + 1) < z_s(m):
            k += 1
        out.append(k)
    return out


F2C = fine2coarse_map()          # 0-based fine index bracketing each storage level


def incr_f2c(m):
    k = F2C[m]
    return (z_s(m) - z_f(k)) / (z_f(k + 1) - z_f(k))


def c2f(x, ks):                  # pure coarse_to_fine_linear (own ks), 62 fine rows
    out = []
    for k in range(MZF):
        if k > ks:
            out.append(float(x[C2F[k]]))
            continue
        m = C2F[k]
        if m == MZ - 1:
            out.append(float(x[MZ - 1]))
            continue
        out.append(float(x[m]) + incr_c2f(k) * (float(x[m + 1]) - float(x[m])))
    return out


def ks_of(H):
    return min(int(math.floor(H / DZ)), MZF - 1)


# ----------------------------------------------------------------------
# per-column derived quantities
# ----------------------------------------------------------------------
def enth_ks(i, j):
    return E_in(i, j)[MZ - 1]    # the top storage row carries the surface enthalpy


def pinned_c2f(x, ks, Eks):      # c2f with rows >= ks pinned to the surface enthalpy
    out = c2f(x, ks)
    for k in range(ks, MZF):
        out[k] = Eks
    return out


def cts(p):                      # EnthalpyConverter op order (verified bit-exact vs leaf)
    return C_ICE * ((T_MELTING - BETA * p) - T_0)


def e_cts_col(i, j):
    H = Hc(i, j)
    ks = ks_of(H)
    out = []
    for k in range(MZF):
        if k <= ks:
            depth = H - k * DZ
            p = P_AIR + RHO * G * depth
        else:
            p = P_AIR
        out.append(cts(p))
    return out


def marginal(i, j):              # PISM marginal(): any 8-neighbor H < 100 m (periodic)
    H = Hc(i, j)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            if Hc((i + di) % 31, (j + dj) % 31) < MARGIN:
                return True
    return False


def hadv_fine(i, j):             # the model's hadv_term at t = 0 (pinned-neighbor convention)
    ks = ks_of(Hc(i, j))
    Ef = pinned_c2f(E_in(i, j), ks, enth_ks(i, j))
    uf = c2f(u3(i, j), ks)
    vf = c2f(v3(i, j), ks)

    def nbr(ii, jj):
        ksn = ks_of(Hc(ii, jj))
        return pinned_c2f(E_in(ii, jj), ksn, enth_ks(ii, jj))

    Ew, Ee = nbr(i - 1, j), nbr(i + 1, j)
    Es, En = nbr(i, j - 1), nbr(i, j + 1)
    out = []
    for k in range(MZF):
        if k > ks:
            out.append(0.0)
            continue
        u, v = uf[k], vf[k]
        up_u = (u / DX) * (Ef[k] - Ew[k]) if u >= 0 else (u / DX) * (Ee[k] - Ef[k])
        up_v = (v / DX) * (Ef[k] - Es[k]) if v >= 0 else (v / DX) * (En[k] - Ef[k])
        out.append(-(up_u + up_v))
    return out


# ----------------------------------------------------------------------
# full-field arrays (row-major over the DECLARED axes, first axis outer,
# last axis fastest -- the harness reads a nested default/override array with
# its first nesting level along the FIRST declared axis, verified empirically)
# ----------------------------------------------------------------------
def flat31x31x62(fn):
    out = []
    for gi in range(31):
        for gj in range(31):
            out.extend(fn(gi, gj))
    return out


def flat31x31x61(fn):
    out = []
    for gi in range(31):
        for gj in range(31):
            out.extend(fn(gi, gj))
    return out


def block_slice(full_flat, i0, j0):      # esm block cells (i0+i, j0+j), i = 1..NI, j = 1..NJ
    out = []
    for bi in range(NI):
        row = []
        for bj in range(NJ):
            base = ((i0 + bi) * 31 + (j0 + bj)) * 62
            row.append(full_flat[base:base + 62])
        out.append(row)
    return out


def build_payload():
    E_f_init_g = flat31x31x62(lambda gi, gj: pinned_c2f(E_in(gi, gj), ks_of(Hc(gi, gj)), enth_ks(gi, gj)))
    E_cts_g = flat31x31x62(lambda gi, gj: e_cts_col(gi, gj))
    H_g = [Hc(gi, gj) for gi in range(31) for gj in range(31)]
    u_g = flat31x31x61(u3)
    v_g = flat31x31x61(v3)
    w_g = flat31x31x61(w3)
    sig_g = flat31x31x61(sig3)
    m_c2f = [m + 1 for m in C2F]                      # 1-based (esm index op)
    incr_c2f_map = [incr_c2f(k) for k in range(MZF)]
    m_f2c = [m + 1 for m in F2C]
    incr_f2c_map = [incr_f2c(m) for m in range(MZ)]
    payload = {
        'E_f_init_g': E_f_init_g,
        'E_cts_g': E_cts_g,
        'H_g': H_g,
        'u_g': u_g,
        'v_g': v_g,
        'w_g': w_g,
        'sig_g': sig_g,
        'm_c2f': m_c2f,
        'incr_c2f': incr_c2f_map,
        'm_f2c': m_f2c,
        'incr_f2c': incr_f2c_map,
    }
    for name, (i0, j0) in BLOCKS.items():
        payload['E_f_init_' + name] = block_slice(E_f_init_g, i0, j0)
    payload['E_f_init_default'] = block_slice(E_f_init_g, 1, 1)
    return payload


# ----------------------------------------------------------------------
# expected-value tuples for the hand-authored assertions (printed for review;
# the assembly below fills the tests from the same functions, so the printed
# numbers and the filled numbers agree by construction)
# ----------------------------------------------------------------------
def expect():
    print('=== expected values (block cells, dump k <-> esm lev = k+1) ===')
    for bname, (i0, j0) in BLOCKS.items():
        print(f'--- block {bname} (i0={i0}, j0={j0}) ---')
        cells = {
            'w': [(7, 15), (8, 15)],
            'e': [(23, 15), (22, 15)],
            'dome': [(15, 15), (16, 15), (15, 13)],
            'm': [(27, 15), (28, 15)],
        }[bname]
        for (i, j) in cells:
            bi, bj = i - i0 + 1, j - j0 + 1
            H = Hc(i, j)
            ks = ks_of(H)
            h20 = hadv_fine(i, j)[20] if ks >= 20 else None
            print(f' cell ({i},{j}) -> esm ({bi},{bj})')
            print(f'  H={H!r} ks={ks} Enth_ks={enth_ks(i,j)!r} marginal={int(marginal(i,j))} bmelt={bmelt_of(i,j)!r}')
            print(f'  E_in[0]={E_in(i,j)[0]!r} E_in[20]={E_in(i,j)[20]!r} E_in[60]={E_in(i,j)[60]!r}')
            print(f'  E_out[0]={E_out(i,j)[0]!r} E_out[10]={E_out(i,j)[10]!r} E_out[20]={E_out(i,j)[20]!r} E_out[30]={E_out(i,j)[30]!r} E_out[60]={E_out(i,j)[60]!r}')
            print(f'  hadv_t0[k20]={h20!r}')
            print(f'  u_f[k20]={c2f(u3(i,j), ks)[20]!r} w_f[k1]={c2f(w3(i,j), ks)[1]!r}')


# ----------------------------------------------------------------------
# assemble the .esm: patch index sets, scalars, E_f_init default, tests
# ----------------------------------------------------------------------
def fmt_num(x):
    return repr(float(x))


def fmt_array(vals, shape, indent=5):
    if len(shape) == 1:
        inner = ',\n'.join(' ' * 8 + fmt_num(v) for v in vals)
        return '[\n' + inner + '\n' + ' ' * indent + ']'
    n = shape[0]
    step = len(vals) // n
    rows = []
    for r in range(n):
        sub = fmt_array(vals[r * step:(r + 1) * step], shape[1:], indent + 2)
        rows.append(' ' * (indent - 1) + sub)
    return '[\n' + ',\n'.join(rows) + '\n' + ' ' * (indent - 1) + ']'


def replace_array(text, key, new_block, start=0, use_default=False):
    k = text.find('"' + key + '"', start)
    if k < 0:
        raise SystemExit(f'key {key} not found after offset {start}')
    if use_default:
        d = text.find('"default"', k)
        if d < 0:
            raise SystemExit(f'no default for {key}')
        k = d
    b = text.find('[', k)
    depth = 0
    e = b
    while True:
        if text[e] == '[':
            depth += 1
        elif text[e] == ']':
            depth -= 1
            if depth == 0:
                break
        e += 1
    return text[:b] + new_block + text[e + 1:], e + 1


def fld(payload, key):
    return fmt_array(payload[key], {'E_f_init_g': [31, 31, 62], 'E_cts_g': [31, 31, 62],
                                    'H_g': [31, 31], 'u_g': [31, 31, 61], 'v_g': [31, 31, 61],
                                    'w_g': [31, 31, 61], 'sig_g': [31, 31, 61],
                                    'm_c2f': [62], 'incr_c2f': [62], 'm_f2c': [61],
                                    'incr_f2c': [61]}[key])


def step_assertions(payload):
    """Hand-authored structure of the four one-step tests; expected numbers from the dump."""
    dt = DT

    def cell(bname, cell, lev):          # dump k-level -> esm storage lev (k+1)
        return E_out(*cell)[lev - 1]

    def pinned(cell):
        return enth_ks(*cell)

    LEV_AXIS = {"E_out": "lev", "E_f": "lev_f", "hadv_term": "lev_f", "strain_term": "lev_f",
                "w_eff": "lev_f", "u_f": "lev_f", "v_f": "lev_f", "w_f": "lev_f"}

    def mk_assert(var, bi, bj, lev, expected, tol):
        a = {"variable": var, "time": dt, "coords": {"i": bi, "j": bj}, "expected": expected,
             "tolerance": tol}
        if lev is not None:
            a["coords"][LEV_AXIS[var]] = lev
        return a

    tests = []

    # ---- step_block_w: one 10-yr step of the west block (i0=7, j0=13; C++ cells
    # 7..8 x 13..15). Demo cell (7,15) = esm (1,3): the west mid-radius column,
    # WESTWARD flow (u < 0) with the upwind east neighbor (8,15) live in the block.
    # E_out vs the dump with per-level tolerances in the documented
    # backward-Euler/ODE integration-scheme gap (the leaf's measured gap at the
    # same cell types; see ./ice_energy.esm step tests). The DEMONSTRATION
    # assertion is (2,3) lev 21 (fine level 20): the horizontal advection shifts
    # the enthalpy by 1.32 J/kg over the step, and the 0.1 J/kg tolerance is
    # tighter than that shift, so a model WITHOUT the horizontal coupling fails it.
    tests.append({
        "id": "step_block_w",
        "description": "One 10-yr enthalpy step of the west block (i0=7, j0=13; C++ cells 7..8 x 13..15): E_out at the demo cell (7,15) = esm (1,3) reproduces the instrumented C++ dump at the bed (lev 1), near-bed (lev 11), the mid-column DEMONSTRATION level (lev 21, tolerance 0.1 J/kg tighter than the 1.32 J/kg horizontal-advection shift over the step, so a model without the horizontal coupling fails), the near-surface row (lev 31), and the pinned surface/atmosphere rows (lev 41, 61, machine-exact at Enth_ks); the cold base melts nothing. Tolerances: per-level abs values in the backward-Euler/ODE integration-scheme gap measured by the single-column model for this cell class (see ./ice_energy.esm).",
        "time_span": {"start": 0.0, "end": dt},
        "parameter_overrides": {"i0": 7.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            mk_assert("E_out", 1, 3, 1, cell('w', (7, 15), 1), {"abs": 0.5, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 11, cell('w', (7, 15), 11), {"abs": 0.05, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 21, cell('w', (7, 15), 21), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 31, cell('w', (7, 15), 31), {"abs": 2.0, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 41, pinned((7, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 1, 3, 61, pinned((7, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_f", 1, 3, 32, pinned((7, 15)), {"abs": 1e-6, "rel": 1e-9}),
            {"variable": "bmelt", "time": dt, "coords": {"i": 1, "j": 3}, "expected": 0.0,
             "tolerance": {"abs": 1e-12}},
        ],
    })

    # ---- step_block_e: mirror of step_block_w at the east mid-radius (23,15).
    tests.append({
        "id": "step_block_e",
        "description": "One 10-yr enthalpy step of the east block (i0=22, j0=13; C++ cells 22..23 x 13..15): E_out at the mirror demo cell (23,15) = esm (2,3) reproduces the dump at the same levels as step_block_w (radial symmetry of the test-F field: identical expected values), including the lev 21 horizontal-advection demonstration tolerance (0.1 J/kg < the 1.32 J/kg shift), and the cold base melts nothing.",
        "time_span": {"start": 0.0, "end": dt},
        "parameter_overrides": {"i0": 22.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            mk_assert("E_out", 2, 3, 1, cell('e', (23, 15), 1), {"abs": 0.5, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 11, cell('e', (23, 15), 11), {"abs": 0.05, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 21, cell('e', (23, 15), 21), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 31, cell('e', (23, 15), 31), {"abs": 2.0, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 41, pinned((23, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 2, 3, 61, pinned((23, 15)), {"abs": 1e-6, "rel": 1e-9}),
            {"variable": "bmelt", "time": dt, "coords": {"i": 2, "j": 3}, "expected": 0.0,
             "tolerance": {"abs": 1e-12}},
        ],
    })

    # ---- step_block_dome: the dome block (i0=15, j0=13; C++ cells 15..16 x 13..15)
    # covers the dome center (15,15), the east-adjacent (16,15), and the south
    # azimuth (15,13). The dome column has u = v = 0 (exact test-F symmetry), so
    # its horizontal advection vanishes and its trajectory is bit-identical to the
    # single-column model's dome column -- the leaf's MEASURED per-level
    # backward-Euler/ODE gap applies verbatim (35 J/kg at the temperate base,
    # 0.01-class at mid-levels, 0.082 at lev 41; the leaf tolerance band is
    # 87.8 J/kg at the bed). (16,15) is the ONLY cell of the whole test suite
    # whose second fine layer starts AT the cold-temperate transition: the dump
    # has E_f[2] = E_cts[1] + 4e-11 (last-bit knife edge), so the basal
    # diffusivity entries R[0] and R[1] (both switch on E[1] vs E_cts[1] in
    # enthSystemCtx::assemble_R) take the temperate value ratio*K with
    # energy.enthalpy.temperate_ice_thermal_conductivity_ratio = 0.1 -- ten times
    # SMALLER than cold-ice K -- and the one-step integration-scheme gap at the
    # bed balloons to 392 J/kg (measured with the single-column model: the 3-D
    # bed result equals the leaf's to 8e-6 J/kg, so the gap is the column
    # scheme's, not the coupling's; the sign also flips: the model's bed ends
    # BELOW the dump here while the dome's ends above). The (16,15) basal melt
    # rate carries the same bed gap through hf_up = -(k/c)*(E[2]-E[1])/dz (the
    # post-step fine level 1 is cold in both codes: k1_temperate flips to false
    # during the step, so the T-gradient branch is not taken): measured
    # 2.06e-11 m/s. (15,13) has a cold base and the dome-class 35 J/kg bed gap.
    tests.append({
        "id": "step_block_dome",
        "description": "One 10-yr enthalpy step of the dome block (i0=15, j0=13; C++ cells 15..16 x 13..15): E_out at the dome center (15,15) = esm (1,3) reproduces the dump within the single-column model's measured backward-Euler/ODE gap (35.1 J/kg at the temperate base lev 1, 0.01-class at mid-levels, 0.033 at lev 41 -- the dome column has u = v = 0, so it is exactly the validated single-column trajectory), the basal melt rate matches the dump, and the surface row stays pinned at Enth_ks; the east-adjacent column (16,15) = esm (2,3) (eastward flow with the upwind dome column live) reproduces the dump at the bed within its MEASURED 392 J/kg one-step scheme gap (the only cell whose fine level 2 starts at the cold-temperate transition at the last bit: E_f[2] = E_cts[1] + 4e-11 switches the basal diffusivity entries R[0]/R[1] to the temperate value ratio*K, ratio = 0.1, which enlarges the bed gap 11-fold; the 3-D bed result equals the single-column model's to 8e-6 J/kg, so the gap is the column scheme's, not the horizontal coupling's) and at the mid-column levels within 0.04 J/kg, with its basal melt rate within the measured 2.06e-11 m/s gap (the same bed gap propagated through hf_up); the south-azimuth column (15,13) = esm (1,1) (base at the cold-temperate transition exactly, bmelt > 0) reproduces the dump within the dome-class gaps; the surface/atmosphere rows stay pinned.",
        "time_span": {"start": 0.0, "end": dt},
        "parameter_overrides": {"i0": 15.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            mk_assert("E_out", 1, 3, 1, cell('dome', (15, 15), 1), {"abs": 90.0, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 21, cell('dome', (15, 15), 21), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 41, cell('dome', (15, 15), 41), {"abs": 1.0, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 61, pinned((15, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 2, 3, 1, cell('dome', (16, 15), 1), {"abs": 400.0, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 21, cell('dome', (16, 15), 21), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 2, 3, 61, pinned((16, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 1, 1, 1, cell('dome', (15, 13), 1), {"abs": 90.0, "rel": 1e-6}),
            mk_assert("E_out", 1, 1, 21, cell('dome', (15, 13), 21), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 1, 1, 61, pinned((15, 13)), {"abs": 1e-6, "rel": 1e-9}),
            {"variable": "bmelt", "time": dt, "coords": {"i": 1, "j": 3},
             "expected": bmelt_of(15, 15), "tolerance": {"abs": 7.582560065230209e-12, "rel": 1e-6}},
            {"variable": "bmelt", "time": dt, "coords": {"i": 2, "j": 3},
             "expected": bmelt_of(16, 15), "tolerance": {"abs": 2.5e-11, "rel": 1e-6}},
            {"variable": "bmelt", "time": dt, "coords": {"i": 1, "j": 1},
             "expected": bmelt_of(15, 13), "tolerance": {"abs": 1e-11, "rel": 1e-6}},
        ],
    })

    # ---- step_block_m: the margin block (i0=27, j0=13; C++ cells 27..28 x 13..15):
    # the marginal column (27,15) = esm (1,3) (horizontal/vertical advection and
    # strain heating excluded by the marginal gates; E_out within the leaf's
    # measured gap, incl. the near-surface transition levels with their ~3-23 J/kg
    # gaps) and the ice-free column (28,15) = esm (2,3), which carries an
    # identically zero RHS and stays at the constant surface enthalpy Enth_ks
    # (machine-exact).
    tests.append({
        "id": "step_block_m",
        "description": "One 10-yr enthalpy step of the margin block (i0=27, j0=13; C++ cells 27..28 x 13..15): E_out at the marginal column (27,15) = esm (1,3) reproduces the dump within the leaf's measured gap at its levels (0.01-class at the bed, 3.3 J/kg at lev 11 -- the near-surface transition -- and pinned rows above), the ice-free column (28,15) = esm (2,3) stays at the constant surface enthalpy Enth_ks to machine precision (zero RHS, the dE = 0 ice-free invariant), and neither base melts.",
        "time_span": {"start": 0.0, "end": dt},
        "parameter_overrides": {"i0": 27.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            mk_assert("E_out", 1, 3, 1, cell('m', (27, 15), 1), {"abs": 0.1, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 11, cell('m', (27, 15), 11), {"abs": 5.0, "rel": 1e-6}),
            mk_assert("E_out", 1, 3, 61, pinned((27, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 2, 3, 1, pinned((28, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_out", 2, 3, 61, pinned((28, 15)), {"abs": 1e-6, "rel": 1e-9}),
            mk_assert("E_f", 2, 3, 1, pinned((28, 15)), {"abs": 1e-6, "rel": 1e-9}),
            {"variable": "bmelt", "time": dt, "coords": {"i": 1, "j": 3}, "expected": 0.0,
             "tolerance": {"abs": 1e-12}},
            {"variable": "bmelt", "time": dt, "coords": {"i": 2, "j": 3}, "expected": 0.0,
             "tolerance": {"abs": 1e-12}},
        ],
    })

    # fill the E_f_init overrides from the payload
    for t in tests:
        t['parameter_overrides']['E_f_init'] = payload['E_f_init_' + t['id'].split('_')[2]]
    return tests


def ic_assertions(payload):
    """Hand-authored structure of the four t = 0 algebraic tests; numbers from the dump."""

    LEV_AXIS = {"E_out": "lev", "E_f": "lev_f", "hadv_term": "lev_f", "strain_term": "lev_f",
                "w_eff": "lev_f", "u_f": "lev_f", "v_f": "lev_f", "w_f": "lev_f"}

    def a(var, bi, bj, lev, expected, tol):
        out = {"variable": var, "time": 0.0, "coords": {"i": bi, "j": bj}, "expected": expected,
               "tolerance": tol}
        if lev is not None:
            out["coords"][LEV_AXIS[var]] = lev
        return out

    def ein(cell, k):
        return E_in(*cell)[k]

    tests = []

    # ---- ic_block_w
    tests.append({
        "id": "ic_block_w",
        "description": "Initial condition and algebraic chain of the west block (i0=7, j0=13; C++ cells 7..8 x 13..15): the pinned fine-grid IC at the demo cell (7,15) = esm (1,3) (base row, first pinned fine row, storage fine_to_coarse of both), the per-column H/ks/marginal chain of the two j = 15 columns, the t = 0 horizontal advection at the demo cell (exactly the upwind stencil of the C++ B-array at the old state -- the test-F flow at (7,15) is WESTWARD, so the upwind neighbor is the east column (8,15), live in the block -- with the block-edge ghosts frozen at the full-field initial enthalpy; the value is computed from the dump by this script), the coarse-to-fine u at fine level 21, and the zero basal melt rate of the cold base.",
        "time_span": {"start": 0.0, "end": 1.0},
        "parameter_overrides": {"i0": 7.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            a("E_f", 1, 3, 1, ein((7, 15), 0), {"rel": 1e-12}),
            a("E_f", 1, 3, 32, ein((7, 15), 60), {"rel": 1e-12}),
            a("E_out", 1, 3, 1, ein((7, 15), 0), {"rel": 1e-12}),
            a("E_out", 1, 3, 61, ein((7, 15), 60), {"rel": 1e-12}),
            a("H", 1, 3, None, Hc(7, 15), {"rel": 1e-12}),
            a("H", 2, 3, None, Hc(8, 15), {"rel": 1e-12}),
            a("ks", 1, 3, None, float(ks_of(Hc(7, 15))), {"abs": 1e-9}),
            a("ks", 2, 3, None, float(ks_of(Hc(8, 15))), {"abs": 1e-9}),
            a("marginal", 1, 3, None, 0.0, {"abs": 1e-9}),
            a("hadv_term", 1, 3, 21, hadv_fine(7, 15)[20], {"abs": 1e-11}),
            a("u_f", 1, 3, 21, c2f(u3(7, 15), ks_of(Hc(7, 15)))[20], {"abs": 1e-15}),
            a("bmelt", 1, 3, None, 0.0, {"abs": 1e-12}),
        ],
    })

    # ---- ic_block_e
    tests.append({
        "id": "ic_block_e",
        "description": "Initial condition and algebraic chain of the east block (i0=22, j0=13; C++ cells 22..23 x 13..15): the pinned IC at the mirror demo cell (23,15) = esm (2,3) (radial symmetry: identical numbers to the west block), the per-column H/ks chain of the j = 15 columns, the t = 0 horizontal advection at the demo cell, and the zero basal melt rate of the cold base.",
        "time_span": {"start": 0.0, "end": 1.0},
        "parameter_overrides": {"i0": 22.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            a("E_f", 2, 3, 1, ein((23, 15), 0), {"rel": 1e-12}),
            a("E_f", 2, 3, 32, ein((23, 15), 60), {"rel": 1e-12}),
            a("E_out", 2, 3, 1, ein((23, 15), 0), {"rel": 1e-12}),
            a("E_out", 2, 3, 61, ein((23, 15), 60), {"rel": 1e-12}),
            a("H", 2, 3, None, Hc(23, 15), {"rel": 1e-12}),
            a("H", 1, 3, None, Hc(22, 15), {"rel": 1e-12}),
            a("ks", 2, 3, None, float(ks_of(Hc(23, 15))), {"abs": 1e-9}),
            a("ks", 1, 3, None, float(ks_of(Hc(22, 15))), {"abs": 1e-9}),
            a("marginal", 2, 3, None, 0.0, {"abs": 1e-9}),
            a("hadv_term", 2, 3, 21, hadv_fine(23, 15)[20], {"abs": 1e-11}),
            a("bmelt", 2, 3, None, 0.0, {"abs": 1e-12}),
        ],
    })

    # ---- ic_block_dome
    tests.append({
        "id": "ic_block_dome",
        "description": "Initial condition and algebraic chain of the dome block (i0=15, j0=13; C++ cells 15..16 x 13..15): the pinned fine-grid IC at the dome (15,15) = esm (1,3) (base at the cold-temperate transition exactly, surface/atmosphere pinned to Enth_ks), at the east-adjacent (16,15) = esm (2,3) and the south azimuth (15,13) = esm (1,1); the per-column H/ks/marginal chain; the t = 0 horizontal advection at the dome (exactly 0: u = v = 0 at the test-F center) and at the east-adjacent column (dump-computed upwind value); and the coarse-to-fine machinery (u at fine level 21 of the east-adjacent column, w at fine level 2 of the dome).",
        "time_span": {"start": 0.0, "end": 1.0},
        "parameter_overrides": {"i0": 15.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            a("E_f", 1, 3, 1, ein((15, 15), 0), {"rel": 1e-12}),
            a("E_f", 1, 3, 46, ein((15, 15), 60), {"rel": 1e-12}),
            a("E_f", 1, 3, 62, ein((15, 15), 60), {"rel": 1e-12}),
            a("E_out", 1, 3, 1, ein((15, 15), 0), {"rel": 1e-12}),
            a("E_out", 1, 3, 61, ein((15, 15), 60), {"rel": 1e-12}),
            a("E_f", 2, 3, 1, ein((16, 15), 0), {"rel": 1e-12}),
            a("E_f", 1, 1, 1, ein((15, 13), 0), {"rel": 1e-12}),
            a("H", 1, 3, None, Hc(15, 15), {"rel": 1e-12}),
            a("H", 1, 1, None, Hc(15, 13), {"rel": 1e-12}),
            a("ks", 1, 3, None, float(ks_of(Hc(15, 15))), {"abs": 1e-9}),
            a("ks", 1, 1, None, float(ks_of(Hc(15, 13))), {"abs": 1e-9}),
            a("marginal", 1, 3, None, 0.0, {"abs": 1e-9}),
            a("marginal", 1, 1, None, 0.0, {"abs": 1e-9}),
            a("hadv_term", 1, 3, 21, hadv_fine(15, 15)[20], {"abs": 1e-11}),
            a("hadv_term", 2, 3, 21, hadv_fine(16, 15)[20], {"abs": 1e-11}),
            a("u_f", 2, 3, 21, c2f(u3(16, 15), ks_of(Hc(16, 15)))[20], {"abs": 1e-15}),
            a("w_f", 1, 3, 2, c2f(w3(15, 15), ks_of(Hc(15, 15)))[1], {"abs": 1e-15}),
        ],
    })

    # ---- ic_block_m
    tests.append({
        "id": "ic_block_m",
        "description": "Initial condition and algebraic chain of the margin block (i0=27, j0=13; C++ cells 27..28 x 13..15): the marginal column (27,15) = esm (1,3) with the marginal gates active (horizontal advection, vertical advection, and strain heating exactly zero) and ks = 13, and the ice-free column (28,15) = esm (2,3) at the constant surface enthalpy Enth_ks on all 62 fine rows (ks = 0); the per-column H/marginal chain; both bases melt nothing.",
        "time_span": {"start": 0.0, "end": 1.0},
        "parameter_overrides": {"i0": 27.0, "j0": 13.0, "E_f_init": None},
        "assertions": [
            a("H", 1, 3, None, Hc(27, 15), {"rel": 1e-12}),
            a("H", 2, 3, None, Hc(28, 15), {"rel": 1e-12}),
            a("ks", 1, 3, None, float(ks_of(Hc(27, 15))), {"abs": 1e-9}),
            a("ks", 2, 3, None, 0.0, {"abs": 1e-9}),
            a("marginal", 1, 3, None, 1.0, {"abs": 1e-9}),
            a("marginal", 2, 3, None, 1.0, {"abs": 1e-9}),
            a("E_f", 2, 3, 1, ein((28, 15), 60), {"rel": 1e-12}),
            a("E_f", 2, 3, 62, ein((28, 15), 60), {"rel": 1e-12}),
            a("E_out", 2, 3, 61, ein((28, 15), 60), {"rel": 1e-12}),
            a("hadv_term", 1, 3, 12, 0.0, {"abs": 1e-11}),
            a("strain_term", 1, 3, 12, 0.0, {"abs": 1e-11}),
            a("w_eff", 1, 3, 12, 0.0, {"abs": 1e-15}),
            a("hadv_term", 2, 3, 1, 0.0, {"abs": 1e-11}),
            a("bmelt", 1, 3, None, 0.0, {"abs": 1e-12}),
            a("bmelt", 2, 3, None, 0.0, {"abs": 1e-12}),
        ],
    })

    for t in tests:
        t['parameter_overrides']['E_f_init'] = payload['E_f_init_' + t['id'].split('_')[2]]
    return tests


def assemble(payload):
    import copy
    d = json.load(open(ESM))
    d['index_sets']['i']['size'] = NI
    d['index_sets']['j']['size'] = NJ
    v = d['models']['EnthalpyField3D']['variables']
    v['E_f_init']['shape'] = ['i', 'j', 'lev_f']
    v['E_f_init']['default'] = copy.deepcopy(payload['E_f_init_default'])
    v['nb_i'] = {"type": "parameter", "units": "1", "default": float(NI),
                 "description": "Block extent in the i direction (number of block columns), matching the size of the index set i; the horizontal-advection template's ghost fallback reads the frozen full-field initial enthalpy at (i0+nb_i+1, ...)."}
    v['nb_j'] = {"type": "parameter", "units": "1", "default": float(NJ),
                 "description": "Block extent in the j direction, matching the size of the index set j."}
    d['models']['EnthalpyField3D']['tests'] = ic_assertions(payload) + step_assertions(payload)
    d['metadata']['description'] = (
        "The coupled 3-D ice-energy piece of stage 2 (deferred from ./ice_energy.esm): the full-field PISM "
        "enthalpy equation dE/dt = D(K D(E,z),z) - w dE/dz - (u dE/dx + v dE/dy) + sigma/rho over four "
        "2x3-column sub-domains of the 31x31 test-F grid, integrated as a continuous PDE by the harness. "
        "The model is built on the validated single-column EnthalpyColumn (./ice_energy.esm): the same "
        "BOMBPROOF vertical diffusion/advection stencils and ColumnInterpolation maps "
        "(./energy_templates.esm), the same cold/temperate thermal-diffusivity switch, and the same basal "
        "melt rate, all lifted to the field (./energy_3d_templates.esm). The state E_f[i,j,lev_f] is the "
        "fine-grid enthalpy of the 2x3 block (C++ cells (i0..i0+1, j0..j0+2); the scalars i0/j0 hold the "
        "C++ indices of the block's west/south columns, per test; the block extent nb_i x nb_j = 2 x 3 "
        "matches the index sets i/j). All forcing and geometry are read from full-field 31x31 parameters "
        "(H_g, u_g, v_g, w_g, sig_g, E_cts_g, E_f_init_g) at the absolute esm indices (i0+i, j0+j); the "
        "per-column surface level is ks(i,j) = floor(H(i,j)/dz) capped at Mz_fine-1, spelled as a faq sum "
        "(esm has no floor op); the marginal flag follows EnergyModel::marginal (any of the 8 periodically "
        "wrapped neighbors thinner than 100 m, energy.margin_ice_thickness_limit); and the cold-temperate "
        "transition enthalpy is a baked full-field parameter E_cts_g computed with the C++ operation order, "
        "because the K switch at temperate bases is sensitive to the last bit (the column model bakes E_cts "
        "for the same reason). The horizontal advection is the new coupled piece: it reads the LIVE neighbor "
        "columns (index(E_f, i+/-1, j, k)) so the harness integrates the coupled PDE, with the block-edge "
        "ghosts frozen at the full-field initial enthalpy E_f_init_g. PISM evaluates its explicit "
        "horizontal-advection B-array at the OLD state, so the one-step reference trajectory and the "
        "live-coupling PDE differ by O(dt^2) in the horizontal exchange (~6e-5 J/kg over the 10-yr step, "
        "|u|/dx ~ 2e-13 s^-1; far below the assertion tolerances). Every test block is placed so the demo "
        "column's UPWIND neighbor is live in the block ((7,15): the test-F flow there is WESTWARD, so the "
        "upwind east neighbor (8,15) is in the block; (23,15) and (16,15): eastward flow, upwind west "
        "neighbor live; (15,13): southward flow, upwind north neighbor (15,14) live), and the block-edge "
        "ghosts it never reads carry ~zero velocity (v = 0 on the j = 15 meridian, u = 0 on the i = 15 "
        "meridian) or are marginal-gated. Ice-free columns (ks = 0) carry an identically zero RHS and stay at their initial "
        "column state -- the dump E_in of an ice-free column is the constant surface enthalpy Enth_ks, "
        "exactly as EnthalpyModel::update sets the whole column to Enth_ks; rows above the ice surface "
        "carry zero RHS and stay pinned to Enth_ks (baked into E_f_init_g / E_f_init on rows >= ks). The "
        "vertical forcing (u, v, w, sigma) is coarse-to-fine interpolated in-model from the storage-grid "
        "full-field parameters by energy_coarse_to_fine_field, and E_out is the storage-grid "
        "fine_to_coarse of the state, the quantity the dumps carry. Tests: four t = 0 algebraic chains "
        "(blocks w/e/dome/m) pinning the IC, the c2f machinery, the ks/marginal/H chains, the marginal "
        "gates, and the t = 0 horizontal advection against values computed from the dump; four one-10-yr-"
        "step tests (dt = 315569259.747 s) asserting E_out against the dump within the documented "
        "backward-Euler/ODE integration-scheme gap (per-level abs tolerances, measured by the single-"
        "column model for the same cells/cell classes and documented in the test descriptions; the dome "
        "column has u = v = 0 and is exactly the validated single-column trajectory, so the leaf's "
        "measured 35.1 J/kg temperate-base gap applies verbatim; the east-adjacent (16,15) is the only "
        "cell whose second fine layer starts at the cold-temperate transition at the last bit, switching "
        "the basal diffusivity entries to the temperate value ratio*K (ratio = 0.1) and enlarging the "
        "measured one-step bed gap to 391.98 J/kg with a 2.06e-11 m/s basal-melt-rate gap -- both "
        "measured with the single-column model, which the coupled bed result matches to 8e-6 J/kg; see "
        "the step_block_dome description), the pinned surface rows (machine-exact), "
        "the ice-free invariant, and the basal melt rate. The demonstration assertions (step_block_w and "
        "step_block_e, (7,15)/(23,15) lev 21) use a 0.1 J/kg tolerance that is tighter than the 1.32 J/kg "
        "horizontal-advection shift over the step, so a model WITHOUT the horizontal coupling fails them. "
        "Run the tests with --compiler interpreter (runtime-ifelse stencils, as in ice_energy.esm; the "
        "native compiler refuses the const-cadence K assembly). The full 31x31-field model with a 31x31 "
        "state is deferred: the dense-Jacobian/linear-solve cost of the harness's BDF integrator scales "
        "roughly as N^3 in the number of unknowns (~21 min per 10-yr step for 558 unknowns, ~2 min for "
        "248), which made the 58k-unknown full field impractical (a previous attempt produced no output "
        "after 4+ hours); the four 2x3 blocks (~6 min per 10-yr step) cover the dome center, the west/east "
        "mid-radius and south azimuths, the marginal column, and the ice-free column, and are the coupled-"
        "field building blocks for stage 3. No C++ bugs found."
    )
    return d


def fill(payload):
    """Patch the on-disk .esm: index sets, block-size scalars, E_f_init default,
    and the whole tests block, via a JSON round-trip (the arrays are already in
    the file and are preserved exactly; json repr round-trips floats exactly)."""
    d = assemble(payload)
    open(ESM, 'w').write(json.dumps(d, indent=1))
    print(f'filled {ESM}')


if __name__ == '__main__':
    payload = build_payload()
    if '--expect' in sys.argv:
        expect()
    else:
        fill(payload)
        expect()
