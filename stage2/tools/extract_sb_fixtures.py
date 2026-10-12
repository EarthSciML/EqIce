#!/usr/bin/env python3
"""Fill the data fields of stage2/stress_balance.esm from the Stage-1 stress-balance dumps.

This script performs ONLY numeric data extraction (AGENTS.md: "Scripts MAY be used
to extract numeric test tuples from C++ dumps into an existing hand-authored tests
block"). It fills:

- the StressBalance default parameters (the sia-run cell (20,15) family):
  center-column E/H/ks, the eight neighbor 3D-velocity columns, tauc;
- the test-1 (sia_cell_20_15) surface-gradient overrides (SIAVel.h_*_* and
  SIAFlux.h_* at the four staggered points around (20,15)) and the dumped
  assertion values (u3/v3/w3/sigma at the asserted levels, q_x);
- the test-2 (sia_margin_27_15) parameter overrides (the cell (27,15) family:
  the four staggered-column thk/ks/E/alpha/theta, the surface gradients, the
  center column, the eight neighbor columns, mask_e, tauc) and its dumped
  assertion values.

The equations, structure, and factoring are hand-authored and are not touched.
Sources: stage1/dumps/sia/columns.csv (cell columns), delta.csv (staggered
diffusivity-chain columns; note the dump prints only the ks+1 in-ice levels, so
the E arrays are zero-padded to 61), staggered_D.csv (h_x, h_y, alpha),
diffusive_flux.csv (q at (20,15,0) and (27,15,0)), basal_strength/tauc_map.csv.
"""
import json
import re
import sys

ESM = 'stage2/stress_balance.esm'
COLS = 'stage1/dumps/sia/columns.csv'
DELTA = 'stage1/dumps/sia/delta.csv'
STAG = 'stage1/dumps/sia/staggered_D.csv'
FLUX = 'stage1/dumps/sia/diffusive_flux.csv'
TAUC = 'stage1/dumps/basal_strength/tauc_map.csv'

NLEV = 61


def load_csv(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            rows.append(line.rstrip('\n').split(','))
    return rows


def fmt_num(x):
    return repr(float(x))


def fmt_arr(values):
    return '[' + ', '.join(fmt_num(v) for v in values) + ']'


# ---------------------------------------------------------------- dump tables
CHDR = ['i', 'j', 'x', 'y', 'r', 'H', 'surface', 'bed', 'mask', 'k', 'z', 'E',
        'T_exact', 'u3', 'v3', 'w3', 'sigma3', 'sigma_exact', 'u_exact',
        'v_exact', 'w_exact']
cells = {}
for r in load_csv(COLS):
    cells.setdefault((int(r[0]), int(r[1])), []).append(dict(zip(CHDR, r)))


def cell_arr(cell, name):
    col = cells[cell]
    assert len(col) == NLEV, f'cell {cell} has {len(col)} rows, expected {NLEV}'
    return [float(row[name]) for row in col]


def cell_val(cell, name):
    return float(cells[cell][0][name])


DHDR = ['i', 'j', 'o', 'k', 'z', 'depth', 'pressure', 'E', 'stress', 'flow',
        'delta', 'D', 'theta']
pts = {}
for r in load_csv(DELTA):
    pts.setdefault((int(r[0]), int(r[1]), int(r[2])), []).append(dict(zip(DHDR, r)))


def stag_point(point):
    col = pts[point]
    ks = int(col[-1]['k'])
    thk = float(col[0]['depth'])  # depth at k=0 = full thickness at the point
    E = [float(row['E']) for row in col] + [0.0] * (NLEV - len(col))
    theta = float(col[0]['theta'])
    assert len(E) == NLEV
    return ks, thk, E, theta


stagd = {(int(r[0]), int(r[1]), int(r[2])): (float(r[5]), float(r[6]), float(r[7]))
         for r in load_csv(STAG)}
fluxq = {(int(r[0]), int(r[1]), int(r[2])): (float(r[8]), float(r[9]))
         for r in load_csv(FLUX)}
tauc_map = {(int(r[0]), int(r[1])): float(r[14]) for r in load_csv(TAUC)}

# sanity: q_x = -D*h_x at the two flux points
for p in [(20, 15, 0), (27, 15, 0)]:
    r = [row for row in load_csv(FLUX) if (int(row[0]), int(row[1]), int(row[2])) == p][0]
    assert abs(float(r[8]) - (-float(r[5]) * float(r[6]))) < 1e-18 * abs(float(r[8])) + 1e-30, p

# ------------------------------------------------------------------- scenarios
# scenario 1: the sia-run grounded mid-radius cell (20,15)
S1 = dict(center=(20, 15),
          nb=dict(w=(19, 15), e=(21, 15), s=(20, 14), n=(20, 16)),
          st=dict(e=(20, 15, 0), w=(19, 15, 0), n=(20, 15, 1), s=(20, 14, 1)))
# scenario 2: the margin cell (27,15) with the ice-free east neighbor (28,15)
S2 = dict(center=(27, 15),
          nb=dict(w=(26, 15), e=(28, 15), s=(27, 14), n=(27, 16)),
          st=dict(e=(27, 15, 0), w=(26, 15, 0), n=(27, 15, 1), s=(27, 14, 1)))


def scen_defaults(S):
    """Default-parameter values for the StressBalance model at scenario S."""
    # ks = kBelowHeight(H): the top in-ice level, read off as the last level with
    # a nonzero dumped strain heating (sigma3 is 0 above ks).
    ks = max(int(r['k']) for r in cells[S['center']] if float(r['sigma3']) != 0.0)
    out = {'E': cell_arr(S['center'], 'E'),
           'H': cell_val(S['center'], 'H'),
           'ks': float(ks),
           'tauc': tauc_map[S['center']]}
    # neighbor velocity columns (esm order: w, e, s, n)
    for d in 'uv':
        for side in 'wens':
            name = f'{d}3_{side}'
            out[name] = cell_arr(S['nb'][side], f'{d}3')
    return out


def stag_defaults(S):
    """Per-staggered-column and gradient values at scenario S (the SIA mounts)."""
    out = {}
    for side, point in S['st'].items():
        ks, thk, E, theta = stag_point(point)
        hx, hy, alpha = stagd[point]
        pref = 'SIA' + {'e': 'East', 'w': 'West', 'n': 'North', 's': 'South'}[side]
        out[f'{pref}.thk'] = thk
        out[f'{pref}.ks'] = float(ks)
        out[f'{pref}.E'] = E
        out[f'{pref}.alpha'] = alpha
        out[f'{pref}.theta'] = theta
        out[f'SIAVel.h_x_{side}'] = hx
        out[f'SIAVel.h_y_{side}'] = hy
    out['SIAFlux.h_x'] = stagd[S['st']['e']][0]
    out['SIAFlux.h_y'] = stagd[S['st']['e']][1]
    return out


# ------------------------------------------------------------------- text edits
def read_text():
    return open(ESM).read()


def find_key(text, start, key):
    """Index just past the colon of '"key": ' at or after start (same line)."""
    m = re.search(r'"' + re.escape(key) + r'":\s*', text[start:])
    if not m:
        raise SystemExit(f'key "{key}" not found after offset {start}')
    return start + m.end()


def value_end(text, start):
    """End (exclusive) of the value starting at start: a bracket-matched array
    or a scalar token ending at ',' or newline."""
    i = start
    if text[i] == '[':
        depth = 0
        while i < len(text):
            if text[i] == '[':
                depth += 1
            elif text[i] == ']':
                depth -= 1
                if depth == 0:
                    return i + 1
            i += 1
        raise SystemExit('unbalanced array')
    while i < len(text) and text[i] not in ',\n':
        i += 1
    return i


def obj_span(text, start):
    """Span of the object whose '{' is at start: (open, close_exclusive)."""
    assert text[start] == '{'
    depth = 0
    i = start
    while i < len(text):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1
        i += 1
    raise SystemExit('unbalanced object')


def replace_value(text, start, key, new):
    """Replace the value of '"key": ' at/after start; returns new text."""
    k = find_key(text, start, key)
    end = value_end(text, k)
    return text[:k] + new + text[end:]


def replace_default(text, var, new):
    """Replace the 'default' of variable var inside the StressBalance variables."""
    m = re.search(r'\n {8}"' + re.escape(var) + r'": \{\n', text)
    if not m:
        raise SystemExit(f'variable "{var}" not found')
    start = m.end()
    # variable object ends at the next 8-space key or the closing of variables
    nxt = re.search(r'\n {8}"', text[start:])
    end = start + (nxt.start() if nxt else len(text) - start)
    return replace_value(text, start, 'default', new)


def test_span(text, test_id):
    m = re.search(r'\n {10}"id": "' + re.escape(test_id) + r'",\n', text)
    if not m:
        raise SystemExit(f'test "{test_id}" not found')
    # the test object opens with '{' at 8-space indent on the line before the id
    braces = list(re.finditer(r'\n {8}\{', text[:m.start()]))
    if not braces:
        raise SystemExit(f'test "{test_id}" object brace not found')
    return obj_span(text, braces[-1].end() - 1)


def replace_override(text, test_id, key, new):
    """Replace the parameter-override value of key in the named test."""
    t0, t1 = test_span(text, test_id)
    po = re.search(r'\n {10}"parameter_overrides": \{', text[t0:t1])
    if not po:
        raise SystemExit(f'no parameter_overrides in test {test_id}')
    start = t0 + po.end() - 1
    po0, po1 = obj_span(text, start)
    k = find_key(text, po0, key)
    end = value_end(text, k)
    assert end <= po1, f'override {key} outside parameter_overrides'
    return text[:k] + new + text[end:]


def replace_expected(text, test_id, variable, coords, new):
    """Replace the 'expected' of the assertion (variable, coords) in the test."""
    t0, t1 = test_span(text, test_id)
    for m in re.finditer(r'\n {12}\{', text[t0:t1]):
        a0 = t0 + m.start()
        a1 = t0 + m.end()
        ao0, ao1 = obj_span(text, a1 - 1)
        seg = text[ao0:ao1]
        vm = re.search(r'"variable": "([^"]+)"', seg)
        cm = re.search(r'"coords": \{([^}]*)\}', seg)
        if not vm or vm.group(1) != variable:
            continue
        if cm is None and coords is None:
            pass
        elif cm is not None and coords is not None:
            got = dict(re.findall(r'"(\w+)":\s*([0-9.eE+-]+)', cm.group(1)))
            want = {str(k): str(v) for k, v in coords.items()}
            if got != want:
                continue
        else:
            continue
        return replace_value(text, ao0, 'expected', new)
    raise SystemExit(f'assertion {variable} {coords} not found in {test_id}')


def check_json(text):
    json.loads(text)
    print('JSON OK')


# ----------------------------------------------------------------------- fill
def main():
    d1 = scen_defaults(S1)
    s1 = stag_defaults(S1)
    d2 = scen_defaults(S2)
    s2 = stag_defaults(S2)

    text = read_text()

    # --- StressBalance defaults (scenario 1)
    for key, val in d1.items():
        text = replace_default(text, key, fmt_arr(val) if isinstance(val, list) else fmt_num(val))

    # --- test 1 overrides: the (20,15) staggered family (columns + gradients)
    for key, val in s1.items():
        text = replace_override(text, 'sia_cell_20_15', key,
                                fmt_arr(val) if isinstance(val, list) else fmt_num(val))

    # --- hybrid_wiring overrides: the same (20,15) staggered family, so the
    # SIA shear at sliding 0 equals the dumped (20,15) value
    for key, val in s1.items():
        text = replace_override(text, 'hybrid_wiring', key,
                                fmt_arr(val) if isinstance(val, list) else fmt_num(val))

    # --- test 2 overrides: the (27,15) family
    over2 = {}
    over2.update(s2)
    over2.update({'StressBalance.H': d2['H'],
                  'StressBalance.ks': d2['ks'],
                  'StressBalance.E': d2['E'],
                  'StressBalance.tauc': d2['tauc'],
                  'StressBalance.mask_e': cell_val(S2['nb']['e'], 'mask')})
    for d in 'uv':
        for side in 'wens':
            over2[f'StressBalance.{d}3_{side}'] = d2[f'{d}3_{side}']
    for key, val in over2.items():
        text = replace_override(text, 'sia_margin_27_15', key,
                                fmt_arr(val) if isinstance(val, list) else fmt_num(val))

    # --- dumped assertion values
    def fill_expected(text, test_id, variable, levels, column):
        c = S1['center'] if test_id == 'sia_cell_20_15' else S2['center']
        col = cell_arr(c, column)
        for lev in levels:
            text = replace_expected(text, test_id, variable, {'lev': lev}, fmt_num(col[lev - 1]))
        return text

    text = fill_expected(text, 'sia_cell_20_15', 'u3', [2, 11, 21, 38], 'u3')
    text = fill_expected(text, 'sia_cell_20_15', 'v3', [2], 'v3')
    text = fill_expected(text, 'sia_cell_20_15', 'w3', [2, 11, 21, 38, 61], 'w3')
    text = fill_expected(text, 'sia_cell_20_15', 'sigma', [2, 11, 21, 38], 'sigma3')
    text = replace_expected(text, 'sia_cell_20_15', 'q_x', None, fmt_num(fluxq[(20, 15, 0)][0]))

    text = fill_expected(text, 'sia_margin_27_15', 'u3', [2, 11, 14], 'u3')
    text = fill_expected(text, 'sia_margin_27_15', 'v3', [2], 'v3')
    text = fill_expected(text, 'sia_margin_27_15', 'w3', [2, 11, 14, 61], 'w3')
    text = fill_expected(text, 'sia_margin_27_15', 'sigma', [2, 14], 'sigma3')
    text = replace_expected(text, 'sia_margin_27_15', 'q_x', None, fmt_num(fluxq[(27, 15, 0)][0]))

    check_json(text)
    open(ESM, 'w').write(text)
    print('filled stress_balance.esm:')
    print('  StressBalance defaults:', ', '.join(d1))
    print('  test 1 overrides:', ', '.join(s1))
    print('  test 2 overrides:', ', '.join(over2))
    print('  hybrid_wiring overrides:', ', '.join(s1))
    print('  assertions: 20 (test 1) + 15 (test 2) expected values')


if __name__ == '__main__':
    main()
