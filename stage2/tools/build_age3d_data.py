#!/usr/bin/env python3
"""Fill the data fields and the inline tests block of stage2/ice_age_3d.esm
from the Stage-1 age dumps.

Numeric-tuple extraction only (AGENTS.md allows scripts for this; the
equations, structure, and factoring of the component are hand-authored in
stage2/ice_age_3d.esm). The script replaces the placeholder "default" arrays
of the block-local parameters u3_b, v3_b, w3_b, and ks_b with the dome-block
slice of the dumped fields (stage1/dumps/age/columns.csv columns u3, v3, w3 at
storage level k, and the fine-grid surface level ks = floor(H/dz)+1 derived
from the dumped thickness H with dz = 4000/61, capped at 62), and generates
the inline tests block: one step test per 7x7 block, with parameter_overrides
supplying that block's slice of u3_b/v3_b/w3_b/ks_b, asserting age_storage
after one dt step against the dumped age_out at hand-picked columns and
storage levels, plus pins of the model's own measured hadv at the end of the
step.

Conventions: the esm grid is 1-based (esm cell (i,j) = C++ cell (i-1,j-1));
the dump rows carry the 1-based C++ (i, j); the nested arrays in the file are
indexed by the dump (i, j) directly, so a 7x7 block at full-field origin
(i0, j0) holds block cell (bi,bj) = C++ cell (i0+bi-1, j0+bj-1) at array
position [i0+bi-1][j0+bj-1]; assertion coords are 1-based index positions
(esm-spec 6.6.5 convention 1), at block-local coordinates.

The velocity fields are stored as 63-entry padded columns (index set lev_p,
one more entry than the 61-level storage grid) so the fine-grid bracket reads
of age_coarse_to_fine are plain in-bounds gathers at every fine level k
(entry 0 = storage level 1, entries 1..61 = storage levels 1..61, entry 62 =
storage level 61; esm 1-based: entry k = storage level max(k-1, 1) for
k = 1..62, entry 63 = storage level 61).

The age_out assertions are judged at rel 1e-3 (deep and mid ice) and rel 1e-2
(the storage levels adjacent to each column's surface, where the C++ backward-
Euler-vs-exact-ODE vertical scheme difference of the single-column component
measured up to ~6e-3, and where the integrated horizontal contribution of the
developing age field is largest) and abs 1e-6 s (the first above-ice level and
the ice-free columns, exactly zero); the hadv pins are the model's own
measured horizontal-advection term at the end of the step (the C++ evaluates
it at the old time, where the whole-domain age-zero initial field freezes it
at zero; the ODE harness integrates it continuously, so the pinned values are
the model's own measured horizontal contribution, rel 1e-2 around the
measured values). The measured hadv literals are read from the JSON file
passed with --hadv-json (written from a first run's -v output); without it,
no hadv pins are generated.
"""
import argparse
import json
import math
import re

ESM = 'stage2/ice_age_3d.esm'
COLS = 'stage1/dumps/age/columns.csv'

DT = 315569259.747
DZ_FINE = 4000.0 / 61.0

# One entry per test block: (id, i0, j0, description, [columns...]) where each
# column is (block i, block j, dump i, dump j, ks, label, surface-adjacent
# storage levels judged at rel 1e-2, first above-ice level asserted at abs 1e-6
# or None for ice-free columns). The blocks are placed so every asserted
# column is interior (or, for the ice-free and thinnest-margin columns where
# the ghost error is zero or negligible, at an edge).
BLOCKS = [
    ('block_dome', 13, 13, 'the dome block (i0 = 13, j0 = 13; C++ cells 13..19 x 13..19)', [
        (3, 3, 15, 15, 46, 'dome_column', [44, 45], 46),
    ]),
    ('block_east', 22, 13, 'the east block (i0 = 22, j0 = 13; C++ cells 22..28 x 13..19)', [
        (3, 4, 24, 16, 28, 'mid_x_column', [27], 28),
        (7, 3, 28, 15, 1, 'ice_free_column', [], None),
    ]),
    ('block_north', 13, 21, 'the north block (i0 = 13, j0 = 21; C++ cells 13..19 x 21..27)', [
        (4, 4, 16, 24, 28, 'mid_y_column', [27], 28),
    ]),
    ('block_west', 1, 13, 'the west block (i0 = 1, j0 = 13; C++ cells 1..7 x 13..19)', [
        (4, 4, 4, 16, 20, 'mid_minus_x_column', [19], 20),
        (1, 4, 1, 16, 1, 'ice_free_column', [], None),
    ]),
    ('block_southwest', 5, 5, 'the south-west block (i0 = 5, j0 = 5; C++ cells 5..11 x 5..11)', [
        (4, 4, 8, 8, 25, 'mid_diag_column', [24], 25),
        (3, 1, 7, 5, 5, 'margin_column_8_6', [4], 5),
    ]),
    ('block_margin', 9, 1, 'the south-margin block (i0 = 9, j0 = 1; C++ cells 9..15 x 1..7)', [
        (3, 3, 11, 3, 8, 'margin_column_12_4', [7], 8),
    ]),
]

# storage levels asserted per column (1-based, esm), by phase of the profile:
# deep (age = dt), mid, just-below-surface (rel 1e-2 band from BLOCKS), top
LEVELS = {
    'dome_column': [2, 11, 21, 31, 39, 41, 43, 44, 45, 46],
    'mid_x_column': [2, 9, 16, 23, 26, 27, 28],
    'mid_y_column': [2, 9, 16, 23, 26, 27, 28],
    'mid_diag_column': [2, 8, 14, 20, 23, 24, 25],
    'mid_minus_x_column': [2, 7, 12, 16, 18, 19, 20],
    'margin_column_12_4': [2, 4, 6, 7, 8],
    'margin_column_8_6': [2, 3, 4, 5],
    'ice_free_column': [2, 20, 45],
}

# Per-block hadv pins: (block id, block i, block j, fine level k, label) cells
# whose hadv at the end of the step is pinned to the model's own measured
# value (see the module docstring). The expected literals come from
# --hadv-json; without it, no pins are generated.
HADV_CELLS = [
    ('block_east', 3, 4, 25, 'mid_x_column_level_25'),
    ('block_east', 3, 4, 26, 'mid_x_column_level_26'),
    ('block_north', 4, 4, 25, 'mid_y_column_level_25'),
    ('block_north', 4, 4, 26, 'mid_y_column_level_26'),
    ('block_west', 4, 4, 17, 'mid_minus_x_column_level_17'),
    ('block_southwest', 4, 4, 22, 'mid_diag_column_level_22'),
    ('block_southwest', 3, 1, 3, 'margin_column_8_6_level_3'),
    ('block_margin', 3, 3, 6, 'margin_column_12_4_level_6'),
    ('block_margin', 3, 3, 7, 'margin_column_12_4_level_7'),
]


def read_rows():
    rows = {}
    with open(COLS) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split(',')
            i, j, k = int(p[0]), int(p[1]), int(p[6])
            rows[(i, j, k)] = (float(p[5]), float(p[8]), float(p[9]), float(p[10]), float(p[12]))
    return rows


def fmt_num(v):
    return repr(float(v))


def fmt_field(values):
    """Nested [gi][gj][lev_p] array (padded columns), row-major, esm 1-based."""
    out = ['         "default": [']
    n_i = len(values)
    for n, irow in enumerate(values):
        out.append('          [')
        n_j = len(irow)
        for m, jrow in enumerate(irow):
            sep = ',' if m < n_j - 1 else ''
            out.append('           [' + ', '.join(fmt_num(v) for v in jrow) + ']' + sep)
        out.append('          ]' + (',' if n < n_i - 1 else ''))
    out.append('         ]')
    return '\n'.join(out)


def fmt_ks(values):
    out = ['         "default": [']
    n_i = len(values)
    for n, irow in enumerate(values):
        sep = ',' if n < n_i - 1 else ''
        out.append('          [' + ', '.join(fmt_num(v) for v in irow) + ']' + sep)
    out.append('         ]')
    return '\n'.join(out)


def replace_default(text, var_name, new_block):
    m = re.search(r'\n {8}"' + re.escape(var_name) + r'": \{\n', text)
    if not m:
        raise SystemExit(f'variable {var_name} not found')
    obj_start = m.end()
    nxt = re.search(r'\n {8}"', text[obj_start:])
    obj_end = obj_start + (nxt.start() if nxt else len(text) - obj_start)
    body = text[obj_start:obj_end]
    dm = re.search(r'"default": \[', body)
    if not dm:
        raise SystemExit(f'no default array in {var_name}')
    start = obj_start + dm.start()
    depth = 0
    i = start
    while i < len(text):
        c = text[i]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                break
        i += 1
    if depth != 0:
        raise SystemExit(f'unbalanced brackets in {var_name}')
    end = i + 1
    return text[:start] + new_block + text[end:]


def replace_tests(text, new_block):
    """Replace the model's whole "tests" array (empty or not) with new_block,
    which must start with '[' and end with ']'."""
    m = re.search(r'\n {6}"tests": \[', text)
    if not m:
        raise SystemExit('tests block not found')
    start = m.end() - 1  # the '['
    depth = 0
    i = start
    while i < len(text):
        c = text[i]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                break
        i += 1
    if depth != 0:
        raise SystemExit('unbalanced brackets in tests block')
    return text[:start] + new_block + text[i + 1:]


def assertion(variable, time, coords, expected, tol):
    return {
        "variable": variable, "time": time, "coords": coords,
        "expected": expected, "tolerance": tol,
    }


def build_tests(rows, hadv_measured, slices):
    tests = []
    for (bid, i0, j0, bdesc, columns) in BLOCKS:
        desc_lines = [
            f"One age step from age = 0 on {bdesc}: "
        ]
        col_desc = []
        for (bi, bj, di, dj, ks, label, _wide, _top) in columns:
            col_desc.append(f"({label}: C++ cell ({di},{dj}) at block cell ({bi},{bj}), ks = {ks})")
        desc_lines.append('; '.join(col_desc) + '. ')
        desc_lines.append(
            "age_storage is asserted against the dumped age_out at the storage levels listed; "
            "the whole block is integrated, so the live horizontal exchange between the block's "
            "columns is active (block-edge ghosts frozen at the initial age, exactly zero). "
            "The C++ takes one backward-Euler step with the horizontal term frozen at the old "
            "time (zero here); the harness integrates the ODE continuously, so the asserted "
            "agreement carries the exact-ODE-vs-BE vertical scheme difference of the single-column "
            "component PLUS the integrated horizontal contribution of the developing age field."
        )
        assertions = []
        for (bi, bj, di, dj, ks, label, wide_band, top_level) in columns:
            age_out = {k: rows[(di, dj, k)][4] for k in range(61)}
            for k in LEVELS[label]:
                exp = age_out[k - 1]
                if top_level is not None and k == top_level:
                    tol = {"abs": 1e-6}
                elif k in wide_band:
                    tol = {"rel": 1e-2}
                elif top_level is None:
                    tol = {"abs": 1e-6}
                else:
                    tol = {"rel": 1e-3}
                assertions.append(assertion(
                    "age_storage", DT, {"i": bi, "j": bj, "lev_s": k}, exp, tol))
        test = {
            "id": bid,
            "description": ' '.join(desc_lines),
            "time_span": {"start": 0.0, "end": DT},
            "parameter_overrides": slices[bid],
            "assertions": assertions,
        }
        # hadv pins of this block: the model's own measured horizontal term
        if hadv_measured is not None:
            for (hbid, hi, hj, k, label) in HADV_CELLS:
                if hbid != bid:
                    continue
                key = f"{bid}:{hi},{hj},{k}"
                if key not in hadv_measured:
                    raise SystemExit(f'missing measured hadv for {key} in --hadv-json')
                test["assertions"].append(assertion(
                    "hadv", DT, {"i": hi, "j": hj, "lev": k},
                    hadv_measured[key], {"rel": 1e-2}))
        tests.append(test)
    return tests


def block_slice(full, i0, j0):
    """7x7 slice of a full 31x31 field (3-D [i][j][k] or 2-D [i][j]) at block
    origin (i0, j0): block cell (bi,bj) (1-based) = C++ cell (i0+bi-1, j0+bj-1)
    = array [i0+bi-1][j0+bj-1]."""
    return [[full[i0 + bi - 1][j0 + bj - 1] for bj in range(1, 8)]
            for bi in range(1, 8)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hadv-json', default=None,
                    help='JSON file {key: value} of measured hadv at the end of the step, '
                         'key = "<block_id>:<i>,<j>,<k>"; without it, no hadv pins are generated')
    args = ap.parse_args()
    rows = read_rows()
    N = 31
    NP = 63  # padded column length: entries 1..61 = storage levels, 0 and 62 = the bed/top pad
    u3 = [[[0.0] * NP for _ in range(N)] for _ in range(N)]
    v3 = [[[0.0] * NP for _ in range(N)] for _ in range(N)]
    w3 = [[[0.0] * NP for _ in range(N)] for _ in range(N)]
    ks = [[1.0] * N for _ in range(N)]
    for (i, j, k), (H, u, v, w, _aout) in rows.items():
        u3[i][j][k + 1] = u  # entries 1..61 = storage levels 1..61
        v3[i][j][k + 1] = v
        w3[i][j][k + 1] = w
        if k == 0:
            u3[i][j][0] = u  # entry 0 = storage level 1 (the bed row)
            v3[i][j][0] = v
            w3[i][j][0] = w
            m_ks = min(int(math.floor(H / DZ_FINE)), 61)
            ks[i][j] = float(m_ks + 1)
    for i in range(N):
        for j in range(N):
            u3[i][j][NP - 1] = u3[i][j][NP - 2]  # entry 62 = storage level 61 (the top)
            v3[i][j][NP - 1] = v3[i][j][NP - 2]
            w3[i][j][NP - 1] = w3[i][j][NP - 2]

    hadv_measured = None
    if args.hadv_json:
        hadv_measured = json.load(open(args.hadv_json))

    # Defaults: the dome-block slice (C++ cells 13..19 x 13..19); every test
    # overrides its own block's slice through parameter_overrides.
    dome = next(b for b in BLOCKS if b[0] == 'block_dome')
    di0, dj0 = dome[1], dome[2]
    slices = {
        bid: {
            'u3_b': block_slice(u3, i0, j0),
            'v3_b': block_slice(v3, i0, j0),
            'w3_b': block_slice(w3, i0, j0),
            'ks_b': block_slice(ks, i0, j0),
        }
        for (bid, i0, j0, _bdesc, _columns) in BLOCKS
    }

    text = open(ESM).read()
    text = replace_default(text, 'u3_b', fmt_field(slices['block_dome']['u3_b']))
    text = replace_default(text, 'v3_b', fmt_field(slices['block_dome']['v3_b']))
    text = replace_default(text, 'w3_b', fmt_field(slices['block_dome']['w3_b']))
    text = replace_default(text, 'ks_b', fmt_ks(slices['block_dome']['ks_b']))

    tests = build_tests(rows, hadv_measured, slices)
    block = json.dumps(tests, indent=2).replace('\n', '\n      ')
    text = replace_tests(text, block)

    doc = json.loads(text)  # verify
    json.dump(doc, open('/dev/null', 'w'))  # ensure serializable
    open(ESM, 'w').write(text)
    n_hadv = 0 if hadv_measured is None else len(HADV_CELLS)
    print(f'filled u3_b, v3_b, w3_b (7x7x{NP} padded), ks_b (7x7) dome-block defaults, '
          f'{len(tests)} tests ({n_hadv} hadv pins); JSON OK')
    print('dome-block ice ks range:', min(min(r) for r in ks), '-', max(max(r) for r in ks))


if __name__ == '__main__':
    main()
