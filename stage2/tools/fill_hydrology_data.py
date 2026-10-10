#!/usr/bin/env python3
"""Fill the data fields of stage2/hydrology.esm from the Stage-1 hydrology dumps.

This script performs ONLY numeric data extraction (AGENTS.md: "Scripts MAY be used
to extract numeric test tuples from C++ dumps into an existing hand-authored tests
block" -- here, field data into a hand-authored component). It replaces the
"default" arrays of the H, surface_input, and basal_melt parameters with the
dumped fields (columns.csv columns H, surface_input_m_s, basal_melt_m_s), and the
W_law / G_law parameters with the full pointwise flux-law sample table
(flux_law.csv columns W, G). The equations, structure, and factoring are
hand-authored and are not touched.

The esm grid is 1-based (esm cell (i,j) = C++ cell (i-1,j-1)); the nested default
arrays are row-major over (i, j) with row r = C++ i = esm i-1.
"""
import json
import re
import sys

ESM = 'stage2/hydrology.esm'
COLS = 'stage1/dumps/hydrology/columns.csv'
LAW = 'stage1/dumps/hydrology/flux_law.csv'

def read_csv(path, columns):
    with open(path) as f:
        lines = [l for l in f if not l.startswith('#')]
    return [dict(zip(columns, l.split(','))) for l in lines]

def fmt_num(x):
    return repr(float(x))

def fmt_field(values31):
    """Format a 31x31 nested array in the file's existing style."""
    rows = []
    for row in values31:
        inner = ',\n'.join(' ' * 8 + fmt_num(v) for v in row)
        rows.append(' ' * 6 + '[\n' + inner + '\n' + ' ' * 6 + ']')
    return '     "default": [\n' + ',\n'.join(rows) + '\n     ]'

def fmt_list(values):
    return '     "default": [\n' + ',\n'.join(' ' * 6 + fmt_num(v) for v in values) + '\n     ]'

def replace_default(text, var_name, new_block):
    """Replace the "default" array of variable var_name (at the model's 5-space
    variable indent) with new_block (a '     "default": [...]' fragment)."""
    # locate the variable object
    m = re.search(r'\n    "' + re.escape(var_name) + r'": \{\n', text)
    if not m:
        raise SystemExit(f'variable {var_name} not found')
    obj_start = m.end()
    # find the "default" key inside this object (before the next 4-space key)
    nxt = re.search(r'\n    "', text[obj_start:])
    obj_end = obj_start + (nxt.start() if nxt else len(text) - obj_start)
    body = text[obj_start:obj_end]
    dm = re.search(r'"default": \[', body)
    if not dm:
        raise SystemExit(f'no default array in {var_name}')
    # bracket-match from the start of the default array
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
    end = i + 1  # include the closing ]
    return text[:start] + new_block + text[end:]

def main():
    text = open(ESM).read()

    cols = read_csv(COLS, ['i', 'j', 'x', 'y', 'r', 'H', 'mask', 'surface_input_rate',
                           'basal_melt_rate', 'sliding_speed', 'surface_input_m_s',
                           'basal_melt_m_s', 'W_till_in', 'W_till_out', 'W_in', 'W_out',
                           'overburden_pressure', 'flux_u', 'flux_v'])
    N = 31
    H = [[0.0] * N for _ in range(N)]
    surf = [[0.0] * N for _ in range(N)]
    basal = [[0.0] * N for _ in range(N)]
    for r in cols:
        i, j = int(r['i']), int(r['j'])
        H[i][j] = float(r['H'])
        surf[i][j] = float(r['surface_input_m_s'])
        basal[i][j] = float(r['basal_melt_m_s'])

    text = replace_default(text, 'H', fmt_field(H))
    text = replace_default(text, 'surface_input', fmt_field(surf))
    text = replace_default(text, 'basal_melt', fmt_field(basal))

    law = read_csv(LAW, ['W', 'G', 'alpha', 'beta', 'k', 'K', 'q_mag'])
    Wl = [float(r['W']) for r in law]
    Gl = [float(r['G']) for r in law]
    assert len(Wl) == 99 and len(Gl) == 99
    text = replace_default(text, 'W_law', fmt_list(Wl))
    text = replace_default(text, 'G_law', fmt_list(Gl))

    # law index set: 42 -> 99
    text = text.replace('"law": {\n   "kind": "interval",\n   "size": 42',
                        '"law": {\n   "kind": "interval",\n   "size": 99')

    json.loads(text)  # verify
    open(ESM, 'w').write(text)
    print('filled H, surface_input, basal_melt (31x31) and W_law, G_law (99 samples); JSON OK')

if __name__ == '__main__':
    main()
