"""Retain one scalar case and conjugate a genuinely nonreal translation-phase entry."""
import json
import sys
from pathlib import Path

import numpy as np
from prepare_wyckoff_characters import complex_values, cfmt

lines = Path(sys.argv[1]).read_text().splitlines()
cursor = 1
for _ in range(int(lines[0])):
    start = cursor
    sg, site, spin, grey, pg, n, nh, ni, orbit = lines[cursor].split()
    spin, grey, n, nh, ni = map(int, (spin, grey, n, nh, ni))
    cursor += 2+n+nh
    chars = []
    names = []
    for _ in range(ni):
        names.append(lines[cursor])
        chars.append(complex_values(lines[cursor+1]))
        cursor += 2
    nq = int(lines[cursor])
    cursor += 1
    candidate = None
    for query in range(nq):
        point = lines[cursor].split()
        count = int(point[-1])
        cursor += 1
        for _ in range(count):
            values = complex_values(lines[cursor])
            if spin == 0 and grey == 0 and query > 0:
                for j, value in enumerate(values):
                    if np.max(np.abs(chars[j].imag)) < 1e-8 and abs(value.imag) > 0.1:
                        candidate = cursor, j, value, query+1
                        break
            cursor += 1
    if candidate is None:
        continue
    at, column, original, query = candidate
    values = complex_values(lines[at])
    values[column] = values[column].conjugate()
    lines[at] = lines[at].split()[0]+" "+cfmt(values)
    Path(sys.argv[2]).write_text("1\n"+"\n".join(lines[start:cursor])+"\n")
    print(json.dumps(dict(sg=int(sg), site=site, spin=spin, grey=grey, query=query,
                          label=names[column], original=[original.real, original.imag],
                          changed=[values[column].real, values[column].imag])))
    break
else:
    raise SystemExit("No translation-phase-only negative case found")
