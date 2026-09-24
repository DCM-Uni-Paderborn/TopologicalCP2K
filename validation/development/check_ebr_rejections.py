"""Check native EBR-reference rejection with dimension-preserving corruptions."""
import argparse
import copy
import json
import os
import subprocess
from pathlib import Path

import numpy as np

from prepare_ebr_characters import cfmt, fmt
from prepare_wyckoff_characters import complex_values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    lines = iter(args.fixture.read_text().splitlines())
    selected = None
    for _ in range(int(next(lines))):
        header = next(lines)
        sg, spin, n, nr, nq, orbit, site = header.split()
        n, nr, nq = map(int, (n, nr, nq))
        rows = [header, next(lines)]
        rows.extend(next(lines) for _ in range(n+nr))
        points = []
        for _ in range(nq):
            point = next(lines)
            count = int(point.split()[-1])
            start = len(rows)
            rows.append(point)
            rows.extend(next(lines) for _ in range(count))
            points.append((start, count))
        if (sg, spin, site) == ("75", "1", "a"):
            selected = rows, n, nr, points
    assert selected is not None
    rows, n, nr, points = selected
    altered = {"baseline": (copy.copy(rows), None)}
    bad = copy.copy(rows)
    index = next(start+2 for start, count in reversed(points) if count > 1)
    prefix = bad[index].split()[0]
    values = complex_values(bad[index])
    values[0] += .25
    bad[index] = prefix+" "+cfmt(values)
    altered["boundary-character"] = bad, "Independent EBR character sequence did not match"
    bad = copy.copy(rows)
    geometry = np.fromstring(bad[1], sep=" ")
    geometry[9:] += [.5, .5, 0]
    bad[1] = fmt(geometry)
    altered["site-shift"] = bad, "Independent EBR character sequence did not match"
    bad = copy.copy(rows)
    for start, count in points:
        for index in range(start+1, start+count+1):
            prefix = bad[index].split()[0]
            values = complex_values(bad[index])
            values[1] = values[0]
            bad[index] = prefix+" "+cfmt(values)
    altered["duplicate-character-column"] = bad, "EBR reference column multiplicity exceeds native count"
    bad = copy.copy(rows)
    bad[3] = bad[3].split("(")[0]+cfmt(np.zeros(4, dtype=complex))
    altered["spin-lift"] = bad, "EBR spin lift mismatch"
    bad = copy.copy(rows)
    start, count = next(p for p in points if p[1] > 1)
    bad[start+2] = "1 "+cfmt(complex_values(bad[start+2]))
    altered["duplicate-operation"] = bad, "Duplicate EBR character record"
    altered["trailing-data"] = rows+["unexpected"], "Unexpected trailing EBR fixture data"
    results = []
    for name, (changed, message) in altered.items():
        fixture = args.output/f"{name}.dat"
        fixture.write_text("1\n"+"\n".join(changed)+"\n")
        result = subprocess.run([str(args.binary.resolve()), f"--ebr-reference={fixture}"],
            env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "2", "OMP_STACKSIZE": "64M"},
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        (args.output/f"{name}.log").write_text(result.stdout)
        if message is None:
            assert result.returncode == 0 and "Site-induced atomic band tests passed." in result.stdout
        else:
            assert result.returncode != 0 and message in result.stdout, (name, result.stdout)
        results.append(dict(test=name, returncode=result.returncode, expected_diagnostic=message, passed=True))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
