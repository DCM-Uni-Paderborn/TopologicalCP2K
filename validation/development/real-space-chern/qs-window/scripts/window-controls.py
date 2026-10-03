"""Local finite-window QS controls, not part of the proposed CP2K patch."""

import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    binary = args.binary.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    base = (root / "tests/QS/regtest-kubo-hall/window-gpw.inp").read_text()
    controls = {
        "reverse-plane": (base.replace("PLANE XZ", "PLANE ZX").replace("1.3 1.8", "1.8 1.3"), None),
        "broad-window": (base.replace("1.3 1.8", "100000 100000"), None),
        "negative-width": (base.replace("1.3 1.8", "-1.3 1.8"), "WINDOW_SIGMA must be positive"),
        "missing-energy": (base.replace("      ENERGY [hartree] -0.15\n", ""), "requires explicit ENERGY"),
        "unresolved-gap": (base.replace("PLANE XZ", "PLANE XZ\n      EPS_GAP 10"), "status -5"),
        "periodic-cell": (
            base.replace("PERIODIC NONE", "PERIODIC XYZ").replace("POISSON_SOLVER MT", "POISSON_SOLVER PERIODIC"),
            "requires CELL/POISSON PERIODIC NONE",
        ),
    }
    env = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1")
    for name, (text, expected_error) in controls.items():
        directory = work / name
        directory.mkdir(exist_ok=True)
        (directory / "input.inp").write_text(text)
        run = subprocess.run([str(binary), "-i", "input.inp"], cwd=directory, env=env, capture_output=True, text=True, timeout=120)
        output = run.stdout + run.stderr
        (directory / "run.out").write_text(output)
        if expected_error is not None:
            assert run.returncode != 0 and expected_error in output, (name, run.returncode, output[-2000:])
            print(name, "rejected as expected")
        else:
            assert run.returncode == 0, (name, output[-2000:])
            weight = float(re.search(r"Window occupied weight:\s+(\S+)", output)[1])
            marker = float(re.search(r"REAL_SPACE_CHERN\| Marker:\s+(\S+)", output)[1])
            expected = 0.684673784725126 if name == "reverse-plane" else 1.0
            assert abs(weight - expected) < 2e-8 and abs(marker) < 1e-10, (name, weight, marker)
            print(name, "weight", weight, "marker", marker)


if __name__ == "__main__":
    main()
