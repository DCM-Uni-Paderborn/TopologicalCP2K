"""Require diagnostic failures for invalid property and SCF symmetry requests."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

from run_property_symmetry_controls import native_libraries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    args.destination.mkdir(parents=True, exist_ok=False)
    binary = root / "build-serial/bin/cp2k.ssmp"
    base = (root / "tests/QS/regtest-property-wilson/neon-common.inc").read_text()
    base = base.replace("${CASE}", "negative").replace("${PROPERTY_SYMMETRY}", "T")
    cases = []
    for backend in ["K290", "SPGLIB"]:
        data = base.replace("${PROPERTY_BACKEND}", backend)
        window = data.replace("ADDED_MOS 8", "ADDED_MOS 0")
        window = window.replace("EXCLUDE_BANDS 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26", "EXCLUDE_BANDS 9 10")
        cases.append((backend.lower() + "-window", window, "Property SOC scalar window cuts a degeneracy"))
        general = data.replace("SCHEME MONKHORST-PACK 2 2 2", "SCHEME GENERAL\n      SYMMETRY_BACKEND " + backend +
                               "\n      KPOINT 0 0 0 0.5\n      KPOINT 0.25 0 0 0.5")
        general = general.replace("FULL_GRID T", "FULL_GRID F").replace("SYMMETRY F", "SYMMETRY T")
        cases.append((backend.lower() + "-closure", general, "to be closed under"))
    invalid = base.replace("${PROPERTY_BACKEND}", "K290").replace("SYMMETRY_BACKEND K290", "SYMMETRY_BACKEND K290\n        EPS_SYMMETRY_OPERATORS 0")
    cases.append(("invalid-tolerance", invalid, "Invalid property symmetry tolerances or memory"))
    records = []
    environment = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2",
                       OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(root / "build-serial/openblas-thread-safe/lib"))
    libraries = native_libraries(root, "build-serial")
    for name, data, diagnostic in cases:
        work = args.destination / name
        work.mkdir()
        (work / "input.inp").write_text(data)
        with (work / "launcher.log").open("w") as stream:
            result = subprocess.run([str(binary), "-i", "input.inp", "-o", "output.out"], cwd=work,
                                    env=environment, stdout=stream, stderr=subprocess.STDOUT, timeout=120, check=False)
        output = (work / "output.out").read_text()
        passed = result.returncode != 0 and diagnostic in output and "PROGRAM ENDED" not in output
        records.append(dict(name=name, returncode=result.returncode, expected=diagnostic, passed=passed,
                            files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in work.iterdir()}))
        print(name, passed, flush=True)
    if libraries != native_libraries(root, "build-serial"):
        raise ValueError("Native library changed during negative controls")
    report = dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  library_sha256=libraries, cases=records,
                  accepted=all(r["passed"] for r in records))
    (args.destination / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    if not report["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
