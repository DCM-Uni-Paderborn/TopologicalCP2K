"""Compare grouped all-electron GAPW property frames with a full diagonalization."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from run_property_symmetry_controls import digest, native_libraries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root, destination = args.root.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    binary = root / "build-mpi/bin/cp2k.psmp"
    libraries, binary_hash = native_libraries(root, "build-mpi"), digest(binary)
    source = root / "tests/QS/regtest-property-wilson/helium-gapw-common.inc"
    nnkp = root / "tests/QS/regtest-topology/helium.nnkp"
    comparator = Path(__file__).with_name("compare_property_wilson.py")
    inputs = [source, nnkp, root / "data/BASIS_MOLOPT_UZH", root / "data/POTENTIAL_UZH"]
    source_hashes = {str(p): digest(p) for p in inputs}
    environment = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OPENBLAS_NUM_THREADS="1",
                       OMP_NUM_THREADS="2", OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(root / "build-serial/openblas-thread-safe/lib"))
    cases, comparisons = [], []
    runs = [("full", 2, 0)] + [(backend, ranks, size) for backend in ("K290", "SPGLIB")
                             for ranks, size in ((2, 1), (4, 2))]
    for backend, ranks, size in runs:
        work = destination / f"{backend.lower()}-r{ranks}-g{size}"
        work.mkdir()
        text = source.read_text().replace("${CASE}", "states")
        text = text.replace("${PROPERTY_SYMMETRY}", "F" if backend == "full" else "T")
        text = text.replace("${PROPERTY_BACKEND}", "K290" if backend == "full" else backend)
        text = text.replace("${PROPERTY_GROUP_SIZE}", str(size))
        text = text.replace("../regtest-topology/helium.nnkp", nnkp.name)
        if "${" in text:
            raise ValueError("Unresolved variable")
        (work / "input.inp").write_text(text)
        (work / nnkp.name).write_bytes(nnkp.read_bytes())
        command = ["mpiexec", "-n", str(ranks), str(binary), "-i", "input.inp", "-o", "output.out"]
        with (work / "launcher.log").open("w") as stream:
            result = subprocess.run(command, cwd=work, env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=180, check=False)
        output = (work / "output.out").read_text()
        accepted = result.returncode == 0 and "PROGRAM ENDED" in output and "[ABORT]" not in output
        groups = re.findall(r"Property scalar validation MPI groups:\s*(\d+)\s+ranks per group:\s*(\d+)", output)
        if backend != "full":
            accepted &= groups == [(str(ranks // size), str(size))]
            accepted &= "Scalar H/S snapshots unavailable" not in output
        record = dict(name=work.name, command=command, returncode=result.returncode, accepted=accepted,
                      validation_groups=groups, files={p.name: digest(p) for p in work.iterdir()})
        cases.append(record)
        (work / "run.json").write_text(json.dumps(record, indent=2) + "\n")
        print(work.name, accepted, flush=True)
        if not accepted:
            raise RuntimeError(work.name)
        if backend != "full":
            with (work / "comparison.log").open("w") as stream:
                subprocess.run([sys.executable, str(comparator), str(root),
                                str(destination / "full-r2-g0/states"), str(work / "states"),
                                str(work / "comparison.json")], stdout=stream,
                               stderr=subprocess.STDOUT, env=environment, timeout=180, check=True)
            comparisons.append(json.loads((work / "comparison.json").read_text()))
    if libraries != native_libraries(root, "build-mpi") or binary_hash != digest(binary):
        raise ValueError("Native implementation changed during validation")
    if any(digest(Path(p)) != value for p, value in source_hashes.items()):
        raise ValueError("Input source changed during validation")
    accepted = all(c["accepted"] for c in comparisons) and len(comparisons) == 4
    report = dict(cases=cases, comparisons=comparisons, accepted=accepted, sources=source_hashes,
                  binary_sha256=binary_hash, library_sha256=libraries,
                  methods={str(Path(__file__)): digest(Path(__file__)), str(comparator): digest(comparator)})
    (destination / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
