"""Compare bounded and collective links for a 56-orbital all-electron GAPW cell."""

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
    common = root / "tests/QS/regtest-property-wilson/helium-gapw-common.inc"
    nnkp = root / "tests/QS/regtest-property-wilson/helium-links.nnkp"
    comparator = Path(__file__).with_name("compare_property_wilson.py")
    source_hashes = {str(p): digest(p) for p in [common, nnkp, root / "data/BASIS_MOLOPT_UZH",
                                                root / "data/POTENTIAL_UZH"]}
    environment = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OPENBLAS_NUM_THREADS="1",
                       OMP_NUM_THREADS="2", OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(root / "build-serial/openblas-thread-safe/lib"))
    runs = [("full", 1024)] + [(b, m) for b in ["K290", "SPGLIB"] for m in [1024, 1]]
    cases, comparisons = [], []
    for backend, memory in runs:
        work = destination / f"{backend.lower()}-m{memory}"
        work.mkdir()
        text = common.read_text().replace("${CASE}", "states")
        text = text.replace("${PROPERTY_SYMMETRY}", "F" if backend == "full" else "T")
        text = text.replace("${PROPERTY_BACKEND}", "K290" if backend == "full" else backend)
        text = text.replace("${PROPERTY_GROUP_SIZE}", "1")
        text = text.replace("${PROPERTY_NNKP}", nnkp.name)
        text = text.replace("${PROPERTY_WILSON}", "T")
        text = text.replace("STATE_EXPORT T", f"STATE_EXPORT T\n        SYMMETRY_MAX_MEMORY_MB {memory}")
        assert text.count("He 0.5 0.5 0.5") == 1
        text = text.replace("He 0.5 0.5 0.5", "He 0.25 0.25 0.25\n      He 0.25 0.75 0.75\n"
                            "      He 0.75 0.25 0.75\n      He 0.75 0.75 0.25")
        assert "${" not in text
        (work / "input.inp").write_text(text)
        (work / nnkp.name).write_bytes(nnkp.read_bytes())
        command = ["mpiexec", "-n", "4", str(binary), "-i", "input.inp", "-o", "output.out"]
        with (work / "launcher.log").open("w") as stream:
            result = subprocess.run(command, cwd=work, env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=300, check=False)
        output = (work / "output.out").read_text()
        accepted = result.returncode == 0 and "PROGRAM ENDED" in output and "[ABORT]" not in output
        fallback = "Directed-link snapshots unavailable; using collective link construction." in output
        scalar_groups = re.findall(r"Property scalar validation MPI groups:\s*(\d+)\s+ranks per group:\s*(\d+)", output)
        link_groups = re.findall(r"Property directed-link MPI groups:\s*(\d+)\s+ranks per group:\s*(\d+)", output)
        if backend != "full":
            accepted &= scalar_groups == [("4", "1")] and "Scalar H/S snapshots unavailable" not in output
            accepted &= fallback == (memory == 1)
            accepted &= link_groups == ([] if memory == 1 else scalar_groups)
        record = dict(name=work.name, command=command, returncode=result.returncode, accepted=accepted,
                      fallback=fallback, scalar_groups=scalar_groups, link_groups=link_groups,
                      files={p.name: digest(p) for p in work.iterdir()})
        cases.append(record)
        (work / "run.json").write_text(json.dumps(record, indent=2) + "\n")
        print(work.name, accepted, fallback, flush=True)
        if not accepted:
            raise RuntimeError(work.name)
        if backend != "full":
            with (work / "comparison.log").open("w") as stream:
                subprocess.run([sys.executable, str(comparator), str(root),
                                str(destination / "full-m1024/states"), str(work / "states"),
                                str(work / "comparison.json")], stdout=stream, stderr=subprocess.STDOUT,
                               env=dict(environment, OMP_NUM_THREADS="1"), timeout=300, check=True)
            comparisons.append(json.loads((work / "comparison.json").read_text()))
    if libraries != native_libraries(root, "build-mpi") or binary_hash != digest(binary):
        raise ValueError("Native implementation changed during validation")
    if any(digest(Path(p)) != value for p, value in source_hashes.items()):
        raise ValueError("Input source changed during validation")
    accepted = len(comparisons) == 4 and all(c["accepted"] for c in comparisons)
    report = dict(accepted=accepted, cases=cases, comparisons=comparisons, sources=source_hashes,
                  binary_sha256=binary_hash, library_sha256=libraries,
                  methods={str(Path(__file__)): digest(Path(__file__)), str(comparator): digest(comparator)})
    (destination / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
