"""Check collective rejection of truncated degenerate SOC windows in MPI property groups."""

import argparse
import json
import os
from pathlib import Path
import subprocess

from run_property_symmetry_controls import digest, native_libraries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root, destination = args.root.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    binary = root / "build-mpi/bin/cp2k.psmp"
    libraries = native_libraries(root, "build-mpi")
    binary_hash = digest(binary)
    source = root / "tests/QS/regtest-property-wilson/neon-common.inc"
    text = source.read_text().replace("${CASE}", "window").replace("${PROPERTY_SYMMETRY}", "T")
    text = text.replace("ADDED_MOS 8", "ADDED_MOS 0")
    text = text.replace("EXCLUDE_BANDS 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26", "EXCLUDE_BANDS 9 10")
    environment = dict(os.environ, CP2K_DATA_DIR=str(root / "data"),
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2", OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(root / "build-serial/openblas-thread-safe/lib"))
    records = []
    for backend in ["K290", "SPGLIB"]:
        for group_size in [1, 2]:
            work = destination / f"{backend.lower()}-g{group_size}"
            work.mkdir()
            data = text.replace("${PROPERTY_BACKEND}", backend)
            data = data.replace("STATE_EXPORT T", f"STATE_EXPORT T\n        SYMMETRY_PARALLEL_GROUP_SIZE {group_size}")
            (work / "input.inp").write_text(data)
            command = ["mpiexec", "-n", "4", str(binary), "-i", "input.inp", "-o", "output.out"]
            with (work / "launcher.log").open("w") as stream:
                result = subprocess.run(command, cwd=work, env=environment, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=120, check=False)
            output = (work / "output.out").read_text()
            accepted = (result.returncode != 0 and "PROGRAM ENDED" not in output
                        and "Property SOC scalar window cuts a degeneracy" in output)
            records.append(dict(name=work.name, command=command, returncode=result.returncode,
                                accepted=accepted, files={p.name: digest(p) for p in work.iterdir()}))
            print(work.name, accepted, flush=True)
    if libraries != native_libraries(root, "build-mpi") or binary_hash != digest(binary):
        raise ValueError("Native code changed during validation")
    report = dict(cases=records, accepted=all(r["accepted"] for r in records),
                  binary_sha256=binary_hash, library_sha256=libraries, source_sha256=digest(source))
    (destination / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    if not report["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
