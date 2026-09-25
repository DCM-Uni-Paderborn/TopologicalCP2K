"""Validate independent property MPI groups against full physical eigenframes and links."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

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
    environment = dict(os.environ, CP2K_DATA_DIR=str(root / "data"),
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2", OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(root / "build-serial/openblas-thread-safe/lib"))
    records, comparisons = [], []
    comparator = Path(__file__).with_name("compare_property_wilson.py")
    sources = {}

    def run(material, backend, ranks, group_size, diagnostic=None):
        name = f"{material}-{backend.lower()}-r{ranks}-g{group_size}"
        work = destination / name
        work.mkdir()
        source = root / f"tests/QS/regtest-property-wilson/{material}-common.inc"
        sources[str(source)] = digest(source)
        text = source.read_text().replace("${CASE}", "states")
        text = text.replace("${PROPERTY_SYMMETRY}", "F" if backend == "full" else "T")
        text = text.replace("${PROPERTY_BACKEND}", "K290" if backend == "full" else backend)
        if "${PROPERTY_GROUP_SIZE}" in text:
            text = text.replace("${PROPERTY_GROUP_SIZE}", str(group_size))
        else:
            text = text.replace("STATE_EXPORT T", f"STATE_EXPORT T\n        SYMMETRY_PARALLEL_GROUP_SIZE {group_size}")
        if material == "helium":
            nnkp = root / "tests/QS/regtest-topology/helium.nnkp"
            sources[str(nnkp)] = digest(nnkp)
            text = text.replace("../regtest-topology/helium.nnkp", str(nnkp))
        if "${" in text:
            raise ValueError("Unresolved input variable")
        (work / "input.inp").write_text(text)
        command = ["mpiexec", "-n", str(ranks), str(binary), "-i", "input.inp", "-o", "output.out"]
        start = time.monotonic()
        with (work / "launcher.log").open("w") as stream:
            result = subprocess.run(command, cwd=work, env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=180, check=False)
        output = (work / "output.out").read_text()
        if diagnostic:
            accepted = result.returncode != 0 and diagnostic in output and "PROGRAM ENDED" not in output
        else:
            accepted = result.returncode == 0 and "PROGRAM ENDED" in output and "[ABORT]" not in output
        groups = re.findall(r"Property scalar MPI groups:\s*(\d+)\s+ranks per group:\s*(\d+)", output)
        diagonalizations = re.findall(r"Property scalar diagonalizations:\s*(\d+)\s*/\s*(\d+)", output)
        if not diagnostic and backend != "full":
            accepted &= bool(groups) and len(groups) == len(diagonalizations)
            for (count, size), (representatives, _) in zip(groups, diagonalizations):
                count, size, representatives = int(count), int(size), int(representatives)
                expected_size = ranks if group_size == 0 else group_size
                if group_size == -1:
                    expected_size = min(d for d in range(1, ranks + 1)
                                        if ranks % d == 0 and ranks // d <= representatives)
                accepted &= size == expected_size and count * size == ranks and count <= representatives
        record = dict(name=name, command=command, returncode=result.returncode, accepted=accepted,
                      elapsed=time.monotonic() - start, expected_diagnostic=diagnostic,
                      groups=groups, diagonalizations=diagonalizations,
                      files={p.name: digest(p) for p in work.iterdir() if p.is_file()})
        (work / "run.json").write_text(json.dumps(record, indent=2) + "\n")
        records.append(record)
        print(name, accepted, f"{record['elapsed']:.3f}s", groups, flush=True)
        if not accepted:
            raise RuntimeError(f"Failed calculation: {name}")
        return work / "states"

    for material in ["helium", "neon", "stanene"]:
        reference = run(material, "full", 2, 0)
        for backend in ["K290", "SPGLIB"]:
            for ranks, group_size in [(2, 0), (2, -1), (4, -1), (4, 2)]:
                candidate = run(material, backend, ranks, group_size)
                comparison = candidate.parent / "comparison.json"
                with (candidate.parent / "comparison.log").open("w") as stream:
                    subprocess.run([sys.executable, str(comparator), str(root), str(reference),
                                    str(candidate), str(comparison)], check=True, timeout=180,
                                   env=dict(environment, OMP_NUM_THREADS="1"), stdout=stream,
                                   stderr=subprocess.STDOUT)
                data = json.loads(comparison.read_text())
                comparisons.append(data)
                print("comparison", candidate.parent.name, data["accepted"],
                      data["subspace_residual"], data["link_covariance_error"], flush=True)
    run("helium", "K290", 2, -2, "SYMMETRY_PARALLEL_GROUP_SIZE must be at least -1")
    run("neon", "K290", 4, 3, "not divisible by the kpoint group size")
    run("stanene", "K290", 4, 1, "Too many kpoint groups")
    run("stanene", "SPGLIB", 4, 1, "Too many kpoint groups")
    if libraries != native_libraries(root, "build-mpi") or binary_hash != digest(binary):
        raise ValueError("Native code changed during validation")
    if any(digest(Path(path)) != value for path, value in sources.items()):
        raise ValueError("Input source changed during validation")
    report = dict(cases=records, comparisons=comparisons, accepted=True, sources=sources,
                  binary_sha256=binary_hash, library_sha256=libraries,
                  methods={str(Path(__file__)): digest(Path(__file__)), str(comparator): digest(comparator)},
                  environment={key: environment[key] for key in ["CP2K_DATA_DIR", "OPENBLAS_NUM_THREADS",
                                                                "OMP_NUM_THREADS", "OMP_STACKSIZE",
                                                                "DYLD_LIBRARY_PATH"]})
    (destination / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print("PASS", len(records), "native cases;", len(comparisons), "physical comparisons", flush=True)


if __name__ == "__main__":
    main()
