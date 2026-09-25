"""Run independent Stanene TRIM symmetry controls without changing a stored reference input."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def native_libraries(root, build):
    paths = {p.resolve() for p in (root / build / "src").glob("libcp2k*")
             if p.is_file() and (".so" in p.name or p.suffix == ".dylib")}
    return {str(p.relative_to(root)): digest(p) for p in sorted(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpi", action="store_true")
    parser.add_argument("--symmetric-scf", action="store_true")
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--scf-grid", type=int, default=4)
    parser.add_argument("--rel-cutoff", type=int)
    parser.add_argument("--wilson", action="store_true")
    parser.add_argument("--wilson-mesh", type=int, nargs=2, default=[8, 5])
    parser.add_argument("--refinements", type=int, default=1)
    parser.add_argument("--no-state-export", action="store_true")
    args = parser.parse_args()
    args.root = args.root.resolve()
    args.destination = args.destination.resolve()
    args.destination.mkdir(parents=True, exist_ok=False)
    version = "psmp" if args.mpi else "ssmp"
    build = "build-mpi" if args.mpi else "build-serial"
    binary = args.root / build / "bin" / f"cp2k.{version}"
    reference = args.root / "tests/QS/regtest-topology/stanene-tqc.inp"
    text = reference.read_text()
    if args.symmetric_scf:
        text = text.replace("FULL_GRID T", "FULL_GRID F").replace("SYMMETRY F", "SYMMETRY T")
        text = text.replace("SCHEME MONKHORST-PACK 4 4 1", "SCHEME MONKHORST-PACK 4 4 1\n      GAMMA_CENTERED T")
        text = text.replace("EPS_SCF 1.0E-8", "EPS_SCF 1.0E-10")
    if args.cutoff:
        text = text.replace("CUTOFF 200", f"CUTOFF {args.cutoff}")
    text = text.replace("SCHEME MONKHORST-PACK 4 4 1", f"SCHEME MONKHORST-PACK {args.scf_grid} {args.scf_grid} 1")
    if args.rel_cutoff:
        text = text.replace("REL_CUTOFF 40", f"REL_CUTOFF {args.rel_cutoff}")
    if args.wilson:
        text = text.replace("INVERSION_TQC T", f"Z2 T\n        WILSON_MESH {args.wilson_mesh[0]} {args.wilson_mesh[1]}\n        WILSON_MAX_REFINEMENT {args.refinements}")
        text = text.replace("KPOINTS_SOURCE TRIM", "KPOINTS_SOURCE WILSON").replace("TQC_DIMENSION 2", "REQUIRE_GLOBAL_GAP T")
    environment = dict(os.environ, CP2K_DATA_DIR=str(args.root / "data"),
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2", OMP_STACKSIZE="64M",
                       DYLD_LIBRARY_PATH=str(args.root / "build-serial/openblas-thread-safe/lib"))
    for backend in ["full", "K290", "SPGLIB"]:
        work = args.destination / backend.lower()
        work.mkdir()
        extra = "STATE_EXPORT " + ("F" if args.no_state_export else "T") + "\n        SYMMETRY " + ("F" if backend == "full" else "T")
        if backend != "full":
            extra += "\n        SYMMETRY_BACKEND " + backend
        data = text.replace("SEED_NAME stanene", "SEED_NAME stanene\n        " + extra)
        if data == text:
            raise ValueError("Expected insertion point missing")
        (work / "input.inp").write_text(data)
        command = [str(binary), "-i", "input.inp", "-o", "output.out"]
        if args.mpi:
            command = ["mpiexec", "-n", str(args.ranks)] + command
        start = time.monotonic()
        libraries = native_libraries(args.root, build)
        with (work / "launcher.log").open("w") as stream:
            result = subprocess.run(command, cwd=work, env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=600, check=False)
        record = dict(command=command, returncode=result.returncode, elapsed=time.monotonic() - start,
                      binary_sha256=digest(binary), reference_sha256=digest(reference),
                      library_sha256=libraries,
                      environment={key: environment[key] for key in ["CP2K_DATA_DIR", "OPENBLAS_NUM_THREADS",
                                                                    "OMP_NUM_THREADS", "OMP_STACKSIZE",
                                                                    "DYLD_LIBRARY_PATH"]})
        if libraries != native_libraries(args.root, build):
            raise ValueError("Native library changed during calculation")
        record["files"] = {p.name: digest(p) for p in work.iterdir() if p.is_file()}
        (work / "run.json").write_text(json.dumps(record, indent=2) + "\n")
        print(backend, result.returncode, record["elapsed"], flush=True)


if __name__ == "__main__":
    main()
