"""Retain independently varied stanene controls and their complete exports."""

import argparse
import ctypes
import json
from pathlib import Path
import re
import sys

import numpy as np
import scipy

from archive_stanene_convergence import digest, pack
from run_stanene_controls import cases, inputs
from replay_stanene_controls import verify


def retain(path, files, metadata):
    if not path.exists():
        return pack(path, files, metadata)
    checksum = digest(path)
    manifest = verify(path, checksum)
    assert manifest == dict(metadata, files={name: dict(bytes=source.stat().st_size, sha256=digest(source))
                                            for name, source in files.items()})
    return dict(file=path.name, bytes=path.stat().st_size, sha256=checksum)


def runtime(root):
    result = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__)
    library = ctypes.CDLL(str(root / "build-mpi/src/libcp2k.2026.2.dylib"))
    result["native_library"] = library._name
    system = ctypes.CDLL(None)
    system._dyld_image_count.restype = ctypes.c_uint32
    system._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    system._dyld_get_image_name.restype = ctypes.c_char_p
    images = {}
    for i in range(system._dyld_image_count()):
        name = system._dyld_get_image_name(i).decode()
        if name.startswith((str(root), "/opt/homebrew", sys.prefix)) and Path(name).is_file():
            images[name] = digest(Path(name))
    assert any("openblas-thread-safe" in name for name in images)
    result["loaded_image_sha256"] = images
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--stage", action="store_true", help="Archive completed cases without finalizing the index")
    args = parser.parse_args()
    root, work, destination = args.root.resolve(), args.work.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    archives, summary = [], []
    shared = {}
    for name in ("run_stanene_controls.py", "archive_stanene_controls.py",
                 "replay_stanene_controls.py", "archive_stanene_convergence.py", "test_stanene_controls.py",
                 "summarize_stanene_controls.py"):
        shared["scripts/" + name] = Path(__file__).with_name(name)
    shared["input-integrity-test.log"] = work / "input-integrity-test.log"
    for case in cases():
        directory = work / case["label"]
        if args.stage and (not (directory / "run.json").exists()
                           or (case["mode"] == "bloch" and not (directory / "analysis.json").exists())):
            continue
        run = json.loads((directory / "run.json").read_text())
        assert run["case"] == case and run["returncode"] == 0
        assert run["scf_converged"] and run["completed"]
        for name, expected in inputs(case, root).items():
            assert (directory / name).read_text() == expected
        for name, expected in run["outputs"].items():
            assert digest(directory / name) == expected
        for name, expected in run["provenance"]["sources"].items():
            assert digest(root / name) == expected
            if not name.endswith(("cp2k.psmp", ".dylib")):
                shared[name] = root / name
        assert digest(Path(__file__).with_name("run_stanene_controls.py")) == run["provenance"]["runner_sha256"]
        names = ["input.inp", "run.json", *run["outputs"]]
        if case["mode"] == "bloch":
            names += ["mesh.nnkp", "analysis.json", "analysis.log"]
            analysis = json.loads((directory / "analysis.json").read_text())
            assert analysis["provenance"] == run["provenance"]
            for name, expected in analysis["inputs"].items():
                assert digest(directory / name) == expected
            entry = dict(case=case, run=run, analysis=analysis)
        else:
            text = (directory / "output.out").read_text()
            assert "Wilson surface sampling converged." in text
            indices = re.findall(r"Converged Z2 invariant:\s+(\d+)", text)
            assert len(indices) == 1
            gaps = re.findall(r"Sampled indirect gap \[eV\]:\s+(\S+)", text)
            entry = dict(case=case, run=run, wilson_index=int(indices[0]),
                         wilson_gap_ev=float(gaps[-1]))
        summary.append(entry)
        archives.append(retain(destination / (case["label"] + ".tar.gz"),
                             {name: directory / name for name in names}, dict(case=case)))
        print(json.dumps(archives[-1]), flush=True)
    if args.stage:
        return
    for label, previous in (("dzvp-200-scf8", "mesh12"), ("tzvp-400-scf8", "mesh12-tzvp")):
        file = root / "build-mpi/delayed-native-materials-verified" / (previous + ".json")
        shared["prior/" + label + ".json"] = file
    for name in ("src/tacho_c_api.cpp", "src/qs_wannier90.F", "src/qs_localizer_torus.F",
                 "src/localizer_sparse_pfaffian.F", "src/spectral_localizer_bloch.F", "build-mpi/CMakeCache.txt"):
        shared[name] = root / name
    archives.append(pack(destination / "methods-results.tar.gz", shared, dict(runtime=runtime(root),
        notes="Complete-band sparse reference; CP2K exports use 2 MPI ranks x 2 OpenMP threads; "
        "reference analyses run sequentially with one BLAS/OpenMP thread. "
        "Ritz gaps are estimates, not certified inertia brackets. "
        "Whole-analysis high-water RSS excludes the earlier MPI export.")))
    for name, data in (("index.json", dict(archives=archives, cases=cases())), ("summary.json", summary)):
        with (destination / name).open("x") as stream:
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write("\n")
    print(json.dumps(archives, indent=2))


if __name__ == "__main__":
    main()
