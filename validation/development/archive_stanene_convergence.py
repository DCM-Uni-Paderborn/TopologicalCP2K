"""Archive larger full-band stanene controls without overwriting earlier evidence."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def pack(path, files, metadata):
    manifest = {**metadata, "files": {}}
    with path.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name, source in sorted(files.items()):
            info = tarfile.TarInfo(name)
            info.mode, info.size = 0o644, source.stat().st_size
            manifest["files"][name] = {"bytes": info.size, "sha256": digest(source)}
            with source.open("rb") as stream:
                archive.addfile(info, stream)
        data = (json.dumps(manifest, indent=2) + "\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root, destination = args.source.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "src"], cwd=root, check=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    cases = [{"label": "mesh3", "directory": "build-mpi/stanene-bloch-localizer/mesh3", "basis": "DZVP",
              "cutoff": 200, "scf_mesh": 8, "size": 3, "energy": -0.1592, "new_export": False}]
    cases += [{"label": f"mesh{n}", "directory": f"build-mpi/stanene-flattening/mesh{n}", "basis": "DZVP",
               "cutoff": 200, "scf_mesh": 8, "size": n, "energy": -0.1592, "new_export": False} for n in (6, 8, 9)]
    cases += [{"label": "mesh6-tzvp", "directory": "build-mpi/stanene-flattening/mesh6-tzvp", "basis": "TZVP",
               "cutoff": 400, "scf_mesh": 8, "size": 6, "energy": -0.16050731, "new_export": False}]
    cases += [{"label": f"mesh{n}", "directory": f"build-mpi/stanene-material-convergence/mesh{n}", "basis": "DZVP",
               "cutoff": 200, "scf_mesh": 8, "size": n, "energy": -0.1592, "new_export": True} for n in (12, 15, 18)]
    cases += [{"label": "mesh12-tzvp", "directory": "build-mpi/stanene-material-convergence/mesh12-tzvp", "basis": "TZVP",
               "cutoff": 400, "scf_mesh": 8, "size": 12, "energy": -0.16050731, "new_export": True},
              {"label": "mesh12-scf12", "directory": "build-mpi/stanene-material-convergence/mesh12-scf12", "basis": "DZVP",
               "cutoff": 200, "scf_mesh": 12, "size": 12, "energy": -0.1592, "new_export": True}]
    sources = ["src/qs_localizer_torus.F", "src/spectral_localizer_bloch.F", "src/qs_spectral_localizer.F",
               "src/spectral_localizer.F", "src/spectral_localizer_unittest.F", "src/spectral_localizer_dbcsr.F", "src/tacho_c_api.cpp",
               "src/localizer_sparse_pfaffian.F", "src/spectral_localizer_flatten_unittest.F",
               "data/BASIS_MOLOPT_UZH", "data/GTH_SOC_POTENTIALS", "build-mpi/CMakeCache.txt"]
    sources += ["build-serial/" + name for name in ("check_sparse_bloch_localizer.py", "check_bloch_localizer.py",
                "check_stanene_localizer.py", "check_sparse_bloch_references.py", "gaussian_states.py",
                "check_complex_cholesky.F90")]
    common = {name: root / name for name in sources}
    for name in ("sparse-reference-comparisons.json", "sparse-reference-comparisons.log",
                 "blas-homebrew-probe.log", "blas-local-probe.log", "local-blas-load.log", "blas-crash-summary.json",
                 "native-mesh8-sample.txt", "native-mesh8-final-sample.txt"):
        common["diagnostics/" + name] = root / "build-mpi/stanene-material-convergence" / name
    for name in ("bloch-flatten-unit-mpi2.log", "bloch-flatten-unit-mpi4.log", "bloch-flatten-regtests.log",
                 "bloch-flatten-sparse-regtests.log"):
        common["build-mpi/" + name] = root / "build-mpi" / name
    for name in ("bloch-flatten-unit.log", "bloch-flatten-regtests.log", "bloch-flatten-pretty-final.log"):
        common["build-serial/" + name] = root / "build-serial" / name
    for build in ("build-serial", "build-mpi"):
        for name in ("aii-permutation-unit.log", "aii-permutation-regtests.log"):
            common[build + "/" + name] = root / build / name
    common["build-serial/aii-permutation-pretty-final.log"] = root / "build-serial/aii-permutation-pretty-final.log"
    common["build-mpi/aii-permutation-sparse-unit.log"] = root / "build-mpi/aii-permutation-sparse-unit.log"
    archives = []
    summary = []
    raw_names = ("input.inp", "output.out", "stdout.log", "run.json", "mesh.nnkp",
                 "stanene.topology", "stanene.mmn", "stanene.eig")
    for case in cases:
        work = root / case["directory"]
        report = json.loads((work / "sparse-final-scan.json").read_text())
        assert report["script_sha256"] == digest(root / "build-serial/check_sparse_bloch_localizer.py")
        for name, expected in report["inputs"].items():
            assert digest(work / name) == expected
        assert "PROGRAM ENDED AT" in (work / "output.out").read_text()
        common["reports/" + case["label"] + ".json"] = work / "sparse-final-scan.json"
        common["reports/" + case["label"] + ".log"] = work / "sparse-final-scan.log"
        summary.append({**case, **report})
        if case["new_export"]:
            archive = pack(destination / (case["label"] + ".tar.gz"),
                           {name: work / name for name in raw_names}, {"case": case, "source_commit": commit})
            case["archive"] = archive["file"]
            archives.append(archive)
    for case in ("mesh3", "mesh6", "mesh8", "mesh8-local-blas", "mesh8-final"):
        work = root / "build-mpi/stanene-bloch-flatten-native" / case
        for name in ("input.inp", "output.out", "stdout.log", "run.json", "summary.json"):
            common["native/" + case + "/" + name] = work / name
    common["native/mesh8-final/binary-provenance.json"] = root / "build-mpi/stanene-bloch-flatten-native/mesh8-final/binary-provenance.json"
    metadata = {"source_commit": commit, "cases": cases,
                "tacho_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                    cwd=root / "build-sparse-deps/Trilinos", text=True).strip(),
                "notes": "Sparse complete-band reference, not a new production localizer basis. "
                "All unresolved Pfaffian attempts and the native mixed-BLAS failure are retained. "
                "The first local-BLAS native N=8 run was terminated after profiling to replace the "
                "cubic dense AII basis conversion by an exactly equivalent quadratic-cost conversion; "
                "only mesh8-final is a completed native N=8 test. It ran with the new conversion, "
                "before the additional finite-value output guard was added. Final unit and regression "
                "logs cover that guard. Reported sparse-library hashes predate the dense-conversion "
                "change; the sparse factorization implementation is unchanged."}
    archives.append(pack(destination / "methods-results.tar.gz", common, metadata))
    index = {**metadata, "archives": archives,
             "previous_reference_archive": "../stanene-flattening.tar.gz",
             "previous_reference_sha256": digest(destination.parent / "stanene-flattening.tar.gz")}
    (destination / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(archives, indent=2))


if __name__ == "__main__":
    main()
