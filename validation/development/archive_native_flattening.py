"""Archive native spectral-flattening inputs, outputs and source provenance."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("archive", type=Path)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()
    root = args.source.resolve()
    tracked = ["src/input_cp2k_properties_dft.F", "src/qs_spectral_localizer.F",
               "src/qs_localizer_torus.F", "src/spectral_localizer.F",
               "src/spectral_localizer_dbcsr.F", "src/spectral_localizer_sparse.F",
               "src/spectral_localizer_flatten.F", "src/spectral_localizer_flatten_unittest.F",
               "src/iterate_matrix.F", "src/localizer_sparse_inertia.F"]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *tracked], cwd=root, check=True)
    files = set(map(Path, tracked))
    cases = [{"label": label, "reference": reference,
              "directory": "build-mpi/stanene-flatten-native/" + label}
             for label, reference in (("mesh3", "mesh3"), ("mesh3-tacho", "mesh3"),
                                      ("mesh6", "mesh6"), ("mesh6-tzvp", "mesh6-tzvp"))]
    for case in cases:
        work = root / case["directory"]
        output = (work / "output.out").read_text()
        assert "PROGRAM ENDED AT" in output and "*** SCF run converged" in output
        assert json.loads((work / "run.json").read_text())["returncode"] == 0
        files.update(Path(case["directory"]) / name for name in
                     ("input.inp", "output.out", "stdout.log", "run.json", "summary.json"))
    finite = Path("build-mpi/flatten-adapter-validation")
    for case in ("gpw", "gapw", "ae", "mumps", "closed-gap", "bad-metric", "bad-scale"):
        files.update(finite / case / name for name in ("input.inp", "output.out", "stdout.log"))
    files.add(finite / "summary.json")
    files.update(Path("build-serial") / name for name in
                 ("check_flatten_adapter.py", "check_stanene_localizer.py"))
    files.update(Path("data") / name for name in
                 ("BASIS_MOLOPT_UZH", "POTENTIAL_UZH", "GTH_SOC_POTENTIALS"))
    for build in ("build-serial", "build-mpi"):
        files.add(Path(build) / "CMakeCache.txt")
        log = Path(build) / "flatten-native-final-regtests.log"
        text = (root / log).read_text()
        assert "Status: OK" in text and "Number of FAILED  tests 0" in text
        files.add(log)
        for directory in re.findall(r"^>>> (\S+)", text, re.MULTILINE):
            path = Path(directory)
            assert path.resolve().is_relative_to(root)
            files.update(p.relative_to(root) for p in path.glob("*.inp.out"))
    files.update(Path(build) / name for build, name in (
        ("build-serial", "flatten-native-pretty.log"),
        ("build-serial", "flatten-native-unit.log"),
        ("build-mpi", "flatten-native-unit-mpi2.log"),
        ("build-mpi", "flatten-native-unit-mpi4.log")))
    for name in ("regtest-spectral-localizer", "regtest-spectral-localizer-mumps",
                 "regtest-spectral-localizer-tacho"):
        files.update(p.relative_to(root) for p in (root / "tests/QS" / name).iterdir() if p.is_file())
    with args.reference.open("rb") as stream:
        reference_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    manifest = {
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "reference_archive": {"name": args.reference.name, "sha256": reference_hash},
        "cases": cases, "files": {},
        "notes": "Native opt-in adapter. Banners and run.json commits may predate incremental rebuilds. "
                 "The source snapshot is the tested implementation. Reference material values come from "
                 "the separately retained full-band archive, not fitted native outputs. Runs overlap in "
                 "wall time and are not performance benchmarks. No binaries or restart files retained."
    }
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in sorted(files):
            path = root / name
            assert path.is_file() and path.resolve().is_relative_to(root), path
            with path.open("rb") as stream:
                size = path.stat().st_size
                manifest["files"][str(name)] = {
                    "bytes": size, "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}
                stream.seek(0)
                info = tarfile.TarInfo(name.as_posix())
                info.mode, info.size = 0o644, size
                archive.addfile(info, stream)
        data = (json.dumps(manifest, indent=2) + "\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    with args.archive.open("rb") as stream:
        print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size,
                          "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
