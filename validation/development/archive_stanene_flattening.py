"""Retain full-band flattening experiments separately from native AO queries."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    root = args.source.resolve()
    cases = [{"label": "mesh3", "directory": "build-mpi/stanene-bloch-localizer/mesh3",
              "basis": "DZVP", "cutoff": 200, "size": 3, "energy": -0.1592}]
    cases += [{"label": f"mesh{n}", "directory": f"build-mpi/stanene-flattening/mesh{n}",
               "basis": "DZVP", "cutoff": 200, "size": n, "energy": -0.1592}
              for n in (6, 8, 9)]
    cases += [{"label": "mesh6-tzvp", "directory": "build-mpi/stanene-flattening/mesh6-tzvp",
               "basis": "TZVP", "cutoff": 400, "size": 6, "energy": -0.16050731}]
    tracked = ["src/qs_localizer_torus.F", "src/qs_spectral_localizer.F", "src/spectral_localizer.F",
               "src/spectral_localizer_dbcsr.F", "src/spectral_localizer_flatten.F",
               "src/spectral_localizer_flatten_unittest.F", "src/iterate_matrix.F",
               "src/qs_wannier90.F", "src/post_scf_bandstructure_utils.F"]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *tracked], cwd=root, check=True)
    files = set(map(Path, tracked))
    for case in cases:
        work = root/case["directory"]
        text = (work/"output.out").read_text()
        assert "PROGRAM ENDED AT" in text and "*** SCF run converged" in text, work
        report = json.loads((work/"flattened-scan.json").read_text())
        assert report["queries"] and all(not row["polar_coordinates"] for row in report["queries"])
        for name in ("input.inp", "output.out", "stdout.log", "run.json", "mesh.nnkp",
                     "stanene.topology", "stanene.mmn", "stanene.eig", "flattened-scan.json"):
            files.add(Path(case["directory"])/name)
    files.update(Path("build-serial")/name for name in (
        "check_stanene_localizer.py", "check_bloch_localizer.py", "gaussian_states.py"))
    files.update(Path("data")/name for name in ("BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS"))
    files.update(Path(build)/"CMakeCache.txt" for build in ("build-serial", "build-mpi"))
    files.update(p.relative_to(root) for p in root.glob("build-*/flatten-validation/*.log"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "cases": cases, "files": {}, "notes":
                "Exact-sign material results are independent NumPy/SciPy analyses of complete CP2K Bloch exports. "
                "No coordinate polar replacement or band truncation is used. The new DBCSR flattening kernel "
                "is tested separately on constructed matrices, not used to obtain these material results. "
                "No production input flag is introduced. CP2K runtime banners can predate incremental rebuilds."}
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in sorted(files):
            path = root/name
            assert path.is_file() and path.resolve().is_relative_to(root), path
            with path.open("rb") as stream:
                size = path.stat().st_size
                manifest["files"][str(name)] = {"bytes": size, "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}
                stream.seek(0)
                info = tarfile.TarInfo(name.as_posix())
                info.mode, info.size = 0o644, size
                archive.addfile(info, stream)
        data = (json.dumps(manifest, indent=2)+"\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    with args.archive.open("rb") as stream:
        print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size,
                          "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
