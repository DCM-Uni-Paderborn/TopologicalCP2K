"""Archive retained periodic-quadratic records, excluding builds and restarts."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="CP2K checkout with retained validation runs")
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    root = args.source.resolve()
    files = set()
    for directory in ["build-serial/periodic-quadratic-smokes", "build-mpi/periodic-quadratic-smokes",
                      "build-serial/periodic-quadratic-variants-checked"]:
        for path in (root / directory).rglob("*"):
            if path.is_file() and path.suffix in {".inp", ".inc", ".out", ".log", ".states", ".json"}:
                files.add(path.relative_to(root))
    for directory in ["tests/QS/regtest-quadratic-pseudospectrum", "tests/QS/regtest-quadratic-pseudospectrum-mumps"]:
        files.update(path.relative_to(root) for path in (root / directory).iterdir() if path.is_file())
    for name in ["quadratic_pseudospectrum.F", "quadratic_pseudospectrum_unittest.F",
                 "qs_quadratic_pseudospectrum.F", "qs_localizer_torus.F", "spectral_localizer_dbcsr.F",
                 "input_cp2k_properties_dft.F"]:
        files.add(Path("src") / name)
    for name in ["run_periodic_quadratic.py", "validate_periodic_quadratic_variants.py"]:
        files.add(Path("build-serial/compatibility-oracle") / name)
    for name in ["BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS"]:
        files.add(Path("data") / name)
    for directory in ["build-serial", "build-mpi"]:
        files.add(Path(directory) / "periodic-quadratic-regtests.log")
        files.add(Path(directory) / "periodic-quadratic-final-build.log")
    for name in ["periodic-quadratic-pretty.log", "periodic-quadratic-debug.log",
                 "periodic-quadratic-debug-ieee.log"]:
        files.add(Path("build-serial") / name)
    files.add(Path("docs/methods/properties/quadratic_pseudospectrum.md"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "files": {}}
    args.archive.parent.mkdir(parents=True, exist_ok=True)
    with args.archive.open("xb") as stream, tarfile.open(fileobj=stream, mode="w:gz") as archive:
        for relative in sorted(files):
            path = root / relative
            assert path.resolve().is_relative_to(root) and path.is_file()
            data = path.read_bytes()
            manifest["files"][relative.as_posix()] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            info = tarfile.TarInfo(relative.as_posix())
            info.size = len(data)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(data))
        data = (json.dumps(manifest, indent=2) + "\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    print(json.dumps({"files": len(files), "uncompressed_bytes": sum(r["bytes"] for r in manifest["files"].values()),
                      "archive_bytes": args.archive.stat().st_size,
                      "sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
