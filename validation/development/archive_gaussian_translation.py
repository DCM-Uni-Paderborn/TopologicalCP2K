"""Archive the finite Gaussian translation checks, without restarts or binaries."""
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
    sources = [Path("src") / name for name in (
        "common/bibliography.F", "input_cp2k_properties_dft.F", "qs_finite_ao.F",
        "qs_quadratic_pseudospectrum.F", "quadratic_pseudospectrum.F",
        "quadratic_pseudospectrum_unittest.F", "quadratic_dbcsr_unittest.F",
        "spectral_localizer_dbcsr.F", "localizer_sparse_inertia.F")]
    sources += [Path("tests/matchers.py"), Path("docs/methods/properties/quadratic_pseudospectrum.md")]
    for name in ("regtest-quadratic-pseudospectrum", "regtest-quadratic-pseudospectrum-mumps"):
        sources += [p.relative_to(root) for p in (root / "tests/QS" / name).iterdir()
                    if p.suffix in {".toml", ".inp", ".inc"}]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *map(str, sources)], cwd=root, check=True)
    files = set(sources)
    for build in ("build-serial", "build-mpi"):
        for name in ("final-build", "final-regtests", "verified"):
            files.add(Path(build) / f"gaussian-translation-{name}.log")
        runs = [root / build / "gaussian-translation-verified"]
        for test in (root / build / "gaussian-translation-final-regtests").glob("TEST-*"):
            for section in ("UNIT", "QS"):
                runs += [path for path in (test / section).iterdir() if
                         ("quadratic" in path.name or "spectral_localizer" in path.name or
                          "spectral-localizer" in path.name)]
        for directory in runs:
            for path in directory.rglob("*"):
                if path.is_file() and (path.suffix in {".out", ".log", ".states", ".inp", ".inc", ".toml", ".json"}
                                       or path.name == "analytic-basis"):
                    files.add(path.relative_to(root))
    for path in (root / "build-mpi/gaussian-translation-extensions").rglob("*"):
        if path.is_file() and path.suffix in {".out", ".log", ".states", ".inp", ".json"}:
            files.add(path.relative_to(root))
    for name in ("check_gaussian_translation.py", "check_translation_extensions.py"):
        files.add(Path("build-serial/gaussian-translation-validation") / name)
    files.update(Path("build-serial") / name for name in (
        "gaussian-translation-debug-build.log", "gaussian-translation-debug-final.log",
        "gaussian-translation-final-pretty.log"))
    files.add(Path("build-mpi/gaussian-translation-extensions-checked.log"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "files": {}}
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for relative in sorted(files):
            path = root / relative
            assert path.is_file() and path.resolve().is_relative_to(root), path
            data = path.read_bytes()
            manifest["files"][str(relative)] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            info = tarfile.TarInfo(relative.as_posix())
            info.mode, info.size = 0o644, len(data)
            archive.addfile(info, io.BytesIO(data))
        data = (json.dumps(manifest, indent=2) + "\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size,
                      "sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
