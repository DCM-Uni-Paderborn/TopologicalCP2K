"""Archive completed complex-query tests without binaries or restart files."""

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
        "quadratic_pseudospectrum.F", "quadratic_pseudospectrum_unittest.F",
        "quadratic_dbcsr_unittest.F", "spectral_localizer_dbcsr.F",
        "localizer_sparse_inertia.F", "CMakeLists.txt")]
    sources += [Path("tests/UNIT_TESTS"),
                Path("docs/methods/properties/quadratic_pseudospectrum.md")]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *map(str, sources)],
                   cwd=root, check=True)
    files = set(sources)
    for build in ("build-serial", "build-mpi"):
        for suffix in ("final-build", "final-unit", "final-regtests"):
            files.add(Path(build) / f"translation-quadratic-{suffix}.log")
        run = root / build / "translation-quadratic-final-regtests"
        assert run.is_dir(), run
        selected = (
            "UNIT/quadratic_pseudospectrum_unittest", "UNIT/quadratic_dbcsr_unittest",
            "UNIT/spectral_localizer_unittest", "UNIT/spectral_localizer_sparse_unittest",
            "QS/regtest-quadratic-pseudospectrum", "QS/regtest-quadratic-pseudospectrum-mumps",
            "QS/regtest-spectral-localizer", "QS/regtest-spectral-localizer-mumps",
            "QS/regtest-spectral-localizer-tacho")
        for test_run in run.glob("TEST-*"):
            for name in selected:
                for path in (test_run / name).rglob("*"):
                    if path.is_file() and path.suffix in {".out", ".inp", ".inc", ".toml"}:
                        files.add(path.relative_to(root))
    files.update(Path("build-serial") / name for name in (
        "translation-quadratic-debug-verified.log", "translation-quadratic-pretty.log"))
    for name in ("one-rank", "two-columns", "four-rows", "four-columns"):
        files.add(Path("build-mpi") / f"translation-dbcsr-{name}.log")
    for name in ("regtest-quadratic-pseudospectrum", "regtest-quadratic-pseudospectrum-mumps"):
        files.update(path.relative_to(root) for path in (root / "tests/QS" / name).iterdir()
                     if path.is_file())
    manifest = {
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "files": {},
    }
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for relative in sorted(files):
            path = root / relative
            assert path.resolve().is_relative_to(root) and path.is_file(), path
            data = path.read_bytes()
            manifest["files"][relative.as_posix()] = {
                "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            info = tarfile.TarInfo(relative.as_posix())
            info.mode = 0o644
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        data = (json.dumps(manifest, indent=2) + "\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    print(json.dumps({"files": len(files), "archive_bytes": args.archive.stat().st_size,
                      "sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
