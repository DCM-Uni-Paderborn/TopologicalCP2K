"""Archive the displaced-Si and frozen-density controls, without build products."""
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
    tracked = ["src/qs_gamma2kp.F", "src/qs_localizer_torus.F", "src/qs_ao_translations.F",
               "src/qs_quadratic_pseudospectrum.F", "src/quadratic_pseudospectrum.F",
               "src/qs_band_structure.F", "tests/matchers.py", "tests/TEST_DIRS",
               "docs/methods/properties/quadratic_pseudospectrum.md"]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *tracked], cwd=root, check=True)
    files = set(map(Path, tracked))
    files.update(path.relative_to(root) for path in (root/"tests/QS/regtest-gamma-frozen-density").iterdir())
    for build, name in (("build-serial", "frozen-density-before"),
                        ("build-serial", "frozen-density-after"),
                        ("build-mpi", "frozen-density-after"),
                        ("build-mpi", "silicon-translation-projector-tight")):
        for path in (root/build/name).rglob("*"):
            if path.is_file() and path.suffix in {".out", ".log", ".inp", ".states", ".bs", ".json", ".csv"}:
                files.add(path.relative_to(root))
    for build in ("build-serial", "build-mpi"):
        files.update(Path(build)/name for name in ("frozen-density-regtests.log", "frozen-density-final-build.log"))
        for directory in (root/build/"frozen-density-regtests").glob("TEST-*/QS/regtest-gamma-frozen-density"):
            files.update(path.relative_to(root) for path in directory.iterdir()
                         if path.suffix in {".inp", ".inc", ".toml", ".out"})
    files.update(Path("build-serial")/name for name in (
        "frozen-density-final-pretty.log", "check_frozen_gamma_density.py", "check_silicon_translation.py",
        "check_periodic_translation.py", "gaussian-translation-validation/check_gaussian_translation.py"))
    files.update(Path("data")/name for name in ("BASIS_MOLOPT_UZH", "POTENTIAL_UZH"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "notes": "Before-control and initial displaced-Si runs used the preceding source version. "
                         "The PBE charge-density path is unchanged; all density fixes are rechecked by paired Gamma regressions. "
                         "Raw SCF matrices, separate shifted-ghost overlaps and complex quadratic states define the references. "
                         "Build products and restart checkpoints are excluded.", "files": {}}
    with args.archive.open("xb") as stream, tarfile.open(fileobj=stream, mode="w:gz") as archive:
        for name in sorted(files):
            path = root/name
            assert path.is_file() and path.resolve().is_relative_to(root), name
            with path.open("rb") as source:
                size = path.stat().st_size
                manifest["files"][str(name)] = {"bytes": size, "sha256": hashlib.file_digest(source, "sha256").hexdigest()}
                source.seek(0)
                info = tarfile.TarInfo(name.as_posix())
                info.mode, info.size = 0o644, size
                archive.addfile(info, source)
        data = (json.dumps(manifest, indent=2)+"\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    with args.archive.open("rb") as stream:
        print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size,
                          "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
