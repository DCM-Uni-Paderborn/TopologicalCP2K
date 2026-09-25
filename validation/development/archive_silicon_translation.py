"""Archive the Si Bloch-band and integral-screening validation."""
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
    args = parser.parse_args()
    root = args.source.resolve()
    files = {Path("src") / name for name in (
        "qs_quadratic_pseudospectrum.F", "spectral_localizer_dbcsr.F",
        "quadratic_dbcsr_unittest.F", "qs_localizer_torus.F",
        "qs_ao_translations.F", "qs_gamma2kp.F", "quadratic_pseudospectrum.F",
        "localizer_sparse_inertia.F", "input_cp2k_properties_dft.F",
        "qs_band_structure.F", "kpoint_methods.F", "qs_finite_ao.F",
        "common/physcon.F")}
    files.add(Path("docs/methods/properties/quadratic_pseudospectrum.md"))
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *map(str, files)], cwd=root, check=True)
    runs = [
        "build-serial/silicon-translation-dense/primitive",
        "build-serial/silicon-translation-dense/pristine",
        "build-mpi/silicon-translation-iterative/primitive",
        "build-serial/silicon-translation-tight/primitive",
        "build-serial/silicon-translation-tight/pristine",
        "build-serial/silicon-translation-orbital-tight/primitive",
        "build-serial/silicon-translation-projector-tight/primitive",
        "build-mpi/silicon-translation-projector-tight/primitive",
        "build-mpi/silicon-translation-projector-tight/pristine",
    ]
    for run in runs:
        directory = root / run
        assert directory.is_dir(), directory
        for path in directory.iterdir():
            if path.suffix in {".out", ".log", ".inp", ".states", ".dat", ".bs", ".json", ".csv"}:
                files.add(path.relative_to(root))
    files.add(Path("build-mpi/silicon-translation-projector-tight/summary.json"))
    files.add(Path("build-mpi/silicon-translation-projector-tight/folded-bands.csv"))
    for build in ("build-serial", "build-mpi"):
        files.update(Path(build)/name for name in (
            "silicon-diagnostic-build.log", "silicon-screening-regtests.log"))
        for test in (root/build/"silicon-screening-regtests").glob("TEST-*"):
            for category in ("UNIT", "QS"):
                for directory in (test/category).iterdir():
                    if not re.fullmatch(r"quadratic.*|spectral_localizer.*|localizer.*|regtest-quadratic-pseudospectrum.*", directory.name):
                        continue
                    files.update(path.relative_to(root) for path in directory.rglob("*")
                                 if path.is_file() and path.suffix in {".out", ".inp", ".inc", ".toml"})
    files.update(Path("build-serial")/name for name in (
        "silicon-diagnostic-pretty.log", "silicon-diagnostic-unit.log",
        "check_silicon_translation.py", "check_periodic_translation.py",
        "gaussian-translation-validation/check_gaussian_translation.py"))
    files.update(Path("build-mpi")/name for name in (
        "silicon-diagnostic-unit-2.log", "silicon-diagnostic-unit-4.log"))
    files.update(Path("data")/name for name in ("BASIS_MOLOPT_UZH", "POTENTIAL_UZH"))
    manifest = {
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "notes": "Restart checkpoints and build products are omitted. The input generator can run from atomic guesses. "
                 "Failure-control operator dumps were collected with a temporary dense dump immediately before the unchanged abort.",
        "files": {},
    }
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for relative in sorted(files):
            path = root/relative
            assert path.is_file() and path.resolve().is_relative_to(root), path
            with path.open("rb") as source:
                digest = hashlib.file_digest(source, "sha256").hexdigest()
                size = path.stat().st_size
                source.seek(0)
                manifest["files"][str(relative)] = {"bytes": size, "sha256": digest}
                info = tarfile.TarInfo(relative.as_posix())
                info.mode, info.size = 0o644, size
                archive.addfile(info, source)
        data = (json.dumps(manifest, indent=2)+"\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    with args.archive.open("rb") as source:
        digest = hashlib.file_digest(source, "sha256").hexdigest()
    print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size, "sha256": digest}, indent=2))


if __name__ == "__main__":
    main()
