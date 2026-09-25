"""Archive periodized Gaussian translation validation without build products."""
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
    sources = [Path("src") / name for name in (
        "CMakeLists.txt", "qs_ao_translations.F", "qs_localizer_torus.F",
        "qs_finite_ao.F", "qs_gamma2kp.F", "qs_quadratic_pseudospectrum.F",
        "input_cp2k_properties_dft.F", "quadratic_pseudospectrum.F",
        "quadratic_pseudospectrum_unittest.F", "quadratic_dbcsr_unittest.F",
        "spectral_localizer_dbcsr.F", "localizer_sparse_inertia.F")]
    sources += [Path("tests/matchers.py"), Path("docs/methods/properties/quadratic_pseudospectrum.md")]
    for name in ("regtest-quadratic-pseudospectrum", "regtest-quadratic-pseudospectrum-mumps"):
        sources += [p.relative_to(root) for p in (root / "tests/QS" / name).iterdir()
                    if p.suffix in {".toml", ".inp", ".inc"}]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *map(str, sources)], cwd=root, check=True)
    files = set(sources)
    runs = ["build-serial/periodic-translation-reference",
            "build-mpi/periodic-translation-reference",
            "build-serial/periodic-translation-skew",
            "build-serial/periodic-translation-gapw-checked",
            "build-mpi/periodic-translation-gapw-checked",
            "build-serial/periodic-translation-covariance",
            "build-mpi/periodic-translation-skew"]  # Retain the pre-fix GAPW failure.
    for build in ("build-serial", "build-mpi"):
        for test in (root / build / "torus-translation-final-regtests").glob("TEST-*"):
            for category in ("UNIT", "QS"):
                for path in (test / category).iterdir():
                    if re.fullmatch(r"(quadratic.*|spectral_localizer.*|localizer.*|"
                                    r"regtest-(quadratic-pseudospectrum.*|spectral-localizer.*|kubo-transport.*|kp-6))",
                                    path.name):
                        runs.append(str(path.relative_to(root)))
        for suffix in ("final-build", "final-regtests"):
            files.add(Path(build) / f"torus-translation-{suffix}.log")
    suffixes = {".out", ".log", ".states", ".inp", ".inc", ".toml", ".json"}
    for run in runs:
        directory = root / run
        assert directory.is_dir(), directory
        for path in directory.rglob("*"):
            if path.is_file() and path.suffix in suffixes:
                files.add(path.relative_to(root))
    files.update(Path("build-serial") / name for name in (
        "check_periodic_translation.py", "check_translation_covariance.py",
        "torus-translation-pretty.log",
        "gaussian-translation-validation/check_gaussian_translation.py"))
    files.update(Path("data") / name for name in (
        "BASIS_MOLOPT_UZH", "POTENTIAL_UZH", "GTH_SOC_POTENTIALS"))
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
