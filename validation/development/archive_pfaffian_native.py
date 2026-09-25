"""Retain native delayed-skew integration, including the diagnosed staging failure."""

import argparse
import json
from pathlib import Path
import subprocess

from archive_stanene_convergence import digest, pack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root, destination = args.source.resolve(), args.destination.resolve()
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "src", "docs", "tests"], cwd=root, check=True)
    results = root / "build-mpi/delayed-native-materials-verified"
    summary = json.loads((results / "summary.json").read_text())
    assert len(summary) == 33 and all(row["after"] is not None for row in summary)
    files = {}
    sources = ["src/tacho_c_api.cpp", "src/localizer_sparse_pfaffian.F", "src/qs_spectral_localizer.F",
               "src/input_cp2k_properties_dft.F", "src/spectral_localizer_sparse_unittest.F",
               "docs/technologies/tacho.md", "tests/QS/regtest-spectral-localizer/neon-common.inc"]
    for name in sources:
        files[name] = root / name
    for name in ("CMakeCache.txt", "delayed-native-final-pretty.log", "delayed-native-pattern-build.log",
                 "delayed-native-complete-unit-mpi2.log", "delayed-native-complete-unit-mpi4.log",
                 "delayed-native-verified-regtests.log", "delayed-native-regtests.log",
                 "delayed-native-materials-verified.log", "delayed-native-bismuth-debug.log",
                 "delayed-native-bismuth-variables.log", "delayed-native-pattern-chain-before.json",
                 "delayed-native-pattern-chain-after.json", "delayed-native-sanitizer-configure.log",
                 "delayed-native-sanitizer-final-build.log", "delayed-native-verified-stress.json",
                 "delayed-native-verified-stress.log", "delayed-native-verified-independent.json",
                 "delayed-native-verified-independent.log"):
        files["build-mpi/" + name] = root / "build-mpi" / name
    for path in sorted(results.iterdir()):
        if path.is_file():
            files[str(path.relative_to(root))] = path
    for name in ("CMakeCache.txt", "delayed-native-build.log", "delayed-native-regtests.log"):
        files["build-serial/" + name] = root / "build-serial" / name
    for name in ("validate_adaptive.py", "validate_delayed.py", "recheck_delayed_materials.py",
                 "check_skew_pattern.py", "native-check/CMakeLists.txt",
                 "build-native-check/CMakeCache.txt"):
        path = root / "build-serial/pfaffian-stability" / name
        files[str(path.relative_to(root))] = path
    for name in ("check_sparse_bloch_localizer.py", "check_bloch_localizer.py", "gaussian_states.py",
                 "check_stanene_localizer.py"):
        files["build-serial/" + name] = root / "build-serial" / name
    fixed = list((root / "build-mpi/delayed-native-regtests").glob("TEST-*/QS/regtest-spectral-localizer-tacho/bismuth-fixed.out"))
    assert len(fixed) == 1
    files["diagnostics/bismuth-fixed.out"] = fixed[0]
    for build, directory in (("build-mpi", "delayed-native-verified-regtests"),
                             ("build-serial", "delayed-native-regtests")):
        for run in (root / build / directory).glob("TEST-*/QS/regtest-spectral-localizer*"):
            for path in run.iterdir():
                if path.suffix in (".inp", ".inc", ".out", ".toml"):
                    files[str(path.relative_to(root))] = path
    metadata = dict(
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        source_sha256={name: digest(root / name) for name in sources},
        previous_index_sha256=digest(destination.parent / "pfaffian-delayed/index.json"),
        material_library_sha256=json.loads((results / "manifest.json").read_text())["library_sha256"],
        final_library_sha256=digest(root / "build-mpi/src/libcp2k.2026.2.dylib"),
        checked_library_sha256=digest(root / "build-serial/pfaffian-stability/build-native-check/libpfaffian_native_checked.dylib"),
        notes="Native serial delayed-skew factors with 4 GiB default numerical-buffer limit. "
        "MPI assembly and MUMPS gap bracketing remain separate. Full 33-query complete-band rerun; "
        "its gaps are independent Ritz estimates, not rigorous MUMPS bounds. The material library "
        "precedes the last unit-driver rebuild, which may change build-info bytes but not numerical sources. "
        "The first integration regtest failed on bismuth because tolerated one-sided roundoff made "
        "the symbolic graph asymmetric. Its log and exact reduced-chain reproducer are retained. "
        "The final implementation uses the real-skew projection with a spectral-norm upper bound "
        "on the discarded part relative to the independent physical gap. No reference values changed. "
        "Sanitizers cover the C++ adapter/extension, not prebuilt dependencies. Timings are not "
        "controlled scaling benchmarks; numerical-buffer limits exclude sparse and dependency storage.")
    destination.mkdir(parents=True, exist_ok=False)
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    with (destination / "index.json").open("x") as handle:
        json.dump({**metadata, "archive": archive}, handle, indent=2)
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
