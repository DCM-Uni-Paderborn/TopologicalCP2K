"""Retain the floating-pivot correction without replacing previous evidence."""

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
    work = root / "build-serial/pfaffian-stability"
    destination.mkdir(parents=True, exist_ok=False)
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "src", "docs", "tools"],
                   cwd=root, check=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    files = {}
    for name in ("src/tacho_c_api.cpp", "src/spectral_localizer_sparse_unittest.F",
                 "docs/technologies/tacho.md", "tools/toolchain/scripts/stage5/tacho-sk-floating-pivot.patch",
                 "build-mpi/CMakeCache.txt"):
        files[name] = root / name
    for name in ("recheck_materials.py", "material-recheck.log", "check_pivot_fixture.py",
                 "resume_materials.py", "material-recheck-resume.log",
                 "pivot-fixture-v2.json", "fixture-v2-baseline.json", "fixture-v2-typed.json",
                 "abs-diagnostic.log", "CMakeLists.txt", "driver.cpp", "tacho_diagnostic.cpp",
                 "headers/Tacho_SkLDL_Internal.hpp", "all-rows-matching-mesh6.log",
                 "max-pivot-mesh6.json", "max-pivot-mesh6.log",
                 "max-pivot-mesh12.json", "max-pivot-mesh12.log",
                 "baseline-local-blas-tzvp.json", "baseline-local-blas-tzvp.log",
                 "patched-homebrew-tzvp.json", "patched-homebrew-tzvp.log",
                 "max-pivot-local-tzvp.json", "max-pivot-local-tzvp.log",
                 "dependency-build.log", "dependency-install.log", "cp2k-build.log",
                 "unit-mpi2.log", "unit-mpi4.log", "regtest-mpi.log", "regtest-sparse-mpi.log", "pretty.log"):
        files["evidence/" + name] = work / name
    reports = sorted((work / "material-recheck").glob("*.json"))
    assert len(reports) == 11
    for path in reports:
        files["reports/" + path.name] = path
    for name in ("check_sparse_bloch_localizer.py", "check_bloch_localizer.py", "gaussian_states.py",
                 "check_stanene_localizer.py"):
        files["scripts/" + name] = root / "build-serial" / name
    files["dependency/Tacho_SkLDL_Internal.hpp"] = (
        root / "build-sparse-deps/Trilinos/packages/shylu/shylu_node/tacho/src/impl/Tacho_SkLDL_Internal.hpp")
    original = destination.parent / "stanene-convergence/index.json"
    metadata = dict(source_commit=commit, original_index_sha256=digest(original),
                    dependency_revision=subprocess.check_output(["git", "rev-parse", "HEAD"],
                        cwd=root / "build-sparse-deps/Trilinos", text=True).strip(),
                    notes="Only typed magnitudes changed in the production dependency. "
                    "Original CP2K matching, pivot rows, threshold and residual checks retained. "
                    "The diagnostic header/candidate-matching experiments are not production code. "
                    "CP2K binaries were built just before this source-only checkpoint commit; "
                    "the archive source and dependency patch identify the actual tested implementation. "
                    "The first batch deliberately failed its monotonic-success assertion on a newly "
                    "unresolved TZVP case. The resume retains that failure, reuses only completed "
                    "hash-checked reports, and finishes the missing cases with all regressions explicit. "
                    "Native material scans use the process-local thread-safe BLAS path; separate "
                    "old/patched BLAS controls retain the near-threshold runtime dependence.")
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    with (destination / "index.json").open("x") as handle:
        json.dump({**metadata, "archive": archive}, handle, indent=2)
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
