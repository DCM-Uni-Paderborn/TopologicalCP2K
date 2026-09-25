"""Archive independent solve probes and a separate adaptive-pivot prototype."""

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
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", "src", "docs"], cwd=root, check=True)
    destination.mkdir(parents=True, exist_ok=False)
    files = {}
    for name in ("src/tacho_c_api.cpp", "src/spectral_localizer_sparse_unittest.F", "docs/technologies/tacho.md",
                 "tools/toolchain/scripts/stage5/tacho-sk-floating-pivot.patch", "build-mpi/CMakeCache.txt"):
        files[name] = root / name
    for name in ("validate_adaptive.py", "find_singular_probe_fixture.py", "singular-probe-fixture.json",
                 "singular-probe-fixture.log", "production-probes-before.json", "production-probes-before.log",
                 "production-probes-after.json", "production-probes-after.log",
                 "adaptive-checked-independent.json", "adaptive-checked-independent.log",
                 "adaptive-checked-mesh6.json", "adaptive-checked-mesh6.log", "adaptive-checked-build.log",
                 "independent-probes-mesh3.json", "independent-probes-mesh3.log",
                 "independent-probes-unit-mpi2-final.log", "independent-probes-unit-mpi4.log",
                 "independent-probes-regtests.log", "independent-probes-sparse-regtests.log",
                 "independent-probes-pretty-final.log", "independent-probes-final-build.log",
                 "CMakeLists.txt", "driver.cpp", "tacho_pair_adaptive.cpp", "tacho_adaptive.cpp", "tacho_diagnostic.cpp",
                 "headers/Tacho_Driver.hpp", "headers/Tacho_Driver_Impl.hpp", "headers/Tacho_NumericTools_Base.hpp",
                 "headers/Tacho_SkLDL_Internal.hpp", "headers/Tacho_SkLDL_Supernodes_Serial.hpp",
                 "headers/research_skew_pivot.hpp", "adaptive-early-mesh12.log"):
        files["build-serial/pfaffian-stability/" + name] = work / name
    for name in ("check_sparse_bloch_localizer.py", "check_bloch_localizer.py", "gaussian_states.py", "check_stanene_localizer.py"):
        files["build-serial/" + name] = root / "build-serial" / name
    # Keep quoted includes inside the private header tree reproducible.
    for path in sorted((work / "headers").rglob("*")):
        if path.is_file():
            files[str(path.relative_to(root))] = path
    metadata = dict(
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        previous_index_sha256=digest(destination.parent / "pfaffian-pivot/index.json"),
        production_library_sha256=digest(root / "build-mpi/src/libcp2k.2026.2.dylib"),
        prototype_library_sha256=digest(work / "build/libpfaffian_pair_checked.dylib"),
        notes="Only independent right-hand sides and their tests enter CP2K. The adaptive header/adapter "
        "is an isolated serial research build, not installed in CP2K. Its N=6 and synthetic results use "
        "the same independent probes. The retained interrupted N=12 trace used the preceding A*x-probe "
        "prototype: it was terminated after more than 100 symbolic retries, without an accepted result. "
        "It is a negative diagnostic, not a completed test or a timing benchmark. No operator entries, "
        "gap tolerance or 1e-10 acceptance threshold are relaxed. Native MPI tests ran before the final "
        "comment-only source annotation and checkpoint commit. Source and binary hashes are retained; "
        "the binary commit banner alone is not used as provenance.")
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    with (destination / "index.json").open("x") as handle:
        json.dump({**metadata, "archive": archive}, handle, indent=2)
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
