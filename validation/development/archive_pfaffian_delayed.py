"""Archive the separate delayed-skew implementation and unchanged material controls."""

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
    summary = json.loads((work / "delayed-material-recheck/summary.json").read_text())
    assert len(summary) == 33
    destination.mkdir(parents=True, exist_ok=False)
    files = {}
    names = ["tacho_delayed.cpp", "CMakeLists.txt", "standalone/CMakeLists.txt", "standalone/check_delayed_memory.cpp",
             "validate_adaptive.py", "validate_delayed.py", "recheck_delayed_materials.py",
             "delayed-public-build.log", "delayed-material-recheck.log",
             "delayed-standalone-configure.log", "delayed-standalone-build.log",
             "delayed-memory-build.log", "delayed-memory-check.log"]
    for stem in ("delayed-independent", "delayed-factor-verified", "delayed-blocked-independent",
                 "delayed-scalar-independent", "delayed-checked-stress", "delayed-public-independent",
                 "delayed-mesh6", "delayed-mesh12", "delayed-blocked-mesh6",
                 "delayed-scalar-mesh6", "delayed-scalar-mesh12", "delayed-scalar-mesh12-tzvp",
                 "delayed-standalone-stress"):
        names.extend((stem + ".json", stem + ".log"))
    for name in names:
        files["build-serial/pfaffian-stability/" + name] = work / name
    for path in sorted((work / "headers").rglob("*")):
        if path.is_file():
            files[str(path.relative_to(root))] = path
    for path in sorted((work / "delayed-material-recheck").iterdir()):
        if path.is_file():
            files[str(path.relative_to(root))] = path
    for name in ("check_sparse_bloch_localizer.py", "check_bloch_localizer.py", "gaussian_states.py",
                 "check_stanene_localizer.py"):
        files["build-serial/" + name] = root / "build-serial" / name
    files["build-serial/pfaffian-stability/build/CMakeCache.txt"] = work / "build/CMakeCache.txt"
    material_index = destination.parent / "stanene-convergence/index.json"
    cases = json.loads(material_index.read_text())["cases"]
    for case in cases:
        files["baseline/" + case["label"] + ".json"] = root / case["directory"] / "sparse-final-scan.json"
    metadata = dict(
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        previous_index_sha256=digest(destination.parent / "pfaffian-probes/index.json"),
        material_index_sha256=digest(material_index),
        libraries={name: digest(work / "build" / ("lib" + name + ".dylib")) for name in
                   ("pfaffian_delayed", "pfaffian_delayed_verify", "pfaffian_delayed_blocked",
                    "pfaffian_delayed_scalar", "pfaffian_delayed_public", "pfaffian_delayed_checked")},
        source_sha256={name: digest(work / name) for name in
                       ("tacho_delayed.cpp", "headers/Tacho_SkLDL_Delayed_Research.hpp",
                        "headers/Tacho_SkLDL_Symbolic_Research.hpp")},
        notes="Separate serial skew-specific numerical extension, not installed in CP2K. "
        "The final public and checked targets use only public Tacho symbolic interfaces; "
        "the preceding targets used private read-only driver accessors. The blocked and "
        "unblocked algorithms are retained separately. Full factor reconstruction in earlier "
        "binaries was limited to order 192, extended to 512 for public/checked targets. "
        "Four independent RHS, absolute pivot cutoff and 1e-10 worst-RHS residual threshold "
        "remain unchanged. No operator perturbations or missing-sign substitutions. "
        "Gaps are independent near-zero Ritz estimates with eigen residuals, not rigorous inertia brackets. "
        "Panel counts are scalar entries, not condensed blocks. peak_entries excludes small work arrays, "
        "maps and RHS storage; whole-process time -l measurements are separate. These single-host, "
        "possibly overlapping research timings are not controlled scalability benchmarks. "
        "The extra standalone TZVP mesh12 control uses -0.1592 hartree; the complete 33-query rerun "
        "uses the archived basis-specific energies, including -0.16050731 for TZVP.")
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    with (destination / "index.json").open("x") as handle:
        json.dump({**metadata, "archive": archive}, handle, indent=2)
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
