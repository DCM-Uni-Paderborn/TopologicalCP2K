"""Archive checked MPI property frames, negative controls and native regression records."""

import argparse
import json
from pathlib import Path
import subprocess

from archive_stanene_convergence import digest, pack
from run_property_symmetry_controls import native_libraries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--soc", action="store_true")
    args = parser.parse_args()
    root, destination = args.root.resolve(), args.destination.resolve()
    prefix = "property-soc-groups" if args.soc else "property-groups"
    expected_comparisons, expected_cases, expected_checks = (32, 44, 29) if args.soc else (24, 35, 27)
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root):
        raise ValueError("Native source must be a clean committed checkpoint")
    files, records = {}, {}
    for label in ["validation", "window-controls"]:
        folder = root / f"build-mpi/{prefix}-{label}"
        record = json.loads((folder / "result.json").read_text())
        if not record["accepted"] or record["library_sha256"] != native_libraries(root, "build-mpi"):
            raise ValueError("Incomplete checks or changed native library")
        if record["binary_sha256"] != digest(root / "build-mpi/bin/cp2k.psmp"):
            raise ValueError("Changed native executable")
        for case in record["cases"]:
            if not case["accepted"]:
                raise ValueError("Rejected control")
            for name, expected in case["files"].items():
                if digest(folder / case["name"] / name) != expected:
                    raise ValueError("Changed native output")
        for name, expected in record.get("sources", {}).items():
            if digest(Path(name)) != expected:
                raise ValueError("Changed input source")
        for name, expected in record.get("methods", {}).items():
            if digest(Path(name)) != expected:
                raise ValueError("Changed validation method")
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                files[str(path.relative_to(root))] = path
        records[label] = record
    for build in ["build-serial", "build-mpi"]:
        logfile = root / build / f"{prefix}-regtests.log"
        if "Status: OK" not in logfile.read_text() or f"correct: {expected_checks} / {expected_checks}" not in logfile.read_text():
            raise ValueError("Incomplete directory regressions")
        files[str(logfile.relative_to(root))] = logfile
        files[build + "/CMakeCache.txt"] = root / build / "CMakeCache.txt"
    for name in ["src/qs_property_symmetry.F", "src/qs_wannier90.F", "src/input_cp2k_print_dft.F",
                 "src/kpoint_methods.F", "docs/methods/properties/inversion_topology.md",
                 "build-serial/gaussian_states.py", f"build-serial/{prefix}-pretty-final.log",
                 f"build-mpi/{prefix}-validation.log", f"build-mpi/{prefix}-window-controls.log"]:
        files[name] = root / name
    if args.soc:
        if not records["validation"].get("require_soc_groups") or not records["validation"].get("include_stanene_loop"):
            raise ValueError("Missing SOC distribution or displaced-loop checks")
        for name in ["src/qs_property_soc.F", "src/qs_kpoint_operators.F", "tests/matchers.py"]:
            files[name] = root / name
        for build in ["build-serial", "build-mpi"]:
            logfile = root / build / f"{prefix}-kp3-regtests.log"
            if "Status: OK" not in logfile.read_text() or "correct: 16 / 16" not in logfile.read_text():
                raise ValueError("Incomplete shared-operator regressions")
            files[str(logfile.relative_to(root))] = logfile
    for path in (root / "tests/QS/regtest-property-wilson").iterdir():
        if path.is_file():
            files[str(path.relative_to(root))] = path
    scripts = ["run_property_groups.py", "check_property_group_windows.py", "compare_property_wilson.py",
               "run_property_symmetry_controls.py", "archive_property_groups.py", "replay_property_groups.py",
               "archive_stanene_convergence.py"]
    for name in scripts:
        files["methods/" + name] = Path(__file__).with_name(name)
    metadata = dict(source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                    comparisons=expected_comparisons, native_cases=expected_cases, regtest_checks_per_build=expected_checks,
                    shared_operator_checks_per_build=16 if args.soc else 0,
                    prefix=prefix,
                    shared_libraries={build: native_libraries(root, build) for build in ["build-serial", "build-mpi"]},
                    scope=("Scalar representatives and SOC orbit assembly/validation distributed; ordered links remain collective."
                           if args.soc else "Scalar representative diagonalization groups only; SOC assembly and ordered links remain collective."))
    if len(records["validation"]["comparisons"]) != expected_comparisons or sum(len(r["cases"]) for r in records.values()) != expected_cases:
        raise ValueError("Missing comparisons or native cases")
    destination.mkdir(parents=True, exist_ok=False)
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    (destination / "index.json").write_text(json.dumps({**metadata, "archive": archive}, indent=2) + "\n")
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
