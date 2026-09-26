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
    parser.add_argument("--targets", action="store_true")
    args = parser.parse_args()
    root, destination = args.root.resolve(), args.destination.resolve()
    soc = args.soc or args.targets
    prefix = "property-target-groups-final" if args.targets else ("property-soc-groups" if soc else "property-groups")
    expected_comparisons, expected_cases, expected_checks = (36, 52, 33) if args.targets else ((32, 44, 29) if soc else (24, 35, 27))
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root):
        raise ValueError("Native source must be a clean committed checkpoint")
    files, records = {}, {}
    for label in ["validation", "window-controls"] + (["gapw"] if args.targets else []):
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
    if soc:
        if not records["validation"].get("require_soc_groups") or not records["validation"].get("include_stanene_loop"):
            raise ValueError("Missing SOC distribution or displaced-loop checks")
        for name in ["src/qs_property_soc.F", "src/qs_kpoint_operators.F", "tests/matchers.py"]:
            files[name] = root / name
        for build in ["build-serial", "build-mpi"]:
            logfile = root / build / f"{prefix}-kp3-regtests.log"
            if "Status: OK" not in logfile.read_text() or "correct: 16 / 16" not in logfile.read_text():
                raise ValueError("Incomplete shared-operator regressions")
            files[str(logfile.relative_to(root))] = logfile
    if args.targets:
        if not records["validation"].get("require_scalar_validation_groups"):
            raise ValueError("Missing distributed target-frame checks")
        for build in ["build-serial", "build-mpi"]:
            name = f"{build}/{prefix}-build.log"
            files[name] = root / name
        folder = root / f"build-mpi/{prefix}-frozen-negative"
        for backend in ["full", "k290", "spglib"]:
            case = folder / backend
            record = json.loads((case / "run.json").read_text())
            output = (case / "output.out").read_text()
            reference = root / "tests/QS/regtest-topology/stanene-tqc.inp"
            if record["reference_sha256"] != digest(reference):
                raise ValueError("Changed frozen-potential reference")
            files[str(reference.relative_to(root))] = reference
            if record["library_sha256"] != native_libraries(root, "build-mpi") or record["binary_sha256"] != digest(root / "build-mpi/bin/cp2k.psmp"):
                raise ValueError("Changed frozen-potential control binary")
            if backend == "full":
                accepted = record["returncode"] == 0 and "PROGRAM ENDED" in output
            else:
                accepted = (record["returncode"] != 0 and "PROGRAM ENDED" not in output
                            and "Property symmetry does not preserve the frozen H/S eigensystem" in output)
            if not accepted:
                raise ValueError("Incorrect frozen-potential rejection")
            for name, expected in record["files"].items():
                if digest(case / name) != expected:
                    raise ValueError("Changed frozen-potential control")
        for path in folder.rglob("*"):
            if path.is_file():
                files[str(path.relative_to(root))] = path
        for label in ["gapw", "frozen-negative"]:
            name = f"build-mpi/{prefix}-{label}.log"
            files[name] = root / name
    for path in (root / "tests/QS/regtest-property-wilson").iterdir():
        if path.is_file():
            files[str(path.relative_to(root))] = path
    scripts = ["run_property_groups.py", "check_property_group_windows.py", "compare_property_wilson.py",
               "run_property_symmetry_controls.py", "archive_property_groups.py", "replay_property_groups.py",
               "archive_stanene_convergence.py"]
    if args.targets:
        scripts.append("run_property_gapw_groups.py")
    for name in scripts:
        files["methods/" + name] = Path(__file__).with_name(name)
    metadata = dict(source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                    comparisons=expected_comparisons, native_cases=expected_cases, regtest_checks_per_build=expected_checks,
                    shared_operator_checks_per_build=16 if soc else 0,
                    validation_collections=["validation", "gapw"] if args.targets else ["validation"],
                    prefix=prefix,
                    shared_libraries={build: native_libraries(root, build) for build in ["build-serial", "build-mpi"]},
                    scope=("Scalar target H/S checks and SOC orbits distributed; ordered links remain collective." if args.targets else
                           ("Scalar representatives and SOC orbit assembly/validation distributed; ordered links remain collective."
                           if soc else "Scalar representative diagonalization groups only; SOC assembly and ordered links remain collective.")))
    if sum(len(r.get("comparisons", [])) for r in records.values()) != expected_comparisons or sum(len(r["cases"]) for r in records.values()) + (3 if args.targets else 0) != expected_cases:
        raise ValueError("Missing comparisons or native cases")
    destination.mkdir(parents=True, exist_ok=False)
    archive = pack(destination / "methods-results.tar.gz", files, metadata)
    (destination / "index.json").write_text(json.dumps({**metadata, "archive": archive}, indent=2) + "\n")
    print(json.dumps(archive, indent=2))


if __name__ == "__main__":
    main()
