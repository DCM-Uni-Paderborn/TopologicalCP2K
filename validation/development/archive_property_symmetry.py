"""Archive property-symmetry inputs, states, rejected controls and final native sources."""

import argparse
import json
from pathlib import Path
import subprocess

from archive_stanene_convergence import digest, pack
from run_property_symmetry_controls import native_libraries


def files_below(root, directory):
    return {str(p.relative_to(root)): p for p in directory.rglob("*") if p.is_file()
            and "__pycache__" not in p.parts and not p.name.endswith((".wfn", ".bak", ".restart"))
            and ".wfn." not in p.name and ".restart." not in p.name}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root, destination = args.root.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip():
        raise ValueError("Commit the complete native source and regression inputs before archiving")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    archives, comparisons, surface = [], [], []
    shared = {}
    selected = {}
    for build in ["build-serial", "build-mpi"]:
        if "\nStatus: OK\n" not in (root / build / "property-final-regtests.log").read_text():
            raise ValueError("Native regression suite did not finish successfully: " + build)
        tests = sorted((root / build / "property-final-regtests").glob("TEST-*/QS/regtest-property-wilson"))
        if len(tests) != 1:
            raise ValueError("Expected one final regression directory")
        directory = tests[0]
        archives.append(pack(destination / (build + "-regtests.tar.gz"), files_below(root, directory),
                             dict(source_commit=commit, stage="final native regression")))
        selected[build] = directory.relative_to(root)
        for name in ["property-final-regtests.log", "property-kpsym-final.log", "CMakeCache.txt"]:
            shared[build + "/" + name] = root / build / name
    for build, directory in selected.items():
        for material in ["helium", "neon", "stanene"]:
            for backend in ["k290", "spglib"]:
                comparisons.append(dict(name=build + "-" + material + "-" + backend,
                                        reference=str(selected["build-serial"] / (material + "-full")),
                                        candidate=str(directory / (material + "-" + backend))))
    for material in ["helium", "neon", "stanene"]:
        comparisons.append(dict(name="mpi-full-" + material,
                                reference=str(selected["build-serial"] / (material + "-full")),
                                candidate=str(selected["build-mpi"] / (material + "-full"))))
    for stage in ["mirrored", "symmetric", "small", "tight", "ci", "wilson", "surface", "frozen-negative"]:
        work = root / "build-mpi" / ("property-stanene-" + stage)
        for backend in ["full", "k290", "spglib"]:
            directory = work / backend
            run = json.loads((directory / "run.json").read_text())
            for name, expected in run["files"].items():
                if digest(directory / name) != expected:
                    raise ValueError("Changed native control: " + str(directory / name))
            if stage == "frozen-negative":
                output = (directory / "output.out").read_text()
                if run["library_sha256"] != native_libraries(root, "build-mpi"):
                    raise ValueError("Frozen-potential control used a different native library")
                if backend == "full":
                    accepted = run["returncode"] == 0 and "PROGRAM ENDED" in output
                else:
                    accepted = run["returncode"] != 0 and "Invalid property orbit" in output
                if not accepted:
                    raise ValueError("Incorrect outcome of frozen-potential control")
            archives.append(pack(destination / ("stanene-" + stage + "-" + backend + ".tar.gz"),
                                 files_below(root, directory),
                                 dict(stage=stage, binary_sha256=run["binary_sha256"], returncode=run["returncode"])))
        if stage in ["tight", "ci", "wilson"]:
            for backend in ["k290", "spglib"]:
                comparisons.append(dict(name="stanene-" + stage + "-" + backend,
                                        reference=str(work.relative_to(root) / "full/stanene"),
                                        candidate=str(work.relative_to(root) / backend / "stanene")))
        if stage == "surface":
            surface = [dict(name="stanene-surface-" + backend,
                            reference=str(work.relative_to(root) / "full"),
                            candidate=str(work.relative_to(root) / backend)) for backend in ["k290", "spglib"]]
        for p in work.glob("*-comparison.*"):
            shared[str(p.relative_to(root))] = p
    negatives = root / "build-serial/property-final-negatives"
    negative_report = json.loads((negatives / "result.json").read_text())
    if not negative_report["accepted"] or len(negative_report["cases"]) != 5:
        raise ValueError("Expected five successful diagnostic controls")
    if negative_report["library_sha256"] != native_libraries(root, "build-serial"):
        raise ValueError("Negative controls used a different native library")
    for case in negative_report["cases"]:
        for name, expected in case["files"].items():
            if digest(negatives / case["name"] / name) != expected:
                raise ValueError("Changed negative control")
    archives.append(pack(destination / "negative-controls.tar.gz", files_below(root, negatives),
                         dict(source_commit=commit, stage="final diagnostic failures")))
    paths = subprocess.check_output(["git", "diff-tree", "--no-commit-id", "--name-only", "-r", "HEAD"],
                                    cwd=root, text=True).splitlines()
    for name in paths:
        shared[name] = root / name
    for name in ["src/kpoint_k_r_trafo_simple.F", "src/kpoint_types.F", "src/qs_basis_rotation_methods.F",
                 "src/topology_symmetry.F", "build-serial/gaussian_states.py",
                 "build-serial/property-comparison-tests.log", "build-serial/property-symmetry-pretty-final.log",
                 "build-serial/property-default-stack-regtests.log", "build-serial/property-diamond-stack.txt",
                 "build-serial/property-input-pretty.log", "data/BASIS_MOLOPT_UZH", "data/POTENTIAL_UZH",
                 "data/GTH_SOC_POTENTIALS", "tests/QS/regtest-topology/helium.nnkp",
                 "tests/QS/regtest-topology/stanene-tqc.inp"]:
        shared[name] = root / name
    for name in ["archive_property_symmetry.py", "archive_stanene_convergence.py", "replay_property_symmetry.py",
                 "compare_property_wilson.py", "compare_property_surface.py", "test_property_wilson.py",
                 "run_property_symmetry_controls.py", "run_property_symmetry_negatives.py"]:
        shared["validation-methods/" + name] = Path(__file__).with_name(name)
    archives.append(pack(destination / "methods.tar.gz", shared, dict(source_commit=commit,
        binaries={b: digest(root / b / "bin" / ("cp2k.ssmp" if b == "build-serial" else "cp2k.psmp"))
                  for b in selected},
        native_libraries={b: native_libraries(root, b) for b in selected},
        regtest_environment=dict(OMP_STACKSIZE="64M", OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1"))))
    report = dict(archives=archives, source_commit=commit, comparisons=comparisons, surfaces=surface,
                  note="Precursor controls retain executable-only hashes, not contemporaneous shared-library hashes. "
                       "Final source, library hashes and regression inputs are recorded separately. "
                       "The first serial suite lacked an explicit OpenMP stack, timed out in grid collocation, "
                       "and was stopped before restarting in a new directory with OMP_STACKSIZE=64M. "
                       "The 144-point Wilson controls test reconstruction but are not sampling-converged.")
    (destination / "index.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(archives, indent=2), flush=True)


if __name__ == "__main__":
    main()
