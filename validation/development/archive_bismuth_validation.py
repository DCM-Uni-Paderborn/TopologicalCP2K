"""Archive completed Bi calculations, retaining their original per-run methods."""

import argparse
import json
from pathlib import Path
import re

from archive_stanene_convergence import digest, pack
from archive_stanene_controls import retain
from replay_stanene_controls import verify


def collect(work):
    run = json.loads((work / "run.json").read_text())
    if not run["completed"] or not run["scf_converged"] or run["returncode"] != 0:
        raise ValueError(f"Incomplete calculation: {work}")
    for name, expected in run["files"].items():
        if digest(work / name) != expected:
            raise ValueError(f"Changed run file: {work / name}")
    if digest(work / "runner.py") != run["provenance"]["runner_sha256"]:
        raise ValueError("Missing original runner provenance")
    files = {name: work / name for name in {"runner.py", "run.json", *run["files"]}}
    record = dict(label=work.name, run=run)
    if run["options"]["mode"] == "wilson":
        text = (work / "output.out").read_text()
        indices = re.findall(r"Converged Z2 invariant:\s+(\d+)", text)
        if "Wilson surface sampling converged." not in text or len(indices) != 1:
            raise ValueError("Wilson sampling not converged")
        record.update(wilson_index=int(indices[0]), wilson_gap_ev=float(re.findall(
            r"Sampled indirect gap \[eV\]:\s+(\S+)", text)[-1]))
    else:
        analysis = json.loads((work / "independent-localizers.json").read_text())
        if analysis["snapshot_sha256"] != run["files"]["bismuth.topology"]:
            raise ValueError("Analysis refers to a different spectrum")
        if run["options"]["mode"] == "localizer" and not analysis["native_comparison"]["accepted"]:
            raise ValueError("Native and independent localizers disagree")
        files.update({name: work / name for name in
                      ("independent-localizers.json", "independent-localizers.log",
                       "independent-ao-moments.npz", "analysis-method.py")})
        if digest(files["analysis-method.py"]) != analysis["source_sha256"]:
            raise ValueError("Missing original analysis method")
        record["analysis"] = analysis
    return files, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--cases", nargs="+", required=True)
    parser.add_argument("--stage", action="store_true")
    args = parser.parse_args()
    root, work, destination = args.root.resolve(), args.work.resolve(), args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    archives, results = [], []
    for label in args.cases:
        files, record = collect(work / label)
        archive = retain(destination / (label + ".tar.gz"), files, dict(label=label))
        verify(destination / archive["file"], archive["sha256"])
        results.append(record)
        archives.append(archive)
        print(json.dumps(archive), flush=True)
    if args.stage:
        return
    shared = {"scripts/" + name: Path(__file__).with_name(name) for name in
              ("run_bismuth_flakes.py", "analyze_bismuth_flakes.py", "test_bismuth_flakes.py",
               "archive_bismuth_validation.py", "archive_stanene_convergence.py",
               "archive_stanene_controls.py", "replay_stanene_controls.py", "run_stanene_controls.py")}
    for name in ("data/BASIS_MOLOPT_UZH", "data/GTH_SOC_POTENTIALS", "build-mpi/CMakeCache.txt",
                 "build-serial/gaussian_states.py", "build-serial/check_bloch_localizer.py",
                 "src/qs_finite_ao.F", "src/qs_spectral_localizer.F", "src/spectral_localizer.F",
                 "src/spectral_localizer_dbcsr.F", "src/localizer_sparse_pfaffian.F",
                 "src/tacho_c_api.cpp", "src/qs_wannier90.F", "src/qs_gamma2kp.F"):
        shared[name] = root / name
    shared["method-tests.log"] = work / "method-tests.log"
    for record in results:
        for name, expected in record["run"]["provenance"]["sources"].items():
            if digest(shared[name]) != expected:
                raise ValueError("Changed basis or pseudopotential")
        if "analysis" in record and digest(shared["build-serial/gaussian_states.py"]) != record["analysis"]["integral_sha256"]:
            raise ValueError("Changed independent Gaussian integrator")
    archives.append(pack(destination / "methods-results.tar.gz", shared,
        dict(notes="Fixed Bi(111) geometry; restricted PBE with post-SCF GTH SOC. "
             "Finite bare flakes are not relaxed or size-converged bulk classifications. "
             "Original runners are retained per run; early runs did not hash the shared library.")))
    verify(destination / archives[-1]["file"], archives[-1]["sha256"])
    for name, data in (("index.json", dict(cases=args.cases, archives=archives)), ("summary.json", results)):
        with (destination / name).open("x") as stream:
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write("\n")


if __name__ == "__main__":
    main()
