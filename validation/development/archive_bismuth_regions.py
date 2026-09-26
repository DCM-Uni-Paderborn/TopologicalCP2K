"""Retain and replay joint finite-Bi parameter regions without rerunning DFT."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

from archive_stanene_convergence import digest, pack
from replay_stanene_controls import verify


def compare(expected, actual, path="", tolerance=1e-10):
    """Check the complete partition, provenance and spectra, except elapsed time."""
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise ValueError("Changed report fields: " + path)
        return max((compare(value, actual[key], path + "/" + key, tolerance)
                    for key, value in expected.items() if key != "elapsed_seconds"), default=0.)
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError("Changed report length: " + path)
        return max((compare(a, b, path + f"/{i}", tolerance)
                    for i, (a, b) in enumerate(zip(expected, actual, strict=True))), default=0.)
    if isinstance(expected, float):
        if not isinstance(actual, (int, float)) or not math.isfinite(actual):
            raise ValueError("Nonfinite or nonnumeric replay: " + path)
        error = abs(expected - actual)
        if not math.isfinite(error) or error > tolerance:
            raise ValueError("Numeric replay mismatch: " + path)
        return error
    if type(expected) is not type(actual) or expected != actual:
        raise ValueError("Changed discrete metadata: " + path)
    return 0.


def replay(bundle):
    index = json.loads((bundle / "index.json").read_text())
    methods_archive = bundle / index["archive"]["file"]
    manifest = verify(methods_archive, index["archive"]["sha256"])
    native_archives = {}
    results = []
    scratch = bundle.parents[2] / ".build"
    scratch.mkdir(exist_ok=True)
    def fresh_member(member, destination):
        # Reproducible archives use epoch mtimes; avoid marking active replay
        # inputs as ancient files to temporary-directory cleanup services.
        return tarfile.data_filter(member, destination).replace(mtime=None)
    with tempfile.TemporaryDirectory(prefix="bismuth-region-replay-", dir=scratch) as temporary:
        root = Path(temporary)
        methods = root / "methods"
        methods.mkdir()
        with tarfile.open(methods_archive) as stream:
            stream.extractall(methods, filter=fresh_member)
        for case in index["cases"]:
            archive = (bundle / case["archive"]).resolve()
            if not archive.is_relative_to(bundle.parent.resolve()):
                raise ValueError("Archive reference outside development evidence")
            if archive not in native_archives:
                native_archives[archive] = verify(archive, case["sha256"])
            expected = json.loads((methods / case["report"]).read_text())
            if expected["case"] != Path(expected["case"]).name:
                raise ValueError("Invalid case name")
            work = root / expected["case"]
            if not work.exists():
                work.mkdir()
                with tarfile.open(archive) as stream:
                    stream.extractall(work, filter=fresh_member)
            actual_file = root / Path(case["report"]).name
            command = [sys.executable, str(methods / "scripts/check_bismuth_region.py"),
                       str(methods), str(work), str(methods / "scans" / (expected["case"] + ".json")),
                       str(methods / "scripts"), str(actual_file), "--kappa",
                       str(expected["lower"][1]), str(expected["upper"][1]),
                       "--energy-half-fraction", str(expected["energy_half_fraction"]),
                       "--position-half-angstrom", str(expected["position_half_angstrom"]),
                       "--max-nodes", str(len(expected["bound"]["evaluations"]))]
            with (root / (actual_file.stem + ".log")).open("x") as stream:
                run = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                     env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"))
            if not actual_file.exists():
                raise RuntimeError("Replay did not create a report:\n" +
                                   (root / (actual_file.stem + ".log")).read_text())
            actual = json.loads(actual_file.read_text())
            if run.returncode != (0 if expected["resolved"] else 1):
                raise ValueError("Unexpected replay exit status")
            results.append(dict(report=case["report"], resolved=actual["resolved"],
                                index=actual["region_z2"],
                                maximum_absolute_difference=compare(expected, actual)))
            print(json.dumps(results[-1]), flush=True)
    return dict(accepted=True, native_scf_rerun=False, results=results,
                verified_members=len(manifest["files"]) + sum(len(m["files"]) for m in native_archives.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--reports", nargs="+", type=Path)
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--output", type=Path, help="Fresh replay-report path instead of bundle/replay.json")
    args = parser.parse_args()
    if args.replay:
        result = replay(args.destination.resolve())
        output = args.output or args.destination / "replay.json"
    else:
        if args.output:
            parser.error("--output is only available with --replay")
        if not args.reports:
            parser.error("--reports is required when collecting evidence")
        root, work, destination = args.root.resolve(), args.work.resolve(), args.destination.resolve()
        destination.mkdir(parents=True, exist_ok=False)
        scripts = Path(__file__).parent
        files = {"scripts/" + name: scripts / name for name in
                 ("check_bismuth_region.py", "localizer_region.py", "test_localizer_region.py",
                  "analyze_bismuth_flakes.py", "archive_bismuth_regions.py",
                  "archive_stanene_convergence.py", "replay_stanene_controls.py")}
        files.update({"build-serial/" + name: root / "build-serial" / name for name in
                      ("gaussian_states.py", "check_bloch_localizer.py")})
        files["method-tests.log"] = work / "method-tests-audited.log"
        cases = []
        for path in args.reports:
            report = json.loads(path.read_text())
            for key, name in (("driver", "check_bismuth_region.py"), ("region", "localizer_region.py")):
                if report["methods"][key] != digest(files["scripts/" + name]):
                    raise ValueError("Changed method: " + name)
            scan = root / "build-mpi/bismuth-neutral-window" / (report["case"].replace("-spectrum", "-frontier") + ".json")
            if digest(scan) != report["scan_sha256"]:
                raise ValueError("Changed underlying scan")
            archive = scripts / "bismuth-size-basis-validation" / (report["case"] + ".tar.gz")
            archived_index = json.loads((archive.parent / "index.json").read_text())
            archive_record = next(item for item in archived_index["archives"] if item["file"] == archive.name)
            native = verify(archive, archive_record["sha256"])
            for name, value in report["inputs"].items():
                if native["files"][name]["sha256"] != value:
                    raise ValueError("Native evidence differs: " + name)
            target = "reports/" + path.name
            if target in files:
                raise ValueError("Duplicate report name")
            files[target] = path
            files[target.replace(".json", ".log")] = path.with_suffix(".log")
            files["scans/" + report["case"] + ".json"] = scan
            cases.append(dict(report=target, archive="../bismuth-size-basis-validation/" + archive.name,
                              sha256=digest(archive)))
        record = pack(destination / "methods-results.tar.gz", files,
                      dict(scope="Frozen finite operators; numerical bounds, not interval arithmetic or size convergence"))
        result = dict(archive=record, cases=cases)
        output = destination / "index.json"
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
