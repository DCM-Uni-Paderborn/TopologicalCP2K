"""Replay complete finite-Bi material evidence, reintegrating its AO moments."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

from archive_bismuth_regions import compare
from archive_stanene_convergence import digest
from replay_stanene_controls import verify


def extract(archive, destination):
    destination.mkdir()

    def fresh_member(member, root):
        return tarfile.data_filter(member, root).replace(mtime=None)

    with tarfile.open(archive) as stream:
        stream.extractall(destination, filter=fresh_member)


def run(command, log):
    with log.open("x") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"))
    if result.returncode:
        raise RuntimeError(log.read_text())


def reference_archive(bundles, label):
    if not label or Path(label).name != label or label in (".", ".."):
        raise ValueError("Invalid reference case name")
    matches = []
    for bundle in bundles:
        index = json.loads((bundle / "index.json").read_text())
        if label in index["cases"]:
            entry, = [item for item in index["archives"] if item["file"] == label + ".tar.gz"]
            matches.append((bundle / entry["file"], entry["sha256"]))
    if len(matches) != 1:
        raise ValueError("Reference case must occur in exactly one declared bundle: " + label)
    return matches[0]


def replay(bundle, scratch, reference_bundles=()):
    index = json.loads((bundle / "index.json").read_text())
    manifests = {}
    for entry in index["archives"]:
        if Path(entry["file"]).name != entry["file"]:
            raise ValueError("Invalid archive filename")
        manifests[entry["file"]] = verify(bundle / entry["file"], entry["sha256"])
    comparisons, scans, controls = [], [], []
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="bismuth-material-replay-", dir=scratch) as temporary:
        root = Path(temporary)
        methods = root / "methods"
        extract(bundle / "methods-results.tar.gz", methods)
        for label in index["cases"]:
            if Path(label).name != label:
                raise ValueError("Invalid case name")
            work = root / label
            extract(bundle / (label + ".tar.gz"), work)
            expected = json.loads((work / "independent-localizers.json").read_text())
            original_cache = work / "independent-ao-moments.npz"
            retained_cache = work / "retained-ao-moments.npz"
            original_cache.rename(retained_cache)
            energies = list(dict.fromkeys(q["energy_hartree"] for q in expected["queries"]))
            kappas = list(dict.fromkeys(q["kappa_hartree_bohr"] for q in expected["queries"]))
            output = work / "recomputed-analysis.json"
            command = [sys.executable, str(work / "analysis-method.py"), str(methods), str(work),
                       "--energy", *map(str, energies), "--kappa", *map(str, kappas),
                       "--offset", "0", "--output", output.name]
            run(command, work / "recomputed-analysis.log")
            actual = json.loads(output.read_text())
            error = compare(expected, actual)
            with np.load(retained_cache, allow_pickle=False) as old, np.load(original_cache, allow_pickle=False) as new:
                if set(old.files) != set(new.files):
                    raise ValueError("Changed moment-cache fields")
                for key in set(old.files) - {"moments"}:
                    if not np.array_equal(old[key], new[key]):
                        raise ValueError("Changed moment-cache provenance")
                moment_error = float(np.max(abs(old["moments"] - new["moments"])))
                if not np.isfinite(moment_error) or moment_error > 1e-12:
                    raise ValueError("Reintegrated Gaussian moments differ")
            # Scan provenance hashes the original compressed cache, whose ZIP
            # timestamps need not match the freshly generated equivalent arrays.
            os.replace(retained_cache, original_cache)
            comparisons.append(dict(case=label, reintegrated_moment_error=moment_error,
                                    maximum_report_difference=error,
                                    native_queries=len(actual.get("native_comparison", {}).get("queries", []))))
            print(json.dumps(comparisons[-1]), flush=True)
            for report_file in sorted((methods / "additional").glob("*.json")):
                reference = json.loads(report_file.read_text())
                if reference.get("comparison_kind") != "independent-neutral-window-scan" or reference["case"] != label:
                    continue
                output = root / report_file.name
                command = [sys.executable, str(methods / "scripts/scan_bismuth_neutral_window.py"),
                           str(methods), str(work), str(output), "--energy-fraction",
                           *map(str, reference["spectrum"]["energy_gap_fractions"]), "--position-fraction",
                           *map(str, reference["spectrum"]["x_half_width_fractions"]), "--kappa",
                           *map(str, reference["native_arguments"]["kappa"])]
                run(command, output.with_suffix(".log"))
                error = compare(reference, json.loads(output.read_text()))
                scans.append(dict(report=report_file.name, queries=len(reference["queries"]),
                                  maximum_report_difference=error))
                print(json.dumps(scans[-1]), flush=True)
        for report_file in sorted((methods / "additional").glob("*.json")):
            expected = json.loads(report_file.read_text())
            if expected.get("comparison_kind") != "controlled-electronic-temperature":
                continue
            if len(expected["cases"]) != 2:
                raise ValueError("A controlled comparison requires two cases")
            for label in expected["cases"]:
                if not label or Path(label).name != label or label in (".", ".."):
                    raise ValueError("Invalid comparison case name")
                if not (root / label).exists():
                    archive, checksum = reference_archive(reference_bundles, label)
                    manifests[str(archive)] = verify(archive, checksum)
                    extract(archive, root / label)
            script = methods / "scripts/compare_bismuth_basis.py"
            if digest(script) != expected["methods"]["comparison"]:
                raise ValueError("Changed controlled-comparison method")
            output = root / ("replayed-" + report_file.name)
            command = [sys.executable, str(script), str(methods),
                       *(str(root / label) for label in expected["cases"]),
                       str(output), "--control", "temperature"]
            run(command, output.with_suffix(".log"))
            controls.append(dict(report=report_file.name,
                maximum_report_difference=compare(expected, json.loads(output.read_text()))))
            print(json.dumps(controls[-1]), flush=True)
    result = dict(accepted=True, native_scf_rerun=False, ao_moments_reintegrated=True,
                verified_members=sum(len(m["files"]) for m in manifests.values()),
                archive_sha256={entry["file"]: entry["sha256"] for entry in index["archives"]},
                analysis=comparisons, scans=scans, replay_method_sha256=digest(Path(__file__)))
    if controls:
        result.update(temperature_comparisons=controls,
                      reference_archive_sha256={name: digest(Path(name)) for name in manifests
                                                if Path(name).is_absolute()})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--reference-bundles", nargs="*", type=Path, default=[],
                        help="Existing retained bundles used by a controlled temperature comparison")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = replay(args.bundle.resolve(), args.scratch.resolve(),
                    [path.resolve() for path in args.reference_bundles])
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
