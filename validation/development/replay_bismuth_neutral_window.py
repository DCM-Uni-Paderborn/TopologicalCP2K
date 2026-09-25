"""Replay a neutral-window scan from retained case and Gaussian-method archives."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

from replay_stanene_controls import digest, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_archive", type=Path)
    parser.add_argument("methods_archive", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("scan_method", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    expected = json.loads(args.report.read_text())
    if digest(args.scan_method) != expected["methods"]["scan"]:
        raise ValueError("Different scan implementation")
    if expected["native_comparison_performed"]:
        raise ValueError("This replay only verifies independent neutral scans")
    if Path(expected["case"]).name != expected["case"]:
        raise ValueError("Invalid case label")
    archives = {str(p.resolve()): digest(p) for p in (args.case_archive, args.methods_archive)}
    manifests = [verify(p, archives[str(p.resolve())]) for p in (args.case_archive, args.methods_archive)]
    for name, value in expected["inputs"].items():
        if manifests[0]["files"][name]["sha256"] != value:
            raise ValueError("Report does not match retained input: " + name)
    with tempfile.TemporaryDirectory(prefix="bismuth-neutral-replay-") as temporary:
        root = Path(temporary)
        work = root / expected["case"]
        work.mkdir()
        with tarfile.open(args.case_archive) as archive:
            archive.extractall(work, filter="data")
        methods = root / "methods"
        methods.mkdir()
        with tarfile.open(args.methods_archive) as archive:
            archive.extractall(methods, filter="data")
        scan = root / "scan.py"
        shutil.copyfile(args.scan_method, scan)
        shutil.copyfile(work / "analysis-method.py", root / "analyze_bismuth_flakes.py")
        output = root / "recomputed.json"
        command = [sys.executable, str(scan), str(methods), str(work), str(output),
                   "--energy-fraction", *map(str, expected["spectrum"]["energy_gap_fractions"]),
                   "--position-fraction", *map(str, expected["spectrum"]["x_half_width_fractions"]),
                   "--kappa", *map(str, expected["native_arguments"]["kappa"])]
        with (root / "scan.log").open("w") as stream:
            subprocess.run(command, env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"),
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        actual = json.loads(output.read_text())
        for name in ("case", "comparison_kind", "native_comparison_performed", "spectrum", "methods", "inputs"):
            if actual[name] != expected[name]:
                raise ValueError("Replay metadata mismatch: " + name)
        gap_error, weight_error = 0., 0.
        for old, new in zip(expected["queries"], actual["queries"], strict=True):
            for key in ("energy_hartree", "kappa_hartree_bohr", "position_bohr", "z2"):
                if old[key] != new[key]:
                    raise ValueError("Query or index changed: " + key)
            gap_error = max(gap_error, abs(old["gap_hartree"] - new["gap_hartree"]))
        for old, new in zip(expected["frontier_localization"]["frontier"],
                            actual["frontier_localization"]["frontier"], strict=True):
            if old["bands_one_based"] != new["bands_one_based"]:
                raise ValueError("Frontier degeneracy group changed")
            weight_error = max(weight_error, float(np.max(abs(np.array(old["atom_weights"]) - new["atom_weights"]))))
        if gap_error > 1e-10 or weight_error > 1e-10:
            raise ValueError("Independent replay differs beyond tolerance")
    result = dict(archives=archives, report_sha256=digest(args.report),
                  script_sha256=digest(Path(__file__)), queries=len(expected["queries"]),
                  verified_archive_files=sum(len(m["files"]) for m in manifests),
                  maximum_gap_error_hartree=gap_error, maximum_frontier_weight_error=weight_error,
                  native_comparison_performed=False)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
