"""Verify hashes and replay independent periodic AO checks without running CP2K."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="torus-translation-replay-") as directory:
        root = Path(directory).resolve()
        with tarfile.open(args.archive, "r:gz") as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root / "manifest.json").read_text())
        for name, item in manifest["files"].items():
            path = root / name
            assert path.resolve().is_relative_to(root)
            data = path.read_bytes()
            assert len(data) == item["bytes"]
            assert hashlib.sha256(data).hexdigest() == item["sha256"], name
        runs = {
            "build-serial/periodic-translation-reference": [],
            "build-mpi/periodic-translation-reference": ["--ranks", "2", "--solver", "ITERATIVE"],
            "build-serial/periodic-translation-skew": ["--skew"],
            "build-serial/periodic-translation-gapw-checked": ["--skew", "--method", "GAPW"],
            "build-mpi/periodic-translation-gapw-checked": ["--skew", "--method", "GAPW", "--ranks", "4", "--solver", "ITERATIVE"],
        }
        results = {}
        for name, options in runs.items():
            output = root / name
            subprocess.run([sys.executable, str(root / "build-serial/check_periodic_translation.py"),
                            str(output), "--replay", *options], check=True)
            results[name] = json.loads((output / "summary.json").read_text())
        diagnosis = subprocess.run(
            [sys.executable, str(root / "build-serial/check_periodic_translation.py"),
             str(root / "build-mpi/periodic-translation-skew"), "--replay", "--skew",
             "--method", "GAPW", "--ranks", "4", "--solver", "ITERATIVE"],
            capture_output=True, text=True)
        assert diagnosis.returncode != 0 and "AssertionError" in diagnosis.stderr
        failure = json.loads(diagnosis.stdout)
        assert 2.6e-5 < failure["spectrum_error"] < 2.8e-5
        assert 1.6e-4 < failure["true_residual"] < 1.8e-4
        results["pre_fix_gapw_failure"] = failure
        output = root / "build-serial/periodic-translation-covariance"
        subprocess.run([sys.executable, str(root / "build-serial/check_translation_covariance.py"),
                        str(output), "--reference", str(root / "build-serial/periodic-translation-reference"),
                        "--replay"], check=True)
        results["covariance"] = json.loads((output / "summary.json").read_text())
        spec = importlib.util.spec_from_file_location("matrix_reference", root / "build-serial/gaussian-translation-validation/check_gaussian_translation.py")
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        directories = []
        for build, suffix in (("build-serial", ""), ("build-mpi", "-mumps")):
            tests = list((root / build / "torus-translation-final-regtests").glob("TEST-*"))
            assert len(tests) == 1
            directories.append(tests[0] / "QS" / f"regtest-quadratic-pseudospectrum{suffix}")
            log = (root / build / "torus-translation-final-regtests.log").read_text()
            assert "Number of FAILED  tests 0" in log and "Number of WRONG   tests 0" in log
            results[build + "-regtests"] = next(line for line in log.splitlines() if line.startswith("Summary:"))
        comparisons = {}
        for name in ("gpw", "gapw", "gapw-supercell", "lattice", "soc", "chain"):
            logs = [(d / f"torus-translation-{name}.inp.out").read_text() for d in directories]
            assert all("PROGRAM ENDED AT" in log for log in logs)
            gaps = [reference.values(log, "Gap [hartree]") for log in logs]
            assert gaps[0].shape == gaps[1].shape == (3, 1)
            difference = float(np.max(abs(gaps[0] - gaps[1])))
            assert difference < 1e-8
            residual = float(np.max(reference.values(logs[1], "Eigenpair residual [hartree^2]")))
            assert residual < 1.01e-10
            comparisons[name] = {"gap_difference": difference, "iterative_residual": residual}
        results["regression_comparisons"] = comparisons
        results["source_commit"] = manifest["source_commit"]
        results["files_verified"] = len(manifest["files"])
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
