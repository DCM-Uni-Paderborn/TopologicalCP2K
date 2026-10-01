"""Replay the basis comparison and corrected control, retaining the original failure."""

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import time

from verify_size_checks import ChunkReader, verify


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def extract(bundle, destination, manifest_name):
    index = json.loads((bundle / "index.json").read_text())
    with io.BufferedReader(ChunkReader([bundle / c["name"] for c in index["chunks"]])) as stream:
        with tarfile.open(fileobj=stream, mode="r|gz") as archive:
            for member in archive:
                if member.name == "manifest.json":
                    member.name = manifest_name
                target = destination / member.name
                if target.exists():
                    assert digest(target) == index["files"][member.name], member.name
                else:
                    archive.extract(member, destination, filter="data")
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_archive", type=Path)
    parser.add_argument("basis_archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    bundles = [args.base_archive.resolve(), args.basis_archive.resolve()]
    integrity = [verify(bundle) for bundle in bundles]
    basis_index = json.loads((bundles[1] / "index.json").read_text())
    assert integrity[0]["archive_sha256"] == basis_index["summary"]["base_archive_sha256"]
    destination = args.output.resolve()
    destination.mkdir(exist_ok=False)
    indices = [extract(bundle, destination, name) for bundle, name in zip(
        bundles, ["size-manifest.json", "basis-manifest.json"])]
    sys.path.insert(0, str(destination))
    import numpy as np
    import scipy
    import analyze_bismuth_flakes as analysis
    from check_bloch_localizer import aii_skew, skew_sign
    from gaussian_states import primitive_integrals, read_snapshot

    checks = []
    for run in ["flake96-spectrum", "flake96-tzvp-spectrum", "flake96-tzvp-roundoff"]:
        work = destination / run
        expected = json.loads((work / "analysis.json").read_text())
        snapshot = read_snapshot(work / "bismuth.topology")
        moments = analysis.moments(snapshot, primitive_integrals)
        with (work / "independent-ao-moments.npz").open("xb") as stream:
            np.savez_compressed(stream, moments=moments, snapshot_sha256=expected["snapshot_sha256"])
        h, xyz, diagnostics = analysis.operators(snapshot, moments)
        reference = expected["queries"]
        assert all(row["energy_hartree"] == reference[0]["energy_hartree"] for row in reference)
        assert all(row["position_bohr"] == reference[0]["position_bohr"] for row in reference)
        actual = analysis.queries(h, xyz, [reference[0]["energy_hartree"]],
                                  [row["kappa_hartree_bohr"] for row in reference],
                                  [np.asarray(reference[0]["position_bohr"])], aii_skew, skew_sign)
        assert len(actual) == len(reference)
        gap_error = max(abs(a["gap_hartree"]-b["gap_hartree"]) for a, b in zip(actual, reference))
        assert gap_error < 1e-10
        assert all(a["z2"] == b["z2"] for a, b in zip(actual, reference))
        result = {"run": run, "queries": len(actual), "maximum_gap_replay_error_hartree": gap_error,
                  "diagnostics": diagnostics}
        if run.endswith("roundoff"):
            comparison = analysis.native_comparison((work / "output.out").read_text(), actual)
            assert comparison["accepted"] and comparison["tolerance_hartree"] == 2e-8
            result["native_comparison"] = comparison
        checks.append(result)
        print(json.dumps(result), flush=True)
        del snapshot, moments, h, xyz

    subprocess.run([sys.executable, "compare_material_basis.py", "flake96-spectrum",
                    "flake96-tzvp-spectrum", "replayed-basis-comparison.json"],
                   cwd=destination, check=True)
    expected = json.loads((destination / "basis-comparison.json").read_text())
    actual = json.loads((destination / "replayed-basis-comparison.json").read_text())
    errors = {}
    for key in ["neutral_gaps_hartree", "combined_moment_errors", "embedding_metric_error",
                "projected_hamiltonian_change_hartree", "mean_projected_hamiltonian_shift_hartree",
                "shift_removed_hamiltonian_change_hartree"]:
        error = float(np.max(abs(np.asarray(actual[key])-expected[key])))
        assert error < 1e-10, (key, error)
        errors[key] = error
    assert len(actual["subspaces"]) == len(expected["subspaces"])
    for a, b in zip(actual["subspaces"], expected["subspaces"]):
        for key in ["name", "left_bands_one_based", "right_bands_one_based"]:
            assert a[key] == b[key]
        for key in ["overlap_singular_values", "projector_frobenius_squared", "projector_spectral_distance"]:
            error = float(np.max(abs(np.asarray(a[key])-b[key])))
            assert error < 5e-7, (a["name"], key, error)
            errors[a["name"]+"/"+key] = error
    subprocess.run([sys.executable, "-m", "unittest", "-v", "test_material_basis",
                    "test_basis_covariance", "test_size_archive"], cwd=destination, check=True)
    for index in indices:
        for name, expected_hash in index["files"].items():
            assert digest(destination / name) == expected_hash, name
    for name in ["flake96-tzvp-native", "flake96-tzvp-diagnostic"]:
        failed = json.loads((destination / "failures" / name / "run-record.json").read_text())
        assert not failed["completed"] and failed["returncode"] != 0
    correction = json.loads((destination / "tzvp-roundoff-state.json").read_text())
    assert correction["completed"] and all(row["returncode"] == 0 for row in correction["steps"])
    assert correction["override_sha256"]["localizer_roundoff.F"] == digest(destination / "localizer_roundoff.F")
    assert correction["override_sha256"]["run_tzvp_roundoff.py"] == digest(destination / "run_tzvp_roundoff.py")
    report = {"accepted": True, "integrity": integrity, "runs": checks,
              "basis_comparison_replay_errors": errors, "elapsed_seconds": time.monotonic()-started,
              "python": sys.version, "numpy_version": np.__version__, "scipy_version": scipy.__version__,
              "openblas_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
              "scope": "Archive-only AO reintegration and cross-basis comparison, not a fresh SCF calculation"}
    with (destination / "replay-report.json").open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"accepted": True, "elapsed_seconds": report["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
