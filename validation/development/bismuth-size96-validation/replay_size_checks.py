"""Reintegrate complete Gaussian AO spaces using only the retained size archive."""

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    source, destination = args.archive.resolve(), args.output.resolve()
    integrity = verify(source)
    destination.mkdir(exist_ok=False)
    index = json.loads((source / "index.json").read_text())
    paths = [source / chunk["name"] for chunk in index["chunks"]]
    with io.BufferedReader(ChunkReader(paths)) as joined:
        with tarfile.open(fileobj=joined, mode="r|gz") as archive:
            for member in archive:
                archive.extract(member, destination, filter="data")
    sys.path.insert(0, str(destination))
    import numpy as np
    from scipy.linalg import eigvalsh, svdvals
    import analyze_bismuth_flakes as analysis
    from check_bloch_localizer import aii_skew, skew_sign
    from gaussian_states import primitive_integrals, read_snapshot
    from localizer_region import cover_region

    results = []
    for name in ["flake96-spectrum", "flake96-native", "flake96-transition"]:
        work = destination / name
        expected = json.loads((work / "analysis.json").read_text())
        snapshot = read_snapshot(work / "bismuth.topology")
        moments = analysis.moments(snapshot, primitive_integrals)
        with (work / "independent-ao-moments.npz").open("xb") as stream:
            np.savez_compressed(stream, moments=moments, snapshot_sha256=expected["snapshot_sha256"])
        h, positions, diagnostics = analysis.operators(snapshot, moments)
        reference = expected["queries"]
        assert all(row["energy_hartree"] == reference[0]["energy_hartree"] for row in reference)
        assert all(row["position_bohr"] == reference[0]["position_bohr"] for row in reference)
        actual = analysis.queries(h, positions, [reference[0]["energy_hartree"]],
            [row["kappa_hartree_bohr"] for row in reference],
            [np.asarray(reference[0]["position_bohr"])], aii_skew, skew_sign)
        assert len(actual) == len(reference)
        gap_error = max(abs(a["gap_hartree"]-b["gap_hartree"]) for a, b in zip(actual, reference))
        assert gap_error < 1e-10
        assert all(a["z2"] == b["z2"] for a, b in zip(actual, reference))
        result = {"run": name, "queries": len(actual), "maximum_gap_replay_error_hartree": gap_error,
                  "diagnostics": diagnostics}
        if name != "flake96-spectrum":
            comparison = analysis.native_comparison((work / "output.out").read_text(), actual)
            assert comparison["accepted"] and comparison["tolerance_hartree"] == 2e-8
            result["native_comparison"] = comparison
        results.append(result)
        print(json.dumps(result), flush=True)
        if name == "flake96-spectrum":
            region = json.loads((destination / "flake96-regions/box-0.json").read_text())
            center = snapshot.positions.mean(axis=0)
            n, eye = len(h), np.eye(len(h))
            z0 = positions[0]-center[0]*eye-1j*(positions[1]-center[1]*eye)
            rho = float(svdvals(z0)[0])
            assert abs(rho-region["position_norm_bohr"]) < 1e-10

            def gap(point):
                energy, kappa, x, y = point
                mass = h-energy*eye
                off_diagonal = kappa*(z0-((x-center[0])-1j*(y-center[1]))*eye)
                matrix = np.block([[mass, off_diagonal], [off_diagonal.conj().T, -mass]])
                values = eigvalsh(matrix, subset_by_index=[n-2, n+1], driver="evr")
                return float(min(abs(values)))

            cover = cover_region(region["lower"], region["upper"], rho, center[:2], gap,
                                 margin=1e-10, max_nodes=63)
            assert cover["resolved"] and len(cover["evaluations"]) == len(region["bound"]["evaluations"])
            for actual_bound, expected_bound in zip(cover["evaluations"], region["bound"]["evaluations"]):
                assert np.max(abs(np.asarray(actual_bound["midpoint"])-expected_bound["midpoint"])) < 1e-12
                assert abs(actual_bound["lower_bound"]-expected_bound["lower_bound"]) < 1e-10
            anchor = (np.asarray(region["lower"])+region["upper"])/2
            anchor_result = analysis.queries(h, positions, [anchor[0]], [anchor[1]],
                [np.r_[anchor[2:], center[2]]], aii_skew, skew_sign)[0]
            assert anchor_result["z2"] == region["region_z2"] == 1
            result["first_parameter_box"] = {"resolved": True, "z2": 1,
                "replayed_bound_evaluations": len(cover["evaluations"]),
                "minimum_gap_lower_bound_hartree": min(row["lower_bound"] for row in cover["covered"]),
                "scope": "Replayed norm, complete cover and anchor index; original 17 full-spectrum audits retained"}
            print(json.dumps(result["first_parameter_box"]), flush=True)
        del moments, h, positions, snapshot

    command = [sys.executable, str(destination / "compare_nested_material.py"),
               str(destination / "reference48"), str(destination / "flake96-spectrum"),
               str(destination / "replayed-shared-comparison")]
    subprocess.run(command, cwd=destination, check=True)
    original = json.loads((destination / "nested-size-comparison/comparison.json").read_text())
    replay = json.loads((destination / "replayed-shared-comparison/comparison.json").read_text())
    differences = {}
    for space, fields in original["shared_subspaces"].items():
        for key, value in fields.items():
            difference = abs(replay["shared_subspaces"][space][key]-value)
            assert difference < 1e-10, (space, key, difference)
            differences[space + "/" + key] = difference
    subprocess.run([sys.executable, "-m", "unittest", "-v", "test_nested_subspace", "test_size_archive"],
                   cwd=destination, check=True)
    for name, expected in index["files"].items():
        with (destination / name).open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected, name
    report = {"accepted": True, "integrity": integrity, "runs": results,
              "shared_subspace_replay_errors": differences, "elapsed_seconds": time.monotonic()-start,
              "numpy_version": np.__version__, "python": sys.version,
              "openblas_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
              "scope": "Archive-only AO reintegration and mathematical replay. No SCF or fresh native CP2K execution."}
    (destination / "replay-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"accepted": True, "elapsed_seconds": report["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
