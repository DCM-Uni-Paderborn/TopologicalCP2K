"""Verify retained larger-torus evidence and optionally repeat sparse analyses."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def unpack(path, destination, expected):
    assert digest(path) == expected, path
    with tarfile.open(path) as archive:
        archive.extractall(destination, filter="data")
    manifest = json.loads((destination / "manifest.json").read_text())
    for name, record in manifest["files"].items():
        file = destination / name
        assert file.stat().st_size == record["bytes"] and digest(file) == record["sha256"], name
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--recompute", nargs="*", default=[])
    args = parser.parse_args()
    index = json.loads((args.bundle / "index.json").read_text())
    if args.recompute and not args.library:
        parser.error("--recompute requires --library pointing to a Tacho-enabled CP2K shared library")
    selected = {case["label"] for case in index["cases"]}
    if set(args.recompute) - selected - {"all"}:
        parser.error("Unknown recomputation case")
    with tempfile.TemporaryDirectory(prefix="stanene-convergence-replay-") as tmp:
        base = Path(tmp)
        previous = base / "previous"
        prior = unpack(args.bundle / index["previous_reference_archive"], previous,
                       index["previous_reference_sha256"])
        count = len(prior["files"])
        for archive in index["archives"]:
            directory = base / Path(archive["file"]).stem.removesuffix(".tar")
            count += len(unpack(args.bundle / archive["file"], directory, archive["sha256"])["files"])
        methods = base / "methods-results"
        summary, unresolved, comparisons, recomputed_errors = [], [], [], []
        for case in index["cases"]:
            label = case["label"]
            work = base / label if case["new_export"] else previous / case["directory"]
            report = json.loads((methods / "reports" / (label + ".json")).read_text())
            assert report["script_sha256"] == digest(methods / "build-serial/check_sparse_bloch_localizer.py")
            for name, expected in report["inputs"].items():
                assert digest(work / name) == expected, (label, name)
            diagnostics = report["diagnostics"]
            assert diagnostics["occupied_spinors"] == 8 * case["size"]**2
            assert diagnostics["valence_max"] < case["energy"] < diagnostics["conduction_min"]
            for row in report["queries"]:
                assert math.isfinite(row["gap"])
                assert row["gap"] > 1e-8 and row["eigen_residual"] < 1e-8
                assert row["aii_residual"] < 1e-9 and row["factor_attempts"]
                if row["nu"] is None:
                    assert not any(attempt["resolved"] for attempt in row["factor_attempts"])
                    unresolved.append({"label": label, "eta": row["eta"], "gap": row["gap"]})
                else:
                    assert row["nu"] in (0, 1) and row["factor_attempts"][-1]["resolved"]
                    assert row["factor_residual"] <= 1e-10
            if not case["new_export"]:
                dense = json.loads((work / "flattened-scan.json").read_text())
                expected = {row["eta"]: row for row in dense["queries"] if row["flatten"] == 1}
                for row in report["queries"]:
                    difference = abs(row["gap"] - expected[row["eta"]]["gap"])
                    assert difference < 1e-9
                    if row["nu"] is not None:
                        assert row["nu"] == expected[row["eta"]]["nu"]
                    comparisons.append(difference)
            if "all" in args.recompute or label in args.recompute:
                command = [sys.executable, str(methods / "build-serial/check_sparse_bloch_localizer.py"), str(work),
                           "--library", str(args.library.resolve()), "--energy", str(case["energy"]),
                           "--eta", *[str(row["eta"]) for row in report["queries"]], "--output", "recomputed.json"]
                if label == "mesh3":
                    command += ["--dense-check"]
                subprocess.run(command, check=True)
                recomputed = json.loads((work / "recomputed.json").read_text())
                for old, new in zip(report["queries"], recomputed["queries"], strict=True):
                    error = abs(old["gap"] - new["gap"])
                    assert error < 1e-9
                    recomputed_errors.append(error)
                    if old["nu"] is not None and new["nu"] is not None:
                        assert old["nu"] == new["nu"]
                    # An unresolved factor is not upgraded to success by the replay summary.
                    assert (old["nu"] is None) == (new["nu"] is None)
            summary.append({**case, **report})
        native_errors = []
        for label, number in (("mesh3", 3), ("mesh6", 2), ("mesh8-final", 3)):
            work = methods / "native" / label
            output = (work / "output.out").read_text()
            assert "PROGRAM ENDED AT" in output and "*** SCF run converged" in output
            gaps = list(map(float, re.findall(r"SPECTRAL_LOCALIZER\| Gap \[hartree\]:\s+(\S+)", output)))
            etas = list(map(float, re.findall(r"SPECTRAL_LOCALIZER\| Eta \[hartree\]:\s+(\S+)", output)))
            indices = list(map(int, re.findall(r"SPECTRAL_LOCALIZER\| Z2 index:\s+(\d+)", output)))
            assert len(gaps) == len(etas) == len(indices) == number
            reference_label = label.removesuffix("-final")
            case = next(c for c in index["cases"] if c["label"] == reference_label)
            dense = json.loads((previous / case["directory"] / "flattened-scan.json").read_text())
            expected = {row["eta"]: row for row in dense["queries"] if row["flatten"] == 1}
            for eta, gap, invariant in zip(etas, gaps, indices, strict=True):
                error = abs(gap - expected[eta]["gap"])
                assert error < 1e-9 and invariant == expected[eta]["nu"]
                native_errors.append(error)
        failed = methods / "native/mesh8"
        assert json.loads((failed / "run.json").read_text())["returncode"] == 139
        assert "SIGSEGV" in (failed / "stdout.log").read_text()
        result = {"verified_files": count, "cases": summary, "unresolved_pfaffians": unresolved,
                  "maximum_dense_reference_gap_difference": max(comparisons),
                  "maximum_native_reference_gap_difference": max(native_errors),
                  "maximum_recomputed_gap_difference": max(recomputed_errors, default=None),
                  "recomputed_library_sha256": digest(args.library) if args.recompute else None,
                  "recomputed": args.recompute}
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        print(f"Verified {count} members, {len(summary)} cases, {len(unresolved)} unresolved Pfaffians.")
        print(f"Largest gap differences: sparse/dense {max(comparisons):.3e}; native/reference {max(native_errors):.3e} Ha.")


if __name__ == "__main__":
    main()
