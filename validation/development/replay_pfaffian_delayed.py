"""Verify delayed-skew evidence against immutable inputs and independent dense references."""

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

from replay_stanene_convergence import digest, unpack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    args = parser.parse_args()
    base = args.bundle.resolve()
    index = json.loads((base / "index.json").read_text())
    assert digest(base.parent / "pfaffian-probes/index.json") == index["previous_index_sha256"]
    material_base = base.parent / "stanene-convergence"
    assert digest(material_base / "index.json") == index["material_index_sha256"]
    material = json.loads((material_base / "index.json").read_text())
    previous = base.parent / "stanene-flattening.tar.gz"
    assert digest(previous) == material["previous_reference_sha256"]
    hashes = {r["file"]: r["sha256"] for r in material["archives"]}
    with tempfile.TemporaryDirectory(prefix="pfaffian-delayed-replay-") as temporary:
        work = Path(temporary)
        manifest = unpack(base / index["archive"]["file"], work, index["archive"]["sha256"])
        evidence = work / "build-serial/pfaffian-stability"
        pivot_base = base.parent / "pfaffian-pivot"
        pivot_index = json.loads((pivot_base / "index.json").read_text())
        probe_index = json.loads((base.parent / "pfaffian-probes/index.json").read_text())
        assert digest(pivot_base / "index.json") == probe_index["previous_index_sha256"]
        pivot_manifest = unpack(pivot_base / pivot_index["archive"]["file"], work / "pivot",
                                pivot_index["archive"]["sha256"])
        synthetic = []
        for label, queries, resolved, singular in (("delayed-public-independent", 600, 468, 132),
                                                   ("delayed-checked-stress", 636, 336, 300)):
            report = json.loads((evidence / (label + ".json")).read_text())
            expected = dict(queries=queries, resolved=resolved, singular_queries=singular,
                            nonsingular_unresolved=0, wrong_accepted_signs=0, singular_accepted=0)
            assert report["summary"] == expected
            assert len(report["results"]) == queries
            for row in report["results"]:
                if row["singular"]:
                    assert row["sign"] is None
                else:
                    assert row["sign"] == row["expected"] and row["residual"] <= 1e-10
            synthetic.append(expected)
        if args.library:
            for script, original in (("validate_adaptive.py", "delayed-public-independent"),
                                     ("validate_delayed.py", "delayed-checked-stress")):
                destination = work / (original + ".json")
                subprocess.run([sys.executable, str(evidence / script), str(args.library.resolve()),
                                str(destination)], check=True)
                old = json.loads((evidence / (original + ".json")).read_text())
                new = json.loads(destination.read_text())
                assert old["summary"] == new["summary"]
                for a, b in zip(old["results"], new["results"], strict=True):
                    for key in ("case", "scale", "expected", "sign", "singular"):
                        assert a[key] == b[key]
        controls, dense_errors, recovered, recovered_corrected, memory = [], [], [], [], []
        for case in material["cases"]:
            label = case["label"]
            report = json.loads((evidence / "delayed-material-recheck" / (label + ".json")).read_text())
            baseline = json.loads((work / "baseline" / (label + ".json")).read_text())
            corrected = json.loads((work / "pivot/reports" /
                                    (label + ".json")).read_text())
            assert report["library_sha256"] == index["libraries"]["pfaffian_delayed_public"]
            assert report["script_sha256"] == digest(work / "build-serial/check_sparse_bloch_localizer.py")
            archive_path = material_base / case["archive"] if case["new_export"] else previous
            if case["new_export"]:
                assert digest(archive_path) == hashes[case["archive"]]
            prefix = "" if case["new_export"] else case["directory"] + "/"
            with tarfile.open(archive_path) as archive:
                for name, value in report["inputs"].items():
                    with archive.extractfile(prefix + name) as stream:
                        assert hashlib.file_digest(stream, "sha256").hexdigest() == value
                    assert baseline["inputs"][name] == value
                dense = {} if case["new_export"] else {
                    row["eta"]: row for row in json.load(archive.extractfile(prefix + "flattened-scan.json"))["queries"]
                    if row["flatten"] == 1}
            for old, new, typed in zip(baseline["queries"], report["queries"], corrected["queries"], strict=True):
                for key in ("eta", "energy", "flatten_scale"):
                    assert old[key] == new[key]
                assert new["energy"] == case["energy"]
                assert abs(old["gap"] - new["gap"]) < 1e-10
                assert new["gap"] > 1e-8 and new["eigen_residual"] < 1e-8
                assert typed["eta"] == new["eta"] and typed["energy"] == new["energy"]
                if typed["nu"] is not None:
                    assert typed["nu"] == new["nu"]
                elif new["nu"] is not None:
                    recovered_corrected.append(dict(case=label, eta=new["eta"], nu=new["nu"]))
                if new["nu"] is not None:
                    assert new["factor_residual"] <= 1e-10
                    assert len(new["factor_attempts"]) == 1
                    if old["nu"] is not None:
                        assert old["nu"] == new["nu"]
                    else:
                        recovered.append(dict(case=label, eta=new["eta"], nu=new["nu"]))
                if dense:
                    reference = dense[new["eta"]]
                    error = abs(reference["gap"] - new["gap"])
                    assert error < 1e-10
                    if new["nu"] is not None:
                        assert new["nu"] == reference["nu"]
                    dense_errors.append(error)
                controls.append(dict(case=label, eta=new["eta"], nu=new["nu"], gap=new["gap"],
                                     residual=new["factor_residual"]))
            log = (evidence / "delayed-material-recheck" / (label + ".log")).read_text()
            maximum = re.search(r"(\d+)\s+peak memory footprint", log)
            memory.append(dict(case=label, peak_process_bytes=int(maximum[1]) if maximum else None))
        reconstruction = []
        for name in ("delayed-public-independent.log", "delayed-checked-stress.log"):
            log = (evidence / name).read_text()
            assert "runtime error:" not in log and "Assertion" not in log
            reconstruction += [float(x) for x in re.findall(r"DELAYED complete factor residual=(\S+)", log)]
        assert reconstruction and all(math.isfinite(r) and r < 1e-10 for r in reconstruction)
        assert "PASS: 676000000-entry front rejected before matrix allocation" in (
            evidence / "delayed-memory-check.log").read_text()
        result = dict(verified_files=len(manifest["files"]), synthetic=synthetic,
                      material_queries=len(controls), resolved=sum(q["nu"] is not None for q in controls),
                      dense_comparisons=len(dense_errors), max_dense_gap_error=max(dense_errors),
                      recovered_from_original=recovered, recovered_from_corrected=recovered_corrected,
                      factor_reconstructions=len(reconstruction), max_factor_reconstruction=max(reconstruction),
                      controls=controls, memory=memory, corrected_baseline_verified_files=len(pivot_manifest["files"]),
                      recomputed_synthetic=bool(args.library))
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
        print(json.dumps({k: v for k, v in result.items() if k != "controls"}, indent=2))


if __name__ == "__main__":
    main()
