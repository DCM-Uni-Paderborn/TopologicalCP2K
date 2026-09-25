"""Verify all separated-control evidence; optionally recompute selected gaps/signs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify(archive, expected):
    assert digest(archive) == expected
    with tarfile.open(archive) as stream:
        manifest = json.load(stream.extractfile("manifest.json"))
        assert set(stream.getnames()) == {"manifest.json", *manifest["files"]}
        for name, record in manifest["files"].items():
            assert stream.getmember(name).size == record["bytes"]
            assert hashlib.file_digest(stream.extractfile(name), "sha256").hexdigest() == record["sha256"]
    return manifest


def member(archive, name, parse=True):
    with tarfile.open(archive) as stream:
        data = stream.extractfile(name).read().decode()
        return json.loads(data) if parse else data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--recompute", nargs="*", default=[])
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    index = json.loads((bundle / "index.json").read_text())
    summaries = json.loads((bundle / "summary.json").read_text())
    methods = bundle / "methods-results.tar.gz"
    if args.recompute and not args.library:
        parser.error("--recompute requires a Tacho-enabled CP2K --library")
    assert not set(args.recompute) - {c["label"] for c in index["cases"] if c["mode"] == "bloch"}
    manifests = {a["file"]: verify(bundle / a["file"], a["sha256"]) for a in index["archives"]}
    comparisons, records, recomputed, residuals, eigen_residuals = [], [], [], [], []
    for case, summary in zip(index["cases"], summaries, strict=True):
        assert case == summary["case"]
        archive = bundle / (case["label"] + ".tar.gz")
        manifest = manifests[archive.name]["files"]
        run = member(archive, "run.json")
        assert run == summary["run"] and run["returncode"] == 0
        assert run["scf_converged"] and run["completed"]
        for name, value in run["outputs"].items():
            assert manifest[name]["sha256"] == value
        for name, value in run["provenance"]["sources"].items():
            if name in manifests[methods.name]["files"]:
                assert manifests[methods.name]["files"][name]["sha256"] == value
        assert manifests[methods.name]["files"]["scripts/run_stanene_controls.py"]["sha256"] == run["provenance"]["runner_sha256"]
        text = member(archive, "output.out", False)
        assert "PROGRAM ENDED AT" in text and "*** SCF run converged" in text
        input_text = member(archive, "input.inp", False)
        assert f"C 0 0 {case['height']}" in input_text
        coords = re.findall(r"Sn\s+0\.(?:3333333333333333|6666666666666667)\s+\S+\s+(\S+)", input_text)
        assert len(coords) == 2 and abs((float(coords[0]) - float(coords[1])) * case["height"] - 0.852) < 1e-12
        if case["mode"] == "reference":
            assert "Wilson surface sampling converged." in text
            assert summary["wilson_index"] == int(re.findall(r"Converged Z2 invariant:\s+(\d+)", text)[0])
            assert summary["wilson_gap_ev"] == float(re.findall(r"Sampled indirect gap \[eV\]:\s+(\S+)", text)[-1])
            continue
        report = member(archive, "analysis.json")
        assert report == summary["analysis"] and report["provenance"] == run["provenance"]
        for name, value in report["inputs"].items():
            assert manifest[name]["sha256"] == value
        d = report["diagnostics"]
        assert d["occupied_spinors"] == 8 * case["size"]**2
        assert abs(d["conduction_min"] - d["valence_max"] - report["sampled_indirect_gap_hartree"]) < 1e-15
        assert len(report["queries"]) == 3
        for q, eta in zip(report["queries"], (.75, 1., 1.25), strict=True):
            assert q["eta"] == eta and q["flatten_scale"] == 1
            assert abs(q["energy"] - (d["valence_max"] + d["conduction_min"]) / 2) < 1e-15
            assert math.isfinite(q["gap"]) and q["gap"] > 1e-8 and q["nu"] in (0, 1)
            assert q["eigen_residual"] < 1e-8 and q["factor_residual"] <= 1e-10 and q["aii_residual"] < 1e-9
            assert q["factor_attempts"][-1]["resolved"]
            residuals.append(q["factor_residual"])
            eigen_residuals.append(q["eigen_residual"])
            records.append(dict(case=case["label"], eta=eta, gap=q["gap"], nu=q["nu"]))
        if case["label"] in ("dzvp-200-scf8", "tzvp-400-scf8"):
            previous = member(methods, "prior/" + case["label"] + ".json")
            for row in report["queries"]:
                old = next(q for q in previous["queries"] if q["eta"] == row["eta"])
                error = abs(row["gap"] - old["gap"])
                assert error < 1e-9 and row["nu"] == old["nu"]
                comparisons.append(error)
        if case["label"] in args.recompute:
            with tempfile.TemporaryDirectory(prefix="stanene-control-replay-") as tmp:
                target = Path(tmp)
                for source, dest in ((methods, target / "methods"), (archive, target / "case")):
                    with tarfile.open(source) as stream:
                        stream.extractall(dest, filter="data")
                command = [sys.executable, str(target / "methods/build-serial/check_sparse_bloch_localizer.py"),
                           str(target / "case"), "--library", str(args.library.resolve()),
                           "--energy", str(report["queries"][0]["energy"]), "--output", "recomputed.json"]
                subprocess.run(command, check=True)
                new = json.loads((target / "case/recomputed.json").read_text())
                for old, row in zip(report["queries"], new["queries"], strict=True):
                    assert old["nu"] == row["nu"] and abs(old["gap"] - row["gap"]) < 1e-9
                    recomputed.append(abs(old["gap"] - row["gap"]))
    result = dict(verified_files=sum(len(m["files"]) for m in manifests.values()),
                  cases=len(index["cases"]), queries=records,
                  maximum_previous_gap_difference=max(comparisons),
                  maximum_factor_residual=max(residuals), maximum_eigen_residual=max(eigen_residuals),
                  recomputed_cases=args.recompute, maximum_recomputed_gap_difference=max(recomputed, default=None))
    with tempfile.TemporaryDirectory(prefix="stanene-table-replay-") as tmp:
        target = Path(tmp)
        shutil.copyfile(bundle / "summary.json", target / "summary.json")
        script = target / "summarize.py"
        script.write_text(member(methods, "scripts/summarize_stanene_controls.py", False))
        subprocess.run([sys.executable, str(script), str(target)], check=True, stdout=subprocess.DEVNULL)
        for name in ("table.tex", "comparisons.json"):
            assert digest(target / name) == digest(bundle / name), name
    result["generated_table_and_comparisons_verified"] = True
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k != "queries"}, indent=2))


if __name__ == "__main__":
    main()
