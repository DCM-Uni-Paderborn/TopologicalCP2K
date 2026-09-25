"""Recheck native flattening against independent Bloch and finite-AO references."""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import tarfile
import tempfile


def unpack(archive_path, directory):
    directory.mkdir()
    with tarfile.open(archive_path) as archive:
        archive.extractall(directory, filter="data")
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        path = directory / name
        assert path.is_file() and path.stat().st_size == expected["bytes"], name
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"], name
    return manifest


def native_queries(text):
    rows = []
    row = {}
    for line in text.splitlines():
        if "SPECTRAL_LOCALIZER| Eta [hartree]:" in line:
            row = {"eta": float(line.split(":")[-1])}
        elif "SPECTRAL_LOCALIZER| Gap [hartree]:" in line:
            gap = float(line.split(":")[-1])
            if "gap_lower" in row:
                assert row["gap_lower"] == gap
            else:
                row["gap_lower"] = row["gap_upper"] = gap
        elif "SPECTRAL_LOCALIZER| Gap bracket [hartree]:" in line:
            row["gap_lower"], row["gap_upper"] = map(float, line.split(":")[-1].split())
        elif "SPECTRAL_LOCALIZER| Z2 index:" in line:
            row["nu"] = int(line.split(":")[-1])
            assert {"eta", "gap_lower", "gap_upper", "nu"} == row.keys(), row
            rows.append(row)
    assert rows and all(r["gap_upper"] >= r["gap_lower"] > 0 for r in rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="native-flattening-replay-") as tmp:
        root, reference_root = Path(tmp) / "native", Path(tmp) / "reference"
        manifest = unpack(args.archive, root)
        with args.reference.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == manifest["reference_archive"]["sha256"]
        references = unpack(args.reference, reference_root)
        reference_cases = {case["label"]: case for case in references["cases"]}
        material, diagnostics = [], []
        for case in manifest["cases"]:
            work = root / case["directory"]
            text = (work / "output.out").read_text()
            assert "PROGRAM ENDED AT" in text and "*** SCF run converged" in text
            report = json.loads((reference_root / reference_cases[case["reference"]]["directory"] /
                                 "flattened-scan.json").read_text())
            expected = {q["eta"]: q for q in report["queries"] if q["flatten"] == 1}
            electronic = list(map(float, re.search(r"Original electronic gap bracket \[hartree\]:([^\n]+)", text)[1].split()))
            residuals = list(map(float, re.search(r"Flattening residuals:([^\n]+)", text)[1].split()))
            assert len(electronic) == 2 and electronic[1] >= electronic[0] > 1e-8
            assert len(residuals) == 4 and max(residuals) < 1e-9
            reference_gap = report["diagnostics"]["sampled_energy_gap_hartree"]
            assert electronic[0] - 1e-9 <= reference_gap <= electronic[1] + 1e-9
            diagnostics.append({"label": case["label"], "electronic_gap_bracket": electronic,
                                "reference_electronic_gap": reference_gap, "residuals": residuals})
            for row in native_queries(text):
                ref = expected[row["eta"]]
                assert not ref["polar_coordinates"] and ref["flatten_scale"] == 1
                assert row["nu"] == ref["nu"]
                assert row["gap_lower"] - 1e-9 <= ref["gap"] <= row["gap_upper"] + 1e-9
                error = max(abs(row["gap_lower"] - ref["gap"]), abs(row["gap_upper"] - ref["gap"]))
                material.append({"label": case["label"], **row, "reference_gap": ref["gap"],
                                 "maximum_endpoint_difference": error})
        spec = importlib.util.spec_from_file_location("finite_reference", root / "build-serial/check_flatten_adapter.py")
        finite_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(finite_module)
        finite_root = root / "build-mpi/flatten-adapter-validation"
        finite = {name: finite_module.reference((finite_root / name / "output.out").read_text())
                  for name in ("gpw", "gapw", "ae", "mumps")}
        for name in finite:
            text = (finite_root / name / "output.out").read_text()
            assert list(map(int, re.findall(r"SPECTRAL_LOCALIZER\| Chern index:\s*(\S+)", text))) == [0] * 8
            if name == "mumps":
                brackets = [list(map(float, row.split())) for row in
                            re.findall(r"SPECTRAL_LOCALIZER\| Gap bracket \[hartree\]:([^\n]+)", text)]
                assert len(brackets) == 8
                for bounds, query in zip(brackets, finite[name], strict=True):
                    assert bounds[0] <= query["reference_gap"] <= bounds[1]
        for name, message in (("closed-gap", "electronic gap unresolved"),
                              ("bad-metric", "invalid original AO metric"),
                              ("bad-scale", "invalid scale or tolerance")):
            text = (finite_root / name / "output.out").read_text()
            assert message in text and "SPECTRAL_LOCALIZER| Gap [hartree]:" not in text
        counts = {}
        for build in ("build-serial", "build-mpi"):
            text = (root / build / "flatten-native-final-regtests.log").read_text()
            assert "Number of FAILED  tests 0" in text and "Number of WRONG   tests 0" in text
            correct, total = map(int, re.search(r"Summary: correct: (\d+) / (\d+)", text).groups())
            assert correct == total
            counts[build] = correct
        result = {"source_commit": manifest["source_commit"], "verified_native_files": len(manifest["files"]),
                  "verified_reference_files": len(references["files"]), "material_queries": material,
                  "material_diagnostics": diagnostics, "finite_queries": finite,
                  "rejected_controls": ["closed-gap", "bad-metric", "bad-scale"], "regtest_counts": counts}
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        with args.output.with_suffix(".csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(material[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(material)
        print(json.dumps({"material_queries": len(material), "finite_queries": sum(map(len, finite.values())),
                          "maximum_material_endpoint_difference": max(r["maximum_endpoint_difference"] for r in material),
                          "regtest_counts": counts}, indent=2))


if __name__ == "__main__":
    main()
