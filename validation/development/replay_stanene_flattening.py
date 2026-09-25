"""Check all retained hashes and optionally recompute complete material scans."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--recompute", nargs="*", default=[],
                        choices=["all", "mesh3", "mesh6", "mesh8", "mesh9", "mesh6-tzvp"])
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="stanene-flattening-replay-") as tmp:
        root = Path(tmp)
        with tarfile.open(args.archive) as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root/"manifest.json").read_text())
        for name, expected in manifest["files"].items():
            path = root/name
            assert path.is_file() and path.stat().st_size == expected["bytes"], name
            with path.open("rb") as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"], name
        records, cases, recomputed = [], [], []
        for case in manifest["cases"]:
            work = root/case["directory"]
            stored = json.loads((work/"flattened-scan.json").read_text())
            queries = stored["queries"]
            if "all" in args.recompute or case["label"] in args.recompute:
                command = [sys.executable, str(root/"build-serial/check_bloch_localizer.py"), str(work),
                           "--energy", str(case["energy"]), "--flatten"]
                command += [str(v) for v in sorted({q["flatten"] for q in queries})]
                command += ["--eta"]+[str(v) for v in sorted({q["eta"] for q in queries})]
                command += ["--scale", "1", "--output", "recomputed.json"]
                subprocess.run(command, check=True)
                replay = json.loads((work/"recomputed.json").read_text())
                assert len(replay["queries"]) == len(queries)
                for original, new in zip(queries, replay["queries"], strict=True):
                    for key in ("flatten", "eta", "nu", "polar_coordinates"):
                        assert original[key] == new[key], (case["label"], key)
                    assert abs(original["gap"]-new["gap"]) < 1e-9, case["label"]
                recomputed.append(case["label"])
            diagnostics = stored["diagnostics"]
            assert diagnostics["fourier_frame_unitarity"] < 1e-11
            assert diagnostics["sampled_energy_gap_hartree"] > 0
            assert diagnostics["occupied_spinors"] == 8*case["size"]**2
            for row in queries:
                assert row["gap"] > 0 and row["aii_residual"] < 1e-9 and not row["polar_coordinates"]
                assert row["flatten_scale"] == 1
                records.append({**case, **row})
            cases.append({**case, "diagnostics": diagnostics})
        report = {"verified_files": len(manifest["files"]), "source_commit": manifest["source_commit"],
                  "recomputed_cases": recomputed, "cases": cases, "queries": records}
        args.output.write_text(json.dumps(report, indent=2)+"\n")
        with args.output.with_suffix(".csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["label", "basis", "cutoff", "size", "energy",
                                                       "flatten", "eta", "gap", "nu", "aii_residual"],
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
        print(f"Verified {len(manifest['files'])} files; recomputed {recomputed}; retained {len(records)} queries.")


if __name__ == "__main__":
    main()
