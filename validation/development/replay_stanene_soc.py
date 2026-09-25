"""Verify corrected stanene evidence and reconstruct the full-band localizer."""

import argparse
import csv
import hashlib
import importlib.util
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
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="stanene-soc-replay-") as tmp:
        root = Path(tmp)
        with tarfile.open(args.archive) as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root/"manifest.json").read_text())
        for name, expected in manifest["files"].items():
            path = root/name
            assert path.is_file() and path.stat().st_size == expected["bytes"], name
            with path.open("rb") as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"], name
        spec = importlib.util.spec_from_file_location("stanene_records", root/"build-serial/summarize_stanene_localizer.py")
        parser_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(parser_module)
        records = []
        for source in sorted(root.glob("build-*/stanene-bloch-localizer/*/input.inp")):
            if "&SPECTRAL_LOCALIZER" not in source.read_text():
                continue
            _, parsed = parser_module.parse_case(source.parent, root)
            for row in parsed:
                row["status"] = "superseded" if source.parent.name == "direct3" else "corrected"
            records.extend(parsed)
        work = root/"build-mpi/stanene-bloch-localizer"
        reconstructed = {}
        for label, name, native in (("before", "mesh3", "direct3"),
                                    ("corrected", "mesh3", "corrected3"),
                                    ("small_regression", "bloch2", None)):
            output = "replayed-"+label+".json"
            command = [sys.executable, str(root/"build-serial/check_bloch_localizer.py"),
                       str(work/name), "--output", output]
            if native:
                command += ["--native", str(work/native/"localizer.bin")]
            subprocess.run(command, check=True)
            reconstructed[label] = json.loads((work/name/output).read_text())
        before = reconstructed["before"]["diagnostics"]
        after = reconstructed["corrected"]["diagnostics"]
        assert before["native_diagonal_difference"] > 0.01
        assert before["native_scalar_difference"] < 1e-10
        assert after["native_diagonal_difference"] < 1e-10
        assert after["native_offdiagonal_difference"] < 1e-12

        def case(path):
            return next(row for row in records if row["case"] == path)

        dense = case("build-mpi/stanene-bloch-localizer/corrected3")
        sparse = case("build-mpi/stanene-bloch-localizer/corrected3-tacho")
        serial = case("build-serial/stanene-bloch-localizer/corrected3")
        small = case("build-mpi/stanene-bloch-localizer/corrected2")
        ref = reconstructed["corrected"]["queries"][0]
        assert dense["nu"] == sparse["nu"] == serial["nu"] == ref["nu"] == 0
        assert abs(dense["gap"]-ref["gap"]) < 1e-10
        assert abs(serial["gap"]-dense["gap"]) < 1e-10
        assert sparse["gap_bounds"][0] <= dense["gap"] <= sparse["gap_bounds"][1]
        assert sparse["solve_residual"] < 1e-10
        ref_small = reconstructed["small_regression"]["queries"][0]
        assert small["nu"] == ref_small["nu"] == 0
        assert abs(small["gap"]-ref_small["gap"]) < 1e-10
        report = {"verified_files": len(manifest["files"]), "reconstructed": reconstructed,
                  "serial_mpi_gap_difference": abs(serial["gap"]-dense["gap"]),
                  "small_regression_gap_difference": abs(small["gap"]-ref_small["gap"]), "queries": records}
        args.output.write_text(json.dumps(report, indent=2)+"\n")
        with args.output.with_suffix(".csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["status", "case", "basis", "cutoff", "ranks", "size", "order",
                                                       "solver", "energy", "eta", "nu", "gap"], extrasaction="ignore",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
        print(f"Verified {len(manifest['files'])} files; reconstructed matrices, Pfaffians and solver checks passed.")


if __name__ == "__main__":
    main()
