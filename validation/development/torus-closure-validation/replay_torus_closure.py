"""Verify the evidence archive and replay eight full-band torus references."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tarfile


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("index", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    index = json.loads(args.index.read_text())
    assert digest(args.archive) == index["archive_sha256"]
    assert args.archive.stat().st_size == index["archive_bytes"]
    args.output.mkdir(exist_ok=False)
    with tarfile.open(args.archive, "r:gz") as archive:
        names = set()
        for member in archive:
            path = Path(member.name)
            assert member.isfile() and not path.is_absolute() and ".." not in path.parts
            assert member.name not in names
            names.add(member.name)
            archive.extract(member, args.output, filter="data")
    manifest = json.loads((args.output / "manifest.json").read_text())
    assert manifest["files"] == index["files"]
    assert manifest["summary"] == index["summary"]
    assert names == set(index["files"]) | {"manifest.json"}
    assert all(digest(args.output / name) == expected for name, expected in index["files"].items())
    sys.path.insert(0, str(args.output.resolve()))
    from validate_native_torus import CELL, ORIGIN, reference

    results = []
    for mode in ["serial", "mpi"]:
        directory = args.output / ("native-torus-" + mode)
        summary = json.loads((directory / "summary.json").read_text())
        assert summary["completed"]
        for case in summary["cases"]:
            name = case["case"]
            origin = ORIGIN + 2 * CELL.sum(axis=0) if name == "shifted" else ORIGIN
            actual = reference(directory / name, origin, name == "flattened")
            error = abs(actual["gap_hartree"] - case["native_gap_hartree"])
            assert error < summary["tolerance_hartree"] == 2e-8
            assert actual["z2"] == case["native_z2"]
            results.append({"mode": mode, "case": name, "native_gap_difference_hartree": error,
                            "reference": actual})
    assert all(digest(args.output / name) == expected for name, expected in index["files"].items())
    report = {"accepted": True, "runs": results, "archive_sha256": index["archive_sha256"],
              "scope": "Archive-only full-band reconstruction, not a new SCF or model calculation"}
    (args.output / "replay-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"accepted": True, "cases": len(results),
                      "maximum_gap_error_hartree": max(r["native_gap_difference_hartree"] for r in results)}))


if __name__ == "__main__":
    main()
