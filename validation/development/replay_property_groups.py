"""Verify archived MPI group evidence and repeat its physical frame/link comparisons."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    index = json.loads((args.archive_directory / "index.json").read_text())
    archive = args.archive_directory / index["archive"]["file"]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != index["archive"]["sha256"]:
        raise ValueError("Archive checksum mismatch")
    comparisons = []
    with tempfile.TemporaryDirectory(prefix="property-groups-replay-") as name:
        root = Path(name)
        with tarfile.open(archive) as stream:
            members = stream.getmembers()
            names = [m.name for m in members]
            if len(names) != len(set(names)):
                raise ValueError("Duplicate archive members")
            manifest = json.load(stream.extractfile("manifest.json"))
            if set(names) != set(manifest["files"]) | {"manifest.json"}:
                raise ValueError("Unlisted archive members")
            for member in members:
                path = Path(member.name)
                if not member.isfile() or path.is_absolute() or ".." in path.parts:
                    raise ValueError("Unsafe archive member")
                data = stream.extractfile(member).read()
                if member.name != "manifest.json":
                    expected = manifest["files"][member.name]
                    if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                        raise ValueError("Member checksum mismatch")
                target = root / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
        for label in index.get("validation_collections", ["validation"]):
            folder = f"build-mpi/{index.get('prefix', 'property-groups')}-{label}"
            record = json.loads((root / folder / "result.json").read_text())
            for i, case in enumerate(record["comparisons"]):
                def seed(key):
                    return root / folder / Path(case[key]).parent.name / "states"
                result = root / f"comparison-{label}-{i}.json"
                subprocess.run([sys.executable, str(root / "methods/compare_property_wilson.py"), str(root),
                                str(seed("reference")), str(seed("candidate")), str(result)], check=True,
                               env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1"), timeout=180)
                comparisons.append(json.loads(result.read_text()))
    accepted = len(comparisons) == index["comparisons"] and all(c["accepted"] for c in comparisons)
    args.output.write_text(json.dumps(dict(accepted=accepted, verified_members=len(manifest["files"]),
                                          comparisons=comparisons), indent=2) + "\n")
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
