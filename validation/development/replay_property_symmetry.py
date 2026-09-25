"""Verify every archived byte and rerun property comparisons without live native outputs."""

import argparse
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
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    index = json.loads((args.archive / "index.json").read_text())
    records, count = [], 0
    with tempfile.TemporaryDirectory(prefix="property-symmetry-replay-") as temporary:
        root = Path(temporary)
        for entry in index["archives"]:
            p = args.archive / entry["file"]
            if hashlib.sha256(p.read_bytes()).hexdigest() != entry["sha256"]:
                raise ValueError("Archive hash mismatch")
            with tarfile.open(p) as archive:
                manifest = json.load(archive.extractfile("manifest.json"))
                members = [m for m in archive.getmembers() if m.name != "manifest.json"]
                if len(members) != len(manifest["files"]) or len({m.name for m in members}) != len(members):
                    raise ValueError("Duplicate or unlisted archive member")
                for member in members:
                    name = Path(member.name)
                    if not member.isfile() or name.is_absolute() or ".." in name.parts:
                        raise ValueError("Unsafe archive member")
                    data = archive.extractfile(member).read()
                    expected = manifest["files"][member.name]
                    if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                        raise ValueError("Member hash mismatch")
                    destination = root / name
                    if destination.exists() and destination.read_bytes() != data:
                        raise ValueError("Conflicting archive members")
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(data)
                    count += 1
        methods = root / "validation-methods"
        subprocess.run([sys.executable, "-m", "unittest", "test_property_wilson", "-v"], cwd=methods, check=True)
        for case in index["comparisons"] + index["surfaces"]:
            output = root / (case["name"] + ".json")
            is_surface = case in index["surfaces"]
            script = "compare_property_surface.py" if is_surface else "compare_property_wilson.py"
            command = [sys.executable, str(methods / script)]
            if not is_surface:
                command += [str(root)]
            command += [str(root / case["reference"]), str(root / case["candidate"]), str(output)]
            subprocess.run(command, cwd=root, check=True, stdout=subprocess.DEVNULL)
            result = json.loads(output.read_text())
            result.pop("sha256")
            result.pop("reference", None)
            result.pop("candidate", None)
            records.append(dict(name=case["name"], **result))
            print(case["name"], result["accepted"], flush=True)
    args.output.write_text(json.dumps(dict(verified_members=count, comparisons=records,
                                          accepted=all(r["accepted"] for r in records)), indent=2) + "\n")


if __name__ == "__main__":
    main()
