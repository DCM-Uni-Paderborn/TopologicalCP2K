"""Verify completed cross-platform evidence without a CP2K installation."""

import hashlib
import json
from pathlib import Path
import sys
import tarfile


root = Path(sys.argv[1]).resolve()
results = []
for host in ["spark", "terok"]:
    path = root / (host + "-platform-checks.tar.gz")
    index = json.loads(path.with_suffix(".index.json").read_text())
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == index["sha256"], path
    assert path.stat().st_size == index["bytes"]
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        assert len(names) == len(set(names)), "Duplicate archive member"
        for member in members:
            relative = Path(member.name)
            assert member.isfile() and not relative.is_absolute() and ".." not in relative.parts
        manifest = json.load(archive.extractfile("manifest.json"))
        assert manifest["summary"] == index["summary"]
        assert set(names) == set(manifest["files"]) | {"manifest.json"}
        for name, expected in manifest["files"].items():
            with archive.extractfile(name) as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == expected, name
        qualification = json.load(archive.extractfile("native-qualification-verified/passed.json"))
        assert len(qualification["runs"]) == 10
        assert {(row["ranks"], row["threads"]) for row in qualification["runs"]} == {(1, 2), (4, 1)}
        for row in qualification["runs"]:
            assert row["returncode"] == 0 and manifest["files"][row["log"]] == row["sha256"]
        assert qualification["source_files_verified"] == len(
            json.load(archive.extractfile("native-qualification-verified/source-hashes.json")))
    results.append({"host": host, "verified_files": len(manifest["files"]),
                    "sha256": index["sha256"], "summary": index["summary"]})
print(json.dumps({"accepted": True, "archives": results,
                  "scope": "Archive integrity, scope and recorded test consistency. No fresh CP2K execution."}, indent=2))
