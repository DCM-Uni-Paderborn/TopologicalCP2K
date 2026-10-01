"""Verify the chunked size-control archive without extracting large files."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile


class ChunkReader(io.RawIOBase):
    def __init__(self, paths):
        super().__init__()
        self.paths = iter(paths)
        self.current = None

    def readable(self):
        return True

    def readinto(self, buffer):
        while True:
            if self.current is None:
                try:
                    self.current = next(self.paths).open("rb")
                except StopIteration:
                    return 0
            size = self.current.readinto(buffer)
            if size:
                return size
            self.current.close()
            self.current = None

    def close(self):
        if self.current is not None:
            self.current.close()
        super().close()


def verify(root):
    index = json.loads((root / "index.json").read_text())
    digest = hashlib.sha256()
    paths = []
    for chunk in index["chunks"]:
        name = chunk["name"]
        assert Path(name).name == name and name not in {".", ".."}, name
        path = root / name
        assert path.is_file() and not path.is_symlink(), name
        assert path.stat().st_size == chunk["bytes"], name
        part_hash = hashlib.sha256()
        with path.open("rb") as stream:
            while data := stream.read(1024**2):
                part_hash.update(data)
                digest.update(data)
        assert part_hash.hexdigest() == chunk["sha256"], name
        paths.append(path)
    assert len(set(paths)) == len(paths)
    assert sum(c["bytes"] for c in index["chunks"]) == index["archive_bytes"]
    assert digest.hexdigest() == index["archive_sha256"]
    hashes, records, names = {}, {}, set()
    with io.BufferedReader(ChunkReader(paths)) as joined:
        with tarfile.open(fileobj=joined, mode="r|gz") as archive:
            for member in archive:
                name = member.name
                relative = Path(name)
                assert member.isfile() and not relative.is_absolute() and ".." not in relative.parts
                assert name not in names, name
                names.add(name)
                with archive.extractfile(member) as stream:
                    if name == "manifest.json":
                        manifest = json.load(stream)
                    elif name.endswith(("/analysis.json", "/run.json")):
                        raw = stream.read()
                        hashes[name] = hashlib.sha256(raw).hexdigest()
                        records[name] = json.loads(raw)
                    else:
                        hashes[name] = hashlib.file_digest(stream, "sha256").hexdigest()
    assert hashes == manifest["files"] == index["files"]
    assert manifest["summary"] == index["summary"]
    for name, record in records.items():
        parent = str(Path(name).parent)
        if name.endswith("/run.json"):
            assert record["completed"] and record["scf_converged"] and record["returncode"] == 0
            for filename, expected in record["files"].items():
                key = parent + "/" + filename
                if key in hashes:
                    assert hashes[key] == expected, key
        else:
            assert record["run_sha256"] == hashes[parent + "/run.json"]
            assert record["snapshot_sha256"] == hashes[parent + "/bismuth.topology"]
            assert all(hashes[method] == expected for method, expected in record["methods"].items())
            if "native_comparison" in record:
                assert record["native_comparison"]["accepted"]
    return {"accepted": True, "verified_files": len(hashes), "archive_sha256": digest.hexdigest(),
            "scope": "Archive integrity and recorded-run consistency, not fresh numerical replay"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory.resolve()), indent=2))
