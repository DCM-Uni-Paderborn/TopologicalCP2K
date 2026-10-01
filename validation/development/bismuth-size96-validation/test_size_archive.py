"""Small positive and corruption controls for the streaming evidence reader."""

import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from verify_size_checks import verify


class ArchiveTests(unittest.TestCase):
    def fixture(self, directory):
        files = {"example.txt": b"finite-size archive fixture\n"}
        hashes = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
        summary = {"scope": "Reader test, not physical validation"}
        files["manifest.json"] = json.dumps({"files": hashes, "summary": summary}).encode()
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
            for name, data in files.items():
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        data = buffer.getvalue()
        chunks = []
        for offset in range(0, len(data), 31):
            block = data[offset:offset+31]
            name = f"part{len(chunks):03d}"
            (directory / name).write_bytes(block)
            chunks.append({"name": name, "bytes": len(block), "sha256": hashlib.sha256(block).hexdigest()})
        index = {"archive_sha256": hashlib.sha256(data).hexdigest(), "archive_bytes": len(data),
                 "chunks": chunks, "summary": summary, "files": hashes}
        (directory / "index.json").write_text(json.dumps(index))
        return index

    def test_arbitrary_chunk_boundaries(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.fixture(root)
            result = verify(root)
            self.assertTrue(result["accepted"])
            self.assertEqual(result["verified_files"], 1)

    def test_changed_byte_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            index = self.fixture(root)
            path = root / index["chunks"][0]["name"]
            data = bytearray(path.read_bytes())
            data[3] ^= 1
            path.write_bytes(data)
            with self.assertRaises(AssertionError):
                verify(root)

    def test_reordered_parts_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            index = self.fixture(root)
            index["chunks"].reverse()
            (root / "index.json").write_text(json.dumps(index))
            with self.assertRaises(AssertionError):
                verify(root)

    def test_manifest_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            index = self.fixture(root)
            index["files"]["example.txt"] = "0"*64
            (root / "index.json").write_text(json.dumps(index))
            with self.assertRaises(AssertionError):
                verify(root)

    def test_unsafe_part_path_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            index = self.fixture(root)
            index["chunks"][0]["name"] = "../outside"
            (root / "index.json").write_text(json.dumps(index))
            with self.assertRaises(AssertionError):
                verify(root)


if __name__ == "__main__":
    unittest.main()
