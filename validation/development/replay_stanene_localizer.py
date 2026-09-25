"""Verify a stanene archive and replay raw Wilson matrices and solver comparisons."""

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
    output = args.output.resolve()
    with tempfile.TemporaryDirectory(prefix="stanene-localizer-replay-") as tmp:
        root = Path(tmp)
        with tarfile.open(args.archive) as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root/"manifest.json").read_text())
        for name, expected in manifest["files"].items():
            path = root/name
            assert path.is_file() and path.stat().st_size == expected["bytes"], name
            with path.open("rb") as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"], name
        subprocess.run([sys.executable, str(root/"build-serial/summarize_stanene_localizer.py"),
                        str(root), str(output), "--basis", manifest["basis"]], check=True)
        print(json.dumps({"verified_files": len(manifest["files"]), "replay": "passed"}, indent=2))


if __name__ == "__main__":
    main()
