"""Verify and replay the finite Gaussian translation archive with NumPy."""
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
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="gaussian-translation-replay-") as directory:
        root = Path(directory).resolve()
        with tarfile.open(args.archive, "r:gz") as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root / "manifest.json").read_text())
        for name, item in manifest["files"].items():
            path = root / name
            assert path.resolve().is_relative_to(root)
            data = path.read_bytes()
            assert len(data) == item["bytes"]
            assert hashlib.sha256(data).hexdigest() == item["sha256"], name
        scripts = root / "build-serial/gaussian-translation-validation"
        results = {}
        for build in ("build-serial", "build-mpi"):
            output = root / build / "gaussian-translation-verified"
            subprocess.run([sys.executable, str(scripts / "check_gaussian_translation.py"),
                            "--root", str(root), "--binary", str(root / build / "bin/cp2k"),
                            "--output", str(output), "--replay"], check=True)
            results[build] = json.loads((output / "summary.json").read_text())
        output = root / "build-mpi/gaussian-translation-extensions"
        subprocess.run([sys.executable, str(scripts / "check_translation_extensions.py"),
                        "--root", str(root), "--output", str(output), "--replay"], check=True)
        results["extensions"] = json.loads((output / "summary.json").read_text())
        for build in ("build-serial", "build-mpi"):
            log = (root / build / "gaussian-translation-final-regtests.log").read_text()
            assert "Number of FAILED  tests 0" in log and "Number of WRONG   tests 0" in log
            results[build + "-regtests"] = next(line for line in log.splitlines() if line.startswith("Summary:"))
        results["source_commit"] = manifest["source_commit"]
        results["files_verified"] = len(manifest["files"])
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
