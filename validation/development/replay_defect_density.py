"""Recompute displaced-Si and frozen-density references from the raw archive."""
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
    with tempfile.TemporaryDirectory(prefix="defect-density-replay-") as tmp:
        root = Path(tmp)
        with tarfile.open(args.archive) as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root/"manifest.json").read_text())
        for name, expected in manifest["files"].items():
            path = root/name
            assert path.is_file() and path.stat().st_size == expected["bytes"], name
            with path.open("rb") as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"], name
        scripts = root/"build-serial"
        for build in ("build-serial", "build-mpi"):
            subprocess.run([sys.executable, str(scripts/"check_frozen_gamma_density.py"),
                            str(root/build/"frozen-density-after"), "--replay", "--check"], check=True)
        subprocess.run([sys.executable, str(scripts/"check_frozen_gamma_density.py"),
                        str(root/"build-serial/frozen-density-before"), "--replay"], check=True)
        before = json.loads((root/"build-serial/frozen-density-before/summary.json").read_text())
        assert before["GPW-TPSS"]["gamma_band_error_hartree"] > 0.004
        assert before["GAPW_XC-TPSS"]["gamma_band_error_hartree"] > 0.015
        subprocess.run([sys.executable, str(scripts/"check_silicon_translation.py"),
                        str(root/"build-mpi/silicon-translation-projector-tight"), "--replay"], check=True)
        for build in ("build-serial", "build-mpi"):
            log = (root/build/"frozen-density-regtests.log").read_text()
            assert "Status: OK" in log
            assert "WRONG RESULT" not in log and "RUNTIME FAIL" not in log
        print(json.dumps({"verified_files": len(manifest["files"]), "replay": "passed"}, indent=2))


if __name__ == "__main__":
    main()
