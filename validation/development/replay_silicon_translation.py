"""Recompute Si band, Gamma-supercell and screening checks from archived raw output."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="silicon-translation-replay-") as tmp:
        root = Path(tmp)
        with tarfile.open(args.archive) as archive:
            archive.extractall(root, filter="data")
        manifest = json.loads((root/"manifest.json").read_text())
        for name, expected in manifest["files"].items():
            path = root/name
            assert path.is_file() and path.stat().st_size == expected["bytes"], name
            with path.open("rb") as source:
                assert hashlib.file_digest(source, "sha256").hexdigest() == expected["sha256"], name
        script = root/"build-serial/check_silicon_translation.py"
        serial = root/"build-serial/silicon-translation-projector-tight"
        mpi = root/"build-mpi/silicon-translation-projector-tight"
        subprocess.run([sys.executable, str(script), str(serial), "--case", "primitive", "--replay"], check=True)
        subprocess.run([sys.executable, str(script), str(mpi), "--pristine-only", "--replay"], check=True)
        sys.path.insert(0, str(script.parent))
        from check_silicon_translation import reference
        a = reference.values((serial/"primitive/output.out").read_text(), "Gap [hartree]").ravel()
        b = reference.values((mpi/"primitive/output.out").read_text(), "Gap [hartree]").ravel()
        assert np.max(abs(a-b)) < 1e-9
        control_errors = {}
        for case in ("silicon-translation-dense", "silicon-translation-tight"):
            directory = root/"build-serial"/case/"primitive"
            h = np.loadtxt(directory/"invalid-operator.dat")
            defect = float(np.max(abs(h-h.T))/np.max(abs(h)))
            assert 1e-8 < defect < 1e-7, defect
            control_errors[case] = defect
        orbital = (root/"build-serial/silicon-translation-orbital-tight/primitive/output.out").read_text()
        assert "4.1143E-08" in orbital and "[ABORT]" in orbital
        for build, total in (("build-serial", 30), ("build-mpi", 54)):
            log = (root/build/"silicon-screening-regtests.log").read_text()
            assert f"Summary: correct: {total} / {total}" in log
        print(json.dumps({"verified_files": len(manifest["files"]),
                          "serial_mpi_gap_difference": float(np.max(abs(a-b))),
                          "failed_control_relative_defects": control_errors}, indent=2))


if __name__ == "__main__":
    main()
