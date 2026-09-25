"""Verify an extracted evidence archive and replay its native Fourier/EBR fixture."""
import argparse
import hashlib
import json
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    binary = args.binary.resolve()
    hashes = dict(line.split("  ", 1)[::-1] for line in
                  (evidence / "fourier-character-files.sha256").read_text().splitlines())
    with tempfile.TemporaryDirectory(prefix="fourier-replay-", dir=args.output.parent) as temporary:
        work = Path(temporary)
        with tarfile.open(evidence / "fourier-character-records.tar.gz") as tar:
            tar.extractall(work, filter="data")
        actual = {str(p.relative_to(work)) for p in work.rglob("*") if p.is_file()}
        assert actual == set(hashes)
        for name, digest in hashes.items():
            assert hashlib.sha256((work / name).read_bytes()).hexdigest() == digest
        result = subprocess.run([str(binary), "--ebr-reference=build-serial/ebr-character-full.dat"], cwd=work,
                                env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "2", "OMP_STACKSIZE": "64M"},
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300)
        log = work / "replay.log"
        log.write_text(result.stdout)
        assert result.returncode == 0 and "Site-induced atomic band tests passed." in result.stdout
        check = subprocess.check_output(["python3", str(evidence / "audit_fourier_characters.py"), str(log)], text=True)
        summary = json.loads(check)["results"][0]["ranks"][0]
        args.output.write_text(json.dumps(dict(files_verified=len(hashes), all_hashes_match=True,
                                              native_exit=result.returncode, audit=summary), indent=2) + "\n")
    print(args.output.read_text())


if __name__ == "__main__":
    main()
