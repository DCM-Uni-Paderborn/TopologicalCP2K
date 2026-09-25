"""Create and verify a bounded evidence archive for the invariant-torus audit."""
import argparse
import hashlib
import json
import platform
import subprocess
import tarfile
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paper", type=Path)
    args = parser.parse_args()
    root = Path.cwd()
    output = args.paper / "validation/development"
    changed = subprocess.check_output(["git", "diff-tree", "--no-commit-id", "--name-only", "-r", "HEAD"], text=True).splitlines()
    files = {root / name for name in changed}
    files.update(root.glob("src/topology_*.F"))
    for build in ("build-serial", "build-mpi"):
        for pattern in ("fourier-*.log", "fourier-*.json", "fourier-*.jsonl", "fourier-*.patch"):
            files.update((root / build).glob(pattern))
        files.update((root / build / "fourier-character-negative").glob("*"))
        test, = (root / build / "fourier-character-regtests").glob("TEST-*")
        for suffix in ("*.inp", "*.inc", "*.nnkp", "*.out", "*.little_group"):
            files.update((test / "QS/regtest-little-group").glob(suffix))
    for name in ("ebr-character-full.dat", "ebr-character-full.json", "character-convention-reference.dat"):
        files.add(root / "build-serial" / name)
    helpers = ("audit_fourier_characters.py", "audit_character_class_export.py", "audit_ebr_characters.py",
               "build_fourier_instrumented.sh", "check_ebr_rejections.py", "prepare_ebr_characters.py",
               "prepare_wyckoff_characters.py", "verify_onsite_export.py", "compare_onsite_exports.py",
               "archive_fourier_validation.py")
    for name in helpers:
        files.add(root / "build-serial/compatibility-oracle" / name)
    environment = root / "build-serial/fourier-character-environment.json"
    commands = {"revision": ["git", "rev-parse", "HEAD"],
                "compiler": ["/opt/homebrew/bin/gfortran", "--version"],
                "mpi": ["mpiexec", "--version"]}
    details = {name: subprocess.check_output(command, text=True).strip() for name, command in commands.items()}
    details.update(platform=platform.platform(), omp_threads=2, openblas_threads=1, omp_stacksize="64M",
                   source_files={name: sha((root / name).read_bytes()) for name in changed},
                   binaries={str(p.relative_to(root)): sha(p.read_bytes()) for p in
                             (root / "build-serial/bin/topology_band_unittest.ssmp",
                              root / "build-mpi/bin/topology_band_unittest.psmp")})
    environment.write_text(json.dumps(details, indent=2) + "\n")
    files.add(environment)
    files = sorted(files)
    assert all(p.is_file() and not p.is_symlink() for p in files)
    hashes = {str(p.relative_to(root)): sha(p.read_bytes()) for p in files}
    manifest = output / "fourier-character-files.sha256"
    manifest.write_text("".join(f"{digest}  {name}\n" for name, digest in hashes.items()))
    archive = output / "fourier-character-records.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for path in files:
            tar.add(path, arcname=str(path.relative_to(root)), recursive=False)
    with tarfile.open(archive) as tar:
        assert {m.name for m in tar} == set(hashes)
        for member in tar:
            assert member.isfile()
            assert sha(tar.extractfile(member).read()) == hashes[member.name]
    result = dict(files=len(files), uncompressed_bytes=sum(p.stat().st_size for p in files),
                  archive_bytes=archive.stat().st_size, archive_sha256=sha(archive.read_bytes()),
                  all_member_hashes_verified=True, revision=details["revision"])
    (output / "fourier-character-archive-check.json").write_text(json.dumps(result, indent=2) + "\n")
    for name in helpers:
        (output / name).write_bytes((root / "build-serial/compatibility-oracle" / name).read_bytes())
    for name in ("fourier-character-comparison.json", "fourier-ebr-comparison.json", "fourier-character-classes.json",
                 "fourier-character-export-comparison.json"):
        (output / name).write_bytes((root / "build-serial" / name).read_bytes())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
