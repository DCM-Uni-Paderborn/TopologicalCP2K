"""Retain TQC transport evidence and replay exact ordered records from the archive."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

from archive_stanene_convergence import digest, pack
from compare_tqc_transport import compare


COLLECTIONS = {
    "baseline-serial": ("build-serial/tqc-transport-baseline-regtests", 135),
    "baseline-mpi": ("build-mpi/tqc-transport-baseline-regtests", 135),
    "serial": ("build-serial/tqc-transport-regtests", 139),
    "mpi4": ("build-mpi/tqc-transport-regtests", 139),
    "mpi3": ("build-mpi/tqc-transport-three-regtests", 139),
    "wilson-serial": ("build-serial/tqc-transport-wilson-regtests", 47),
    "wilson-mpi": ("build-mpi/tqc-transport-wilson-regtests", 47),
}
PAIRS = [("baseline-serial", "serial"), ("baseline-mpi", "mpi4"),
         ("serial", "mpi4"), ("serial", "mpi3")]


def replay(path):
    index = json.loads((path.parent / "index.json").read_text())
    if digest(path) != index["archive"]["sha256"]:
        raise ValueError("Changed archive checksum")
    with tarfile.open(path, "r:gz") as archive:
        entries = {m.name: m for m in archive.getmembers()}
        if len(entries) != len(archive.getmembers()):
            raise ValueError("Duplicate archive member")
        manifest = json.loads(archive.extractfile("manifest.json").read())
        if set(entries) != set(manifest["files"]) | {"manifest.json"}:
            raise ValueError("Unlisted archive member")
        files = {}
        for name, expected in manifest["files"].items():
            data = archive.extractfile(name).read()
            if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                raise ValueError("Changed evidence: " + name)
            files[name] = data
    results = []
    for left, right in PAIRS:
        a = {Path(n).name: b for n, b in files.items()
             if n.startswith(left + "/") and n.endswith(".little_group")}
        b = {Path(n).name: data for n, data in files.items()
             if n.startswith(right + "/") and n.endswith(".little_group")}
        extra = {"helium-unitary.little_group"} if left.startswith("baseline") else set()
        if len(a) != (16 if extra else 17) or set(b) - set(a) != extra or set(a) - set(b):
            raise ValueError("Missing or extra TQC records")
        checked = [compare(a[n].decode(), b[n].decode()) for n in sorted(a)]
        results.append(dict(reference=left, candidate=right, files=len(checked),
                            records=sum(x["records"] for x in checked),
                            max_error=max(x["maximum_absolute_difference"] for x in checked)))
    return dict(payload_files=len(files), comparisons=results, accepted=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    archive = args.destination / "methods-results.tar.gz"
    if args.replay:
        print(json.dumps(replay(archive), indent=2))
        return
    root = args.root.resolve()
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root):
        raise ValueError("Native sources are not a clean checkpoint")
    files = {}
    for name in ("build-serial/tqc-transport-debug/all-settings.log", "build-mpi/tqc-transport-all-settings.log"):
        text = (root / name).read_text()
        if "530 Hall settings, 12720 cases passed." not in text or "12720 passed." not in text:
            raise ValueError("Incomplete all-setting algebra tests")
    probes = json.loads((root / "build-mpi/tqc-transport-probes/results.json").read_text())
    if {(x["ranks"], x["mode"]) for x in probes} != {(n, m) for n in (1, 3, 4) for m in ("valid", "incompatible", "invalid")}:
        raise ValueError("Missing collective control")
    if not all(x["passed"] and ((x["returncode"] != 0) == (x["mode"] == "invalid")) for x in probes):
        raise ValueError("Failed collective control")
    for label, (name, expected) in COLLECTIONS.items():
        log = root / (name + ".log")
        if f"correct: {expected} / {expected}" not in log.read_text() or "Status: OK" not in log.read_text():
            raise ValueError("Incomplete regtest collection: " + name)
        files[label + "/driver.log"] = log
        directories = list((root / name).glob("TEST-*/QS/regtest-*"))
        target = "regtest-property-wilson" if label.startswith("wilson") else "regtest-little-group"
        directories = [p for p in directories if p.name == target]
        if len(directories) != 1:
            raise ValueError("Ambiguous native run")
        for path in directories[0].iterdir():
            if path.is_file() and path.suffix in {".inp", ".inc", ".nnkp", ".out", ".little_group", ".toml", ".wilson"}:
                files[label + "/" + path.name] = path
    for name in ("src/qs_little_group.F", "src/qs_wannier90.F", "src/topology_symmetry.F",
                 "src/topology_symmetry_unittest.F", "docs/methods/properties/inversion_topology.md",
                 "build-serial/tqc-transport-pretty-final.log", "build-serial/tqc-transport-baseline-binaries.sha256",
                 "build-serial/tqc-transport-debug/all-settings.log",
                 "build-mpi/tqc-transport-all-settings.log"):
        files[name] = root / name
    for build in ("build-serial", "build-mpi"):
        for name in ("CMakeCache.txt", "tqc-transport-build.log"):
            files[build + "/" + name] = root / build / name
    for path in (root / "build-mpi/tqc-transport-probes").rglob("*"):
        if path.suffix in {".f90", ".json", ".txt", ".log"} and path.is_file():
            files[str(path.relative_to(root))] = path
    for name in ("archive_tqc_transport.py", "compare_tqc_transport.py", "archive_stanene_convergence.py"):
        files["methods/" + name] = Path(__file__).parent / name
    for path in (root / "build-serial").glob("tqc-transport-*.json"):
        files[str(path.relative_to(root))] = path
    metadata = dict(source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                    baseline_commit=subprocess.check_output(["git", "rev-parse", "HEAD^"], cwd=root, text=True).strip(),
                    regtest_checks={label: count for label, (_, count) in COLLECTIONS.items()},
                    native_binaries={}, native_libraries={})
    for build, version in (("build-serial", "ssmp"), ("build-mpi", "psmp")):
        metadata["native_binaries"][build] = digest(root / build / "bin" / ("cp2k." + version))
        metadata["native_libraries"][build] = digest(root / build / "src/libcp2k.2026.2.dylib")
    args.destination.mkdir(parents=True, exist_ok=True)
    metadata["archive"] = pack(archive, files, metadata)
    (args.destination / "index.json").write_text(json.dumps(metadata, indent=2) + "\n")
    result = replay(archive)
    (args.destination / "replay.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
