"""Archive completed stanene controls and the matching operator source snapshot."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--basis", choices=["DZVP", "TZVP"], required=True)
    args = parser.parse_args()
    root = args.source.resolve()
    tracked = ["src/qs_spectral_localizer.F", "src/qs_localizer_torus.F",
               "src/qs_finite_ao.F", "src/qs_gamma2kp.F", "src/qs_wannier90.F",
               "src/qs_moments.F", "src/spectral_localizer.F", "src/spectral_localizer_dbcsr.F",
               "src/spectral_localizer_sparse.F", "src/localizer_sparse_inertia.F",
               "src/localizer_sparse_pfaffian.F", "src/tacho_c_api.cpp", "src/topology_wilson.F",
               "src/input_cp2k_properties_dft.F", "src/input_cp2k_print_dft.F"]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *tracked], cwd=root, check=True)
    files = set(map(Path, tracked))
    for source in sorted(root.glob("build-*/stanene-localizer-validation/*/input.inp")):
        if f"BASIS_SET {args.basis}-" not in source.read_text():
            continue
        assert "PROGRAM ENDED AT" in source.with_name("output.out").read_text(), source
        for path in source.parent.iterdir():
            if path.suffix in {".out", ".inp", ".log", ".json", ".mmn", ".wilson", ".eig"}:
                files.add(path.relative_to(root))
    files.update(Path("build-serial")/name for name in (
        "check_stanene_localizer.py", "summarize_stanene_localizer.py"))
    files.update(Path("data")/name for name in ("BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS"))
    files.update(Path(build)/"CMakeCache.txt" for build in ("build-serial", "build-mpi"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "basis": args.basis, "files": {}, "executable_sha256": {},
                "notes": "Fresh 25 September 2026 two-dimensional XY calculations, not the earlier 3D-supercell Wilson data. "
                         "All localizer queries and full final Wilson overlaps are retained. No band flattening or occupied-space "
                         "projection was used for the localizer. Source snapshot and CMake caches document the local builds; "
                         "the embedded executable git string can predate the last incremental rebuild."}
    for name in ("build-serial/bin/cp2k.ssmp", "build-mpi/bin/cp2k.psmp"):
        with (root/name).open("rb") as stream:
            manifest["executable_sha256"][name] = hashlib.file_digest(stream, "sha256").hexdigest()
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in sorted(files):
            path = root/name
            assert path.is_file() and path.resolve().is_relative_to(root), name
            with path.open("rb") as stream:
                size = path.stat().st_size
                manifest["files"][str(name)] = {"bytes": size, "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}
                stream.seek(0)
                info = tarfile.TarInfo(name.as_posix())
                info.mode, info.size = 0o644, size
                archive.addfile(info, stream)
        data = (json.dumps(manifest, indent=2)+"\n").encode()
        info = tarfile.TarInfo("manifest.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    with args.archive.open("rb") as stream:
        print(json.dumps({"basis": args.basis, "files": len(files), "bytes": args.archive.stat().st_size,
                          "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
