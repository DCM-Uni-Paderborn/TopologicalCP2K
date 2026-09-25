"""Archive the SOC image-sign diagnosis, corrected runs and full-band references."""

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
    args = parser.parse_args()
    root = args.source.resolve()
    tracked = ["src/qs_localizer_torus.F", "src/qs_spectral_localizer.F", "src/spectral_localizer.F",
               "src/spectral_localizer_dbcsr.F", "src/spectral_localizer_sparse.F", "src/core_ppnl.F",
               "src/soc_pseudopotential_methods.F", "src/qs_wannier90.F", "src/kpoint_methods.F",
               "src/post_scf_bandstructure_utils.F"]
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *tracked], cwd=root, check=True)
    files = set(map(Path, tracked))
    for source in sorted(root.glob("build-*/stanene-bloch-localizer/*/input.inp")):
        assert "PROGRAM ENDED AT" in source.with_name("output.out").read_text(), source
        for path in source.parent.iterdir():
            if path.suffix in {".out", ".inp", ".log", ".json", ".mmn", ".eig", ".topology", ".nnkp", ".bin"}:
                files.add(path.relative_to(root))
    files.update(Path("build-serial")/name for name in (
        "check_stanene_localizer.py", "check_stanene_soc_correction.py", "check_bloch_localizer.py", "gaussian_states.py",
        "summarize_stanene_localizer.py"))
    files.update(root.joinpath("tests/QS").glob("regtest-spectral-localizer*/stanene-*"))
    files = {p.relative_to(root) if p.is_absolute() else p for p in files}
    files.update(Path("data")/name for name in ("BASIS_MOLOPT_UZH", "GTH_SOC_POTENTIALS"))
    files.update(Path(build)/"CMakeCache.txt" for build in ("build-serial", "build-mpi"))
    manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                "files": {}, "notes": "direct3 is the before-fix negative control; corrected* use the source snapshot. "
                "Matrix dumps are little-endian int32 dimension, then complex128 L and B in column-major order. "
                "Dumps used temporary source-rank stream output immediately before localizer_matrix_evaluate; "
                "this instrumentation was removed from the final source. No flattening or polar-coordinate replacement "
                "was used. Runtime commit banners can predate the incremental rebuild."}
    with args.archive.open("xb") as output, tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in sorted(files):
            path = root/name
            assert path.is_file() and path.resolve().is_relative_to(root)
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
        print(json.dumps({"files": len(files), "bytes": args.archive.stat().st_size,
                          "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
