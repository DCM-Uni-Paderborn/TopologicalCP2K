"""Verify the archive and recompute manuscript comparisons from raw CP2K output."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tarfile

import numpy as np


def analyze_output(text):
    assert "PROGRAM ENDED AT" in text
    def values(pattern):
        return np.array([float(v) for v in re.findall(pattern, text)])
    gaps = values(r"QUADRATIC_PSEUDOSPECTRUM\| Gap \[hartree\]:\s+(\S+)")
    energy = values(r"QUADRATIC_PSEUDOSPECTRUM\| Energy residual \[hartree\]:\s+(\S+)")
    ritz = values(r"QUADRATIC_PSEUDOSPECTRUM\| Eigenpair residual \[hartree\^2\]:\s+(\S+)")
    chord = values(r"QUADRATIC_PSEUDOSPECTRUM\| Chord residual axis \d \[bohr\]:\s+(\S+)")
    assert len(gaps) == len(energy) == len(ritz) == 2
    assert np.isfinite(np.concatenate([gaps, energy, ritz, chord])).all()
    assert max(ritz) < 1e-10
    chord = chord.reshape(2, -1)
    decomposition = max(abs(gaps**2 - energy**2 - 0.1**2 * np.sum(chord**2, axis=1)))
    assert decomposition < 1e-13
    total, = values(r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+(\S+)")
    dimension, = re.findall(r"QUADRATIC_PSEUDOSPECTRUM\| AO/spinor dimension:\s+(\d+)", text)
    return dict(gaps=gaps.tolist(), energy_residuals=energy.tolist(), ritz_residuals=ritz.tolist(),
                chord_residuals=chord.tolist(), decomposition_error=float(decomposition),
                scf_energy=float(total), dimension=int(dimension))


def difference(a, b, field="gaps"):
    return float(np.max(np.abs(np.asarray(a[field]) - b[field])))


def summarize(path):
    with tarfile.open(path) as archive:
        contents = {}
        for member in archive.getmembers():
            assert member.isfile() and member.name not in contents
            assert not Path(member.name).is_absolute() and ".." not in Path(member.name).parts
            contents[member.name] = archive.extractfile(member).read()
    manifest = json.loads(contents.pop("manifest.json"))
    assert set(manifest["files"]) == set(contents)
    for name, raw in contents.items():
        assert len(raw) == manifest["files"][name]["bytes"]
        assert hashlib.sha256(raw).hexdigest() == manifest["files"][name]["sha256"]
    def output(directory, test):
        return analyze_output(contents[f"{directory}/{test}/result.out"].decode())
    tests = ["periodic-chain", "periodic-gapw", "periodic-helium", "periodic-reduced",
             "periodic-soc", "periodic-supercell", "periodic-wrapped"]
    serial = {t: output("build-serial/periodic-quadratic-smokes", t) for t in tests}
    mpi = {t: output("build-mpi/periodic-quadratic-smokes", t) for t in tests}
    gaps = {t: difference(serial[t], mpi[t]) for t in tests}
    energy = max(difference(serial[t], mpi[t], "energy_residuals") for t in tests)
    assert max(gaps.values()) < 6e-15 and energy < 2e-13
    variant_dir = "build-serial/periodic-quadratic-variants-checked"
    variants = {t: output(variant_dir, t) for t in ["chain-x", "chain-z", "gapw-yz", "macdonald",
                "general", "explicit-gamma-supercell", "skew", "skew-rotated"]}
    symmetry = max(difference(build["periodic-helium"], build["periodic-reduced"]) for build in [serial, mpi])
    wrapped = max(difference(build["periodic-helium"], build["periodic-wrapped"]) for build in [serial, mpi])
    schemes = max(difference(variants[t], serial["periodic-helium"]) for t in ["macdonald", "general"])
    rotation = difference(variants["skew"], variants["skew-rotated"])
    axes = max(difference(variants["chain-x"], serial["periodic-chain"]),
               difference(variants["chain-z"], serial["periodic-chain"]),
               difference(variants["gapw-yz"], serial["periodic-gapw"]))
    assert max(symmetry, wrapped, schemes, rotation, axes) < 3e-15
    supercell = difference(serial["periodic-helium"], serial["periodic-supercell"])
    explicit = difference(serial["periodic-supercell"], variants["explicit-gamma-supercell"])
    scf_supercell = abs(serial["periodic-helium"]["scf_energy"] - serial["periodic-supercell"]["scf_energy"] / 4)
    assert supercell < 1e-9 and explicit < 2.2e-12 and scf_supercell < 1e-9
    rejected = {}
    for name, message in {
        "guard-open-grid": "Quadratic MP_GRID must be one in open directions.",
        "guard-finite-grid": "Quadratic MP_GRID requires FORMULATION PERIODIC.",
        "guard-torus-size": "Quadratic full AO torus exceeds MAX_AO.",
        "guard-periodicity": "Quadratic torus needs matching, nonempty CELL/POISSON PERIODIC.",
    }.items():
        text = contents[f"{variant_dir}/{name}/result.out"].decode()
        assert "[ABORT]" in text and message in text and "PROGRAM ENDED AT" not in text
        rejected[name] = message
    regressions = {}
    for build, expected in [("build-serial", 52), ("build-mpi", 87)]:
        text = contents[f"{build}/periodic-quadratic-regtests.log"].decode()
        assert "Status: OK" in text and f"Summary: correct: {expected} / {expected}" in text
        regressions[build] = expected
    assert "Status: OK" in contents["build-serial/periodic-quadratic-pretty.log"].decode()
    records = list(serial.values()) + list(mpi.values()) + list(variants.values())
    return dict(source_commit=manifest["source_commit"], verified_files=len(contents),
                archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                serial_dense=serial, two_rank_iterative=mpi, serial_variants=variants,
                dense_iterative_gap_differences=gaps, max_energy_residual_difference=energy,
                max_iterative_ritz_residual=max(max(r["ritz_residuals"]) for r in mpi.values()),
                max_decomposition_error=max(r["decomposition_error"] for r in records),
                symmetry_difference=symmetry, query_wrapping_difference=wrapped,
                scheme_difference=schemes, rotation_difference=rotation, axis_permutation_difference=axes,
                supercell_gap_difference=supercell, supercell_energy_per_cell_difference=scf_supercell,
                explicit_implicit_gamma_difference=explicit, rejected_inputs=rejected, regressions=regressions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = summarize(args.archive)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      {"serial_dense", "two_rank_iterative", "serial_variants"}}, indent=2))


if __name__ == "__main__":
    main()
