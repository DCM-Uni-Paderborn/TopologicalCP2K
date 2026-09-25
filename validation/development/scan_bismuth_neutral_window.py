"""Compare finite Bi localizers at energies inside each case's own neutral gap.

This is an independent postprocessing scan, not a replacement for native
CP2K comparisons. Bare finite flakes are not bulk topological predictions.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import block_diag, eigh

import analyze_bismuth_flakes as analysis


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def window(snapshot, fractions, offsets):
    energies = np.asarray(snapshot.energies)
    positions = np.asarray(snapshot.positions)
    occupied = 5 * len(positions)
    if (energies.ndim != 2 or energies.shape[0] != 1 or occupied < 2 or
            occupied % 2 or energies.shape[1] <= occupied or
            positions.ndim != 2 or positions.shape[1] != 3):
        raise ValueError("A neutral finite Bi spectrum with conduction states is required")
    if (not np.isfinite(energies).all() or not np.isfinite(positions).all() or
            np.any(np.diff(energies[0]) < -1e-12)):
        raise ValueError("Nonfinite geometry or unordered spectrum")
    if (not fractions or not offsets or
            not np.isfinite([*fractions, *offsets]).all() or
            max(abs(f) for f in fractions) >= .5):
        raise ValueError("Energy fractions must be strictly inside the neutral gap")
    lower, upper = energies[0, occupied - 1:occupied + 1]
    gap = upper - lower
    if gap <= 1e-9:
        raise ValueError("The neutral spectral gap is unresolved")
    center = np.mean(positions, axis=0)
    half_width = np.ptp(positions[:, 0]) / 2
    if half_width <= 0:
        raise ValueError("The finite patch has no Cartesian x extent")
    selected = (upper + lower) / 2 + gap * np.asarray(fractions)
    points = [center + [f * half_width, 0., 0.] for f in offsets]
    return selected, points, dict(
        occupied_spinors=occupied, gap_hartree=float(gap),
        midpoint_hartree=float((upper + lower) / 2),
        x_half_width_bohr=float(half_width), energy_gap_fractions=fractions,
        x_half_width_fractions=offsets)


def frontier_localization(snapshot, overlap):
    """Trace normalized Lowdin atom projectors over complete frontier eigenspaces."""
    eigenvalues = snapshot.energies[0]
    occupied = 5 * len(snapshot.positions)
    gap = eigenvalues[occupied] - eigenvalues[occupied - 1]
    tolerance = min(1e-8, gap / 10)
    if tolerance <= 0:
        raise ValueError("Frontier states are not spectrally separated")
    w, v = eigh(overlap)
    if w[0] <= 1e-10 * w[-1]:
        raise ValueError("Unresolved AO metric")
    root = (v * np.sqrt(w)) @ v.conj().T
    frame = block_diag(root, root) @ snapshot.coefficients[0].reshape(2 * len(w), -1)
    if np.max(abs(frame.conj().T @ frame - np.eye(frame.shape[1]))) > 1e-7:
        raise ValueError("Nonorthonormal frontier eigenframe")
    distances = np.linalg.norm(snapshot.positions[:, None] - snapshot.positions[None, :], axis=2)
    # Bi nearest bonds are 3.043 Angstrom; second neighbors are 4.33 Angstrom.
    cutoff = 3.3 / analysis.BOHR_ANGSTROM
    coordination = np.sum((distances > 1e-8) & (distances < cutoff), axis=1)
    if np.max(coordination) > 3:
        raise ValueError("Geometry is not the assumed three-coordinated Bi patch")
    edge = coordination < 3
    boundaries = np.r_[0, np.cumsum(snapshot.atom_sizes)]
    bands = []
    for name, index in (("occupied", occupied - 1), ("unoccupied", occupied)):
        group = np.flatnonzero(abs(eigenvalues - eigenvalues[index]) <= tolerance)
        density = np.sum(abs(frame[:, group])**2, axis=1).reshape(2, -1).sum(axis=0) / len(group)
        weights = [float(np.sum(density[start:end])) for start, end in zip(boundaries[:-1], boundaries[1:])]
        if abs(sum(weights) - 1) > 1e-7:
            raise ValueError("Non-normalized frontier weights")
        bands.append(dict(side=name, bands_one_based=list(map(int, group + 1)),
                          degeneracy_tolerance_hartree=float(tolerance),
                          atom_weights=weights, edge_weight=float(np.sum(np.asarray(weights)[edge]))))
    return dict(definition="Lowdin atom weights averaged over complete frontier eigenspaces; basis dependent",
                bond_cutoff_angstrom=3.3, coordination=list(map(int, coordination)),
                edge_atoms_one_based=list(map(int, np.flatnonzero(edge) + 1)),
                edge_atom_fraction=float(np.mean(edge)), frontier=bands)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--energy-fraction", nargs="+", type=float, default=[-.25, 0., .25])
    parser.add_argument("--position-fraction", nargs="+", type=float, default=[0., .75, 1.5])
    parser.add_argument("--kappa", nargs="+", type=float, default=[.0001, .0003, .001, .003, .01])
    args = parser.parse_args()
    root, work = args.root.resolve(), args.work.resolve()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not args.kappa or not np.isfinite(args.kappa).all() or min(args.kappa) <= 0:
        raise ValueError("Positive finite scales are required")
    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import read_snapshot
    from check_bloch_localizer import aii_skew, skew_sign

    run = json.loads((work / "run.json").read_text())
    original = json.loads((work / "independent-localizers.json").read_text())
    if not run["completed"] or not run["scf_converged"] or run["returncode"] != 0:
        raise ValueError("Only a completed, converged calculation can be scanned")
    if run["options"]["size"] < 1 or run["options"]["mode"] == "wilson":
        raise ValueError("Only finite Bi patches are supported")
    snapshot_hash = digest(work / "bismuth.topology")
    method_hash = digest(Path(analysis.__file__))
    integral_hash = digest(root / "build-serial/gaussian_states.py")
    if (snapshot_hash != run["files"]["bismuth.topology"] or
            snapshot_hash != original["snapshot_sha256"] or
            method_hash != original["source_sha256"] or
            integral_hash != original["integral_sha256"]):
        raise ValueError("Snapshot or independent analysis provenance has changed")
    with np.load(work / "independent-ao-moments.npz", allow_pickle=False) as archive:
        for key, expected in (("snapshot_sha256", snapshot_hash), ("method_sha256", method_hash),
                              ("integral_sha256", integral_hash)):
            if str(archive[key]) != expected:
                raise ValueError("Incompatible AO moment cache")
        integrated = archive["moments"]
    started = time.monotonic()
    snapshot = read_snapshot(work / "bismuth.topology")
    energies, points, spectrum = window(snapshot, args.energy_fraction, args.position_fraction)
    h, xyz, diagnostics = analysis.operators(snapshot, integrated)
    frontier = frontier_localization(snapshot, integrated[0])
    rows = analysis.queries(h, xyz, energies, args.kappa, points, aii_skew, skew_sign)
    temperature = run["options"]["temperature"]
    # Same CODATA-2006 constants as CP2K's physcon.F.
    thermal = temperature * 1.3806504e-23 / (1.602176487e-19 * analysis.HARTREE_EV)
    spectrum.update(temperature_kelvin=temperature, thermal_energy_hartree=thermal,
                    gap_over_kbt=spectrum["gap_hartree"] / thermal if thermal else None)
    files = ("run.json", "bismuth.topology", "independent-localizers.json", "independent-ao-moments.npz")
    report = dict(case=work.name, comparison_kind="independent-neutral-window-scan",
                  native_comparison_performed=False, spectrum=spectrum, diagnostics=diagnostics,
                  frontier_localization=frontier,
                  run_provenance=run["provenance"], options=run["options"],
                  inputs={f: digest(work / f) for f in files},
                  methods={"scan": digest(Path(__file__)), "operators": method_hash,
                           "gaussian_integrals": integral_hash,
                           "pfaffian_reference": digest(root / "build-serial/check_bloch_localizer.py")},
                  native_arguments=dict(energy=list(map(float, energies)), kappa=args.kappa,
                      offset=[f * spectrum["x_half_width_bohr"] * analysis.BOHR_ANGSTROM
                              for f in args.position_fraction]),
                  queries=rows, elapsed_seconds=time.monotonic() - started)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(case=work.name, spectrum=spectrum,
                          nontrivial=sum(q["z2"] == 1 for q in rows), queries=len(rows)), indent=2))


if __name__ == "__main__":
    main()
