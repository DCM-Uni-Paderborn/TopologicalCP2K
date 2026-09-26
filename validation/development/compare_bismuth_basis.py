"""Controlled basis or temperature comparisons of complete finite SOC exports.

These diagnostics measure basis sensitivity, not agreement with a reference
material. In particular, the two SCF Hamiltonians need not coincide.
"""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.linalg import block_diag, eigh, eigvalsh, solve, svdvals

import analyze_bismuth_flakes as analysis


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cross_moments(left, right, primitive_integrals):
    """Integrate both AO bases in the same real-space frame, without padding MOs."""
    if (left.positions.shape != right.positions.shape or
            not np.allclose(left.positions, right.positions, atol=1e-11, rtol=0) or
            not np.allclose(left.cell, right.cell, atol=1e-11, rtol=0) or
            not np.array_equal(left.kinds, right.kinds) or
            np.any(left.periodic) or np.any(right.periodic) or
            np.max(abs(left.kpoints)) > 1e-14 or np.max(abs(right.kpoints)) > 1e-14):
        raise ValueError("Basis comparisons require the same isolated geometry")
    nl = int(sum(left.atom_sizes))
    combined = SimpleNamespace(periodic=np.zeros(3, int), kpoints=np.zeros((1, 3)),
        atom_sizes=[*left.atom_sizes, *right.atom_sizes],
        shells=[*left.shells, *(replace(s, indices=s.indices + nl) for s in right.shells)])
    integrated = analysis.moments(combined, primitive_integrals)
    return integrated[:, :nl, :nl], integrated[:, nl:, nl:], integrated[:, :nl, nl:]


def inverse_root(metric):
    if not np.isfinite(metric).all() or np.max(abs(metric - metric.conj().T)) > 1e-9:
        raise ValueError("Invalid AO metric")
    w, v = eigh(metric)
    if w[0] <= 1e-10 * w[-1]:
        raise ValueError("Unresolved AO metric")
    return (v / np.sqrt(w)) @ v.conj().T


def subspace_overlap(left, right, sl, sr, cross):
    """Gauge-invariant principal angles of physical, individually normalized states."""
    if (left.shape[0] != sl.shape[0] or right.shape[0] != sr.shape[0] or
            cross.shape != (sl.shape[0], sr.shape[0]) or
            min(left.shape[1], right.shape[1]) < 1 or
            not all(np.isfinite(a).all() for a in (left, right, sl, sr, cross))):
        raise ValueError("Invalid cross-basis subspace dimensions or values")
    residuals = [float(np.max(abs(c.conj().T @ s @ c - np.eye(c.shape[1]))))
                 for c, s in ((left, sl), (right, sr))]
    if max(residuals) > 1e-7:
        raise ValueError("States are not orthonormal in their physical AO metrics")
    singular = svdvals(left.conj().T @ cross @ right)
    if singular[0] > 1 + 1e-7:
        raise ValueError("Cross overlap is not a contraction between physical subspaces")
    squared = np.minimum(singular**2, 1.)
    ranks = [left.shape[1], right.shape[1]]
    return dict(ranks=ranks, metric_residuals=residuals,
        overlap_singular_values=list(map(float, singular)),
        projector_frobenius_squared=float(sum(ranks) - 2 * np.sum(squared)),
        projector_spectral_distance=1. if ranks[0] != ranks[1] else float(np.sqrt(1-squared[-1])),
        mean_retained_left_weight=float(np.sum(squared) / ranks[0]),
        mean_retained_right_weight=float(np.sum(squared) / ranks[1]))


def controlled_options(runs, control):
    """Permit only the selected physical parameter and different restart guesses."""
    if control not in ("basis", "temperature"):
        raise ValueError("Unknown comparison control")
    keys = ["source_commit", "executable_sha256", "library_sha256", "sources", "runtime"]
    if control == "temperature":
        keys.append("runner_sha256")
    for key in keys:
        if runs[0]["provenance"][key] != runs[1]["provenance"][key]:
            raise ValueError("Uncontrolled change in native provenance: " + key)
    ignored = {"root", "output", "restart", control}
    for key in set(runs[0]["options"]) | set(runs[1]["options"]):
        if key not in ignored and runs[0]["options"].get(key) != runs[1]["options"].get(key):
            raise ValueError("Uncontrolled change in calculation options: " + key)


def temperature_inputs(texts, temperatures):
    """Audit the actual generated inputs, not only their command-line options."""
    normalized = []
    for text, expected in zip(texts, temperatures):
        lines, count = [], 0
        for line in text.splitlines():
            tokens = line.split()
            if tokens and tokens[0].upper() == "ELECTRONIC_TEMPERATURE":
                if len(tokens) != 2 or float(tokens[1]) != expected:
                    raise ValueError("Input temperature disagrees with run options")
                count += 1
                lines.append("ELECTRONIC_TEMPERATURE <controlled>")
            else:
                lines.append(line)
        if count != 1:
            raise ValueError("Exactly one electronic temperature is required")
        normalized.append(lines)
    if len(normalized) != 2 or normalized[0] != normalized[1]:
        raise ValueError("Inputs differ beyond electronic temperature")


def centered_hamiltonian_change(left, right, metric, midpoints):
    """Measure the physical operator difference after aligning neutral midgaps.

    AO matrices are covariant: subtracting c*S, not c*I, removes an energy
    shift. Congruence with S**(-1/2) gives the physical operator norm.
    """
    if (left.shape != metric.shape or right.shape != metric.shape or
            np.shape(midpoints) != (2,) or
            not all(np.isfinite(a).all() for a in (left, right, midpoints)) or
            max(np.max(abs(a-a.conj().T)) for a in (left, right)) > 1e-9):
        raise ValueError("Invalid Hermitian Hamiltonians or energy centers")
    inv = inverse_root(metric)
    delta = inv @ (right-left) @ inv
    shift = float(midpoints[1]-midpoints[0])
    centered = inv @ (right-left-shift*metric) @ inv
    return dict(energy_zero_shift_hartree=shift,
        absolute_operator_change_hartree=float(np.max(abs(eigvalsh(delta)))),
        midpoint_aligned_operator_change_hartree=float(np.max(abs(eigvalsh(centered)))))


def load_case(work, read_snapshot, integral_hash):
    run = json.loads((work / "run.json").read_text())
    record = json.loads((work / "independent-localizers.json").read_text())
    if not run["completed"] or not run["scf_converged"] or run["returncode"] != 0:
        raise ValueError("Incomplete native calculation")
    for name, expected in run["files"].items():
        if digest(work / name) != expected:
            raise ValueError("Changed native evidence: " + name)
    snapshot_hash = digest(work / "bismuth.topology")
    if (record["snapshot_sha256"] != snapshot_hash or
            record["source_sha256"] != digest(Path(analysis.__file__)) or
            record["integral_sha256"] != integral_hash):
        raise ValueError("Changed independent-analysis provenance")
    with np.load(work / "independent-ao-moments.npz", allow_pickle=False) as cache:
        if (str(cache["snapshot_sha256"]) != snapshot_hash or
                str(cache["method_sha256"]) != record["source_sha256"] or
                str(cache["integral_sha256"]) != integral_hash):
            raise ValueError("Changed integral cache provenance")
        moments = cache["moments"]
    snapshot = read_snapshot(work / "bismuth.topology")
    analysis.operators(snapshot, moments)
    return run, snapshot, moments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--control", choices=["basis", "temperature"], default="basis")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root = args.root.resolve()
    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import primitive_integrals, read_snapshot

    integral_hash = digest(root / "build-serial/gaussian_states.py")
    runs, snapshots, own = zip(*(load_case(p.resolve(), read_snapshot, integral_hash)
                               for p in (args.left, args.right)))
    controlled_options(runs, args.control)
    if args.control == "temperature":
        temperature_inputs([(p / "input.inp").read_text() for p in (args.left, args.right)],
                           [run["options"]["temperature"] for run in runs])
    left, right = snapshots
    ml, mr, cross = cross_moments(left, right, primitive_integrals)
    errors = [float(np.max(abs(a-b))) for a, b in zip(own, (ml, mr))]
    if max(errors) > 1e-10:
        raise ValueError("Combined integration differs from independently cached integrals")
    if args.control == "temperature" and (ml.shape != mr.shape or
            max(np.max(abs(ml-mr)), np.max(abs(ml-cross))) > 1e-10):
        raise ValueError("Temperature comparison requires identical physical AO bases and positions")
    sl, sr = ml[0], mr[0]
    il, ir = inverse_root(sl), inverse_root(sr)
    embedding = solve(sr, cross[0].conj().T, assume_a="pos")
    embedding_error = il @ (sl - embedding.conj().T @ sr @ embedding) @ il
    nesting = dict(singular_values=list(map(float, svdvals(il @ cross[0] @ ir))),
        left_metric_projection_error=float(np.linalg.norm(embedding_error, 2)),
        position_projection_errors_bohr=[float(np.linalg.norm(
            il @ (ml[d] - embedding.conj().T @ mr[d] @ embedding) @ il, 2)) for d in (1, 2, 3)])
    cl, cr = [s.coefficients[0].reshape(2 * len(m[0]), -1) for s, m in zip(snapshots, own)]
    spin_sl, spin_sr, spin_cross = [block_diag(a, a) for a in (sl, sr, cross[0])]
    occupied = 5 * len(left.positions)
    if any(s.energies[0, occupied] - s.energies[0, occupied-1] <= 1e-9 for s in snapshots):
        raise ValueError("Unresolved neutral spectral gap")
    groups = [("occupied", np.arange(occupied), np.arange(occupied))]
    for name, index in (("frontier_occupied", occupied-1), ("frontier_unoccupied", occupied)):
        selected = [np.flatnonzero(abs(s.energies[0] - s.energies[0, index]) <= min(
            1e-8, (s.energies[0, occupied] - s.energies[0, occupied-1]) / 10)) for s in snapshots]
        groups.append((name, *selected))
    subspaces = []
    for name, a, b in groups:
        subspaces.append(dict(name=name, left_bands_one_based=list(map(int, a+1)),
            right_bands_one_based=list(map(int, b+1)),
            **subspace_overlap(cl[:, a], cr[:, b], spin_sl, spin_sr, spin_cross)))
    left_frame = spin_sl @ cl
    right_frame = spin_sr @ cr
    hl = (left_frame * left.energies[0]) @ left_frame.conj().T
    hr = (right_frame * right.energies[0]) @ right_frame.conj().T
    emb, inv = block_diag(embedding, embedding), block_diag(il, il)
    projected_h_change = float(np.linalg.norm(inv @ (hl - emb.conj().T @ hr @ emb) @ inv, 2))
    report = dict(cases=[args.left.name, args.right.name],
        interpretation="Cross-basis sensitivity at identical geometry; not a material convergence certificate",
        methods=dict(comparison=digest(Path(__file__)), operators=digest(Path(analysis.__file__)),
                     gaussian_integrals=integral_hash),
        inputs=[{name: digest(p / name) for name in ("run.json", "bismuth.topology",
            "independent-localizers.json", "independent-ao-moments.npz")} for p in (args.left, args.right)],
        scalar_ao_dimensions=[len(sl), len(sr)], combined_integral_errors=errors,
        ao_nesting=nesting, subspaces=subspaces,
        projected_scf_hamiltonian_change_hartree=projected_h_change,
        caution="Different self-consistent potentials: a nonzero projected Hamiltonian change is not an export error")
    if args.control == "temperature":
        gaps = [float(s.energies[0, occupied]-s.energies[0, occupied-1]) for s in snapshots]
        midpoints = [float(np.mean(s.energies[0, occupied-1:occupied+1])) for s in snapshots]
        report.update(comparison_kind="controlled-electronic-temperature",
            interpretation="Temperature sensitivity at identical geometry and basis; not zero-temperature convergence",
            temperatures_kelvin=[run["options"]["temperature"] for run in runs],
            neutral_gaps_hartree=gaps, neutral_midpoints_hartree=midpoints,
            aligned_hamiltonian_change=centered_hamiltonian_change(hl, hr, spin_sl, midpoints),
            actual_inputs_differ_only_in_temperature=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ("inputs", "methods", "ao_nesting")}, indent=2))


if __name__ == "__main__":
    main()
