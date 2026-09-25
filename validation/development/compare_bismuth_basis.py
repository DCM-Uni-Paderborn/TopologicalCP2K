"""Physical cross-basis comparisons of complete finite SOC exports.

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
from scipy.linalg import block_diag, eigh, solve, svdvals

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
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root = args.root.resolve()
    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import primitive_integrals, read_snapshot

    integral_hash = digest(root / "build-serial/gaussian_states.py")
    runs, snapshots, own = zip(*(load_case(p.resolve(), read_snapshot, integral_hash)
                               for p in (args.left, args.right)))
    for key in ("source_commit", "executable_sha256", "library_sha256", "sources", "runtime"):
        if runs[0]["provenance"][key] != runs[1]["provenance"][key]:
            raise ValueError("Uncontrolled change in native provenance: " + key)
    ignored = {"root", "output", "basis", "restart"}
    for key in set(runs[0]["options"]) | set(runs[1]["options"]):
        if key not in ignored and runs[0]["options"].get(key) != runs[1]["options"].get(key):
            raise ValueError("Uncontrolled change in calculation options: " + key)
    left, right = snapshots
    ml, mr, cross = cross_moments(left, right, primitive_integrals)
    errors = [float(np.max(abs(a-b))) for a, b in zip(own, (ml, mr))]
    if max(errors) > 1e-10:
        raise ValueError("Combined integration differs from independently cached integrals")
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
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(ao_nesting=nesting, subspaces=subspaces,
                         projected_scf_hamiltonian_change_hartree=projected_h_change), indent=2))


if __name__ == "__main__":
    main()
