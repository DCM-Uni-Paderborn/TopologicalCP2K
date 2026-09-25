"""Independent finite AO localizers from complete exported SOC eigenstates.

Reintegrates overlap and first moments from Cartesian Gaussian contractions.
No orbital-centre approximation or band truncation is used. The independently
reconstructed H can then be checked against the native finite CP2K pathway.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import block_diag, eigh, eigvalsh


# Match this CP2K branch's CODATA-2006 constants, not a different library's units.
BOHR_ANGSTROM = 0.52917720859
HARTREE_EV = 2 * 10973731.568527 * 6.62606896e-34 * 299792458 / 1.602176487e-19


def moments(snapshot, primitive_integrals):
    if np.any(snapshot.periodic) or np.max(np.abs(snapshot.kpoints)) > 1e-14:
        raise ValueError("Only an isolated Gamma snapshot is supported")
    nao = sum(snapshot.atom_sizes)
    result = np.zeros((4, nao, nao), complex)
    for sa in snapshot.shells:
        for sb in snapshot.shells:
            if sa.indices[0] > sb.indices[0]:
                continue
            distance = np.linalg.norm(sa.center - sb.center)
            if distance > sa.radius + sb.radius:
                continue
            block = np.zeros((4, len(sa.indices), len(sb.indices)), complex)
            for ia, alpha in enumerate(sa.exponents):
                for ib, beta in enumerate(sb.exponents):
                    if distance > sa.radii[ia] + sb.radii[ib]:
                        continue
                    overlap = primitive_integrals(sa.center, sb.center, alpha, beta,
                                                  sa.powers, sb.powers, np.zeros(3))
                    block[0] += sa.contraction[ia].T @ overlap @ sb.contraction[ib]
                    for axis in range(3):
                        shifted = sa.powers.copy()
                        shifted[:, axis] += 1
                        moment = primitive_integrals(sa.center, sb.center, alpha, beta,
                                                     shifted, sb.powers, np.zeros(3))
                        moment += sa.center[axis] * overlap
                        block[axis + 1] += sa.contraction[ia].T @ moment @ sb.contraction[ib]
            index = np.ix_(sa.indices, sb.indices)
            for axis in range(4):
                result[axis][index] += block[axis]
                if sa.indices[0] != sb.indices[0]:
                    result[axis][np.ix_(sb.indices, sa.indices)] += block[axis].conj().T
    if np.max(np.abs(result - result.swapaxes(1, 2).conj())) > 1e-9:
        raise ValueError("Non-Hermitian integrated moments")
    return result


def operators(snapshot, integrated):
    nk, ns, nao, nb = snapshot.coefficients.shape
    if nk != 1 or ns != 2 or nb != 2 * nao or not np.array_equal(snapshot.bands, np.arange(nb)):
        raise ValueError("A complete finite AO spinor eigensystem is required")
    overlap, *xyz = integrated
    w, v = eigh(overlap)
    if w[0] <= 1e-10 * w[-1]:
        raise ValueError("Unresolved or indefinite AO metric")
    root = (v * np.sqrt(w)) @ v.conj().T
    root_inverse = (v / np.sqrt(w)) @ v.conj().T
    coefficients = snapshot.coefficients[0].reshape(nb, nb)
    frame = block_diag(root, root) @ coefficients
    error = float(np.max(np.abs(frame.conj().T @ frame - np.eye(nb))))
    if error > 1e-7:
        raise ValueError(f"Exported full frame is not S-orthonormal: {error}")
    h = (frame * snapshot.energies[0]) @ frame.conj().T
    positions = [block_diag(root_inverse @ x @ root_inverse,
                           root_inverse @ x @ root_inverse) for x in xyz]
    # Keep the physical, canonical spin basis; do not repair its TR structure.
    t = np.block([[np.zeros((nao, nao)), np.eye(nao)],
                  [-np.eye(nao), np.zeros((nao, nao))]])
    tr_residual = float(np.max(np.abs(t @ h.conj() @ t.T - h)))
    if tr_residual > 1e-8:
        raise ValueError(f"Broken reconstructed time reversal: {tr_residual}")
    return h, positions, dict(metric_minimum=float(w[0]), metric_ratio=float(w[0]/w[-1]),
                             frame_residual=error, time_reversal_residual=tr_residual)


def queries(h, positions, energies, kappas, points, aii_skew, skew_sign):
    n = len(h)
    identity = np.eye(n)
    trivial = np.diag(np.r_[np.ones(n), -np.ones(n)]).astype(complex)
    reference_sign = skew_sign(aii_skew(trivial)[0])
    result = []
    for energy in energies:
        mass = h - energy * identity
        for point in points:
            x, y = [positions[d] - point[d] * identity for d in range(2)]
            for kappa in kappas:
                cross = kappa * (x - 1j * y)
                localizer = np.block([[mass, cross], [cross.conj().T, -mass]])
                values = eigvalsh(localizer, subset_by_index=[n-2, n+1], driver="evr")
                gap = float(np.min(np.abs(values)))
                skew, skew_residual = aii_skew(localizer)
                sign = skew_sign(skew) if gap > 1e-9 else None
                row = dict(energy_hartree=float(energy), kappa_hartree_bohr=float(kappa),
                           position_bohr=list(map(float, point)), gap_hartree=gap,
                           z2=None if sign is None else (1-sign*reference_sign)//2,
                           skew_residual=skew_residual)
                result.append(row)
                print(json.dumps(row), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("--energy-offset-ev", nargs="+", type=float, default=[0.])
    parser.add_argument("--kappa", nargs="+", type=float, default=[.001, .003, .01])
    parser.add_argument("--offset", nargs="+", type=float, default=[0.])
    parser.add_argument("--output", default="independent-localizers.json")
    args = parser.parse_args()
    sys.path.insert(0, str(args.root / "build-serial"))
    from gaussian_states import read_snapshot, primitive_integrals
    from check_bloch_localizer import aii_skew, skew_sign

    work = args.work.resolve()
    snapshot_file = work / "bismuth.topology"
    snapshot = read_snapshot(snapshot_file)
    nelectron = 5 * len(snapshot.positions)
    if nelectron % 2 or snapshot.energies.shape[1] <= nelectron:
        raise ValueError("Missing conduction states or odd number of electrons")
    valence, conduction = snapshot.energies[0, nelectron-1:nelectron+1]
    neutral_midpoint = (valence + conduction) / 2
    energies = neutral_midpoint + np.asarray(args.energy_offset_ev) / HARTREE_EV
    center = np.mean(snapshot.positions, axis=0)
    points = [center + np.array([x / BOHR_ANGSTROM, 0., 0.]) for x in args.offset]
    start = time.monotonic()
    cache = work / "independent-ao-moments.npz"
    with snapshot_file.open("rb") as stream:
        snapshot_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    method_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    integral_hash = hashlib.sha256((args.root / "build-serial/gaussian_states.py").read_bytes()).hexdigest()
    if cache.exists():
        with np.load(cache, allow_pickle=False) as archive:
            if (str(archive["snapshot_sha256"]) != snapshot_hash or
                    str(archive["method_sha256"]) != method_hash or
                    str(archive["integral_sha256"]) != integral_hash):
                raise ValueError("Changed snapshot for moment cache")
            integrated = archive["moments"]
    else:
        integrated = moments(snapshot, primitive_integrals)
        with cache.open("xb") as stream:
            np.savez_compressed(stream, moments=integrated, snapshot_sha256=snapshot_hash,
                                method_sha256=method_hash, integral_sha256=integral_hash)
    h, positions, diagnostics = operators(snapshot, integrated)
    report = dict(snapshot_sha256=snapshot_hash,
                  source_sha256=method_hash, integral_sha256=integral_hash,
                  natoms=len(snapshot.positions), nao=int(sum(snapshot.atom_sizes)),
                  occupied_spinors=nelectron, spectral_gap_hartree=float(conduction-valence),
                  neutral_midpoint_hartree=float(neutral_midpoint), diagnostics=diagnostics,
                  queries=queries(h, positions, energies, args.kappa, points, aii_skew, skew_sign),
                  elapsed_seconds=time.monotonic()-start)
    with (work / args.output).open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
