"""Bound a joint energy/position/scale region for a retained finite Bi snapshot."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigvalsh, svdvals

from localizer_region import cover_region, variation_bound


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("scan", type=Path)
    parser.add_argument("methods", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--kappa", nargs=2, type=float, required=True)
    parser.add_argument("--energy-half-fraction", type=float, default=.02)
    parser.add_argument("--position-half-angstrom", type=float, default=.05)
    parser.add_argument("--max-nodes", type=int, default=255)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if (not np.isfinite([args.energy_half_fraction, args.position_half_angstrom]).all()
            or not 0 <= args.energy_half_fraction < .5 or args.position_half_angstrom < 0):
        raise ValueError("Invalid physical region")
    root, work = args.root.resolve(), args.work.resolve()
    sys.path.insert(0, str(args.methods.resolve()))
    sys.path.insert(0, str(root / "build-serial"))
    import analyze_bismuth_flakes as analysis
    from gaussian_states import read_snapshot
    from check_bloch_localizer import aii_skew, skew_sign
    report = json.loads(args.scan.read_text())
    for name, expected in report["inputs"].items():
        if digest(work / name) != expected:
            raise ValueError("Changed native source: " + name)
    for path, key in ((Path(analysis.__file__), "operators"),
                      (root / "build-serial/gaussian_states.py", "gaussian_integrals"),
                      (root / "build-serial/check_bloch_localizer.py", "pfaffian_reference")):
        if digest(path) != report["methods"][key]:
            raise ValueError("Changed mathematical reference: " + key)
    snapshot = read_snapshot(work / "bismuth.topology")
    with np.load(work / "independent-ao-moments.npz", allow_pickle=False) as cache:
        h, xyz, diagnostics = analysis.operators(snapshot, cache["moments"])
    center = np.mean(snapshot.positions, axis=0)
    n = len(h)
    eye = np.eye(n)
    e0 = report["spectrum"]["midpoint_hartree"]
    de = args.energy_half_fraction * report["spectrum"]["gap_hartree"]
    dr = args.position_half_angstrom / analysis.BOHR_ANGSTROM
    lo = np.r_[e0 - de, args.kappa[0], center[:2] - dr]
    hi = np.r_[e0 + de, args.kappa[1], center[:2] + dr]
    z0 = xyz[0] - center[0]*eye - 1j*(xyz[1] - center[1]*eye)
    rho = float(svdvals(z0)[0])

    def localizer(p):
        e, k, x, y = p
        a = h - e*eye
        b = k*(z0 - ((x-center[0])-1j*(y-center[1]))*eye)
        return np.block([[a, b], [b.conj().T, -a]])

    def gap(p):
        values = eigvalsh(localizer(p), subset_by_index=[n-2, n+1], driver="evr")
        return float(np.min(abs(values)))

    started = time.monotonic()
    bounds = cover_region(lo, hi, rho, center[:2], gap, max_nodes=args.max_nodes)
    anchor = (lo + hi)/2
    anchor_matrix = localizer(anchor)
    trivial = np.diag(np.r_[np.ones(n), -np.ones(n)]).astype(complex)
    reference = skew_sign(aii_skew(trivial)[0])
    variation, _ = variation_bound(lo, hi, rho, center[:2])
    checks = []
    points = [tuple(anchor)] + sorted(set(product(*zip(lo, hi))))
    for point in dict.fromkeys(points):
        matrix = localizer(point)
        # Full divide-and-conquer spectrum is independent of the middle-four
        # eigenvalue path used by the adaptive subdivision.
        spectrum = eigvalsh(matrix, driver="evd")
        exact_gap = float(np.min(abs(spectrum)))
        pairing = float(np.max(abs(spectrum + spectrum[::-1])))
        partial_error = abs(gap(point) - exact_gap)
        distance = float(np.max(abs(eigvalsh(matrix-anchor_matrix, driver="evd"))))
        skew, skew_residual = aii_skew(matrix)
        index = (1-skew_sign(skew)*reference)//2 if exact_gap > 1e-9 else None
        if max(pairing, partial_error, distance-variation) > 1e-10:
            raise ValueError("Independent spectral or perturbation audit failed")
        checks.append(dict(parameters=list(point), gap_hartree=exact_gap, z2=index,
                           spectral_pairing_error_hartree=pairing,
                           middle_four_gap_error_hartree=partial_error,
                           actual_perturbation_norm_hartree=distance,
                           skew_residual=skew_residual))
    anchor_gap, anchor_index = checks[0]["gap_hartree"], checks[0]["z2"]
    if bounds["resolved"] and any(row["z2"] != anchor_index for row in checks):
        raise ValueError("Covered region has inconsistent independently sampled indices")
    import localizer_region
    result = dict(case=work.name, definition="Frozen finite AO localizer in a fixed physical metric",
                  parameter_order=["energy_hartree", "kappa_hartree_per_bohr", "x_bohr", "y_bohr"],
                  lower=lo.tolist(), upper=hi.tolist(), reference_xy_bohr=center[:2].tolist(),
                  position_norm_bohr=rho, energy_half_fraction=args.energy_half_fraction,
                  position_half_angstrom=args.position_half_angstrom, diagnostics=diagnostics,
                  anchor=anchor.tolist(), anchor_gap_hartree=anchor_gap, anchor_z2=anchor_index,
                  resolved=bounds["resolved"] and anchor_index is not None,
                  region_z2=anchor_index if bounds["resolved"] else None, bound=bounds,
                  corner_audit=dict(checks=checks, variation_bound_hartree=variation,
                                    tolerance_hartree=1e-10, accepted=True),
                  minimum_gap_lower_bound_hartree=min((r["lower_bound"] for r in bounds["covered"]), default=None),
                  inputs=report["inputs"], original_methods=report["methods"],
                  scan_sha256=digest(args.scan),
                  methods={"driver":digest(Path(__file__)), "region":digest(Path(localizer_region.__file__))},
                  elapsed_seconds=time.monotonic()-started)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k not in
                     ("bound", "corner_audit", "inputs", "original_methods", "methods")}, indent=2))
    if not result["resolved"]:
        raise SystemExit("Region remains unresolved; no common-index conclusion")


if __name__ == "__main__":
    main()
