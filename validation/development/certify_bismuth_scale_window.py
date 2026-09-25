"""Bound the finite localizer gap between sampled kappa values.

For L(kappa)=L(0)+kappa*D, the gap is Lipschitz with constant ||D||.
These numerical bounds include a fixed roundoff margin; they are not
interval-arithmetic certificates or statements about the infinite material.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigvalsh, svdvals

import analyze_bismuth_flakes as analysis


def cover_gap(lower, upper, lipschitz, evaluate, margin=1e-10, max_nodes=1023):
    if not np.isfinite([lower, upper, lipschitz, margin]).all() or not 0 < lower < upper or lipschitz < 0 or margin <= 0:
        raise ValueError("Invalid scale interval or norm bound")
    pending, covered, unresolved, evaluations = [(lower, upper)], [], [], []
    while pending:
        left, right = pending.pop()
        if len(evaluations) >= max_nodes:
            unresolved.append([left, right])
            continue
        middle = (left + right) / 2
        gap = float(evaluate(middle))
        if not np.isfinite(gap) or gap < 0:
            raise ValueError("Invalid evaluated gap")
        bound = gap - lipschitz * (right - left) / 2 - margin
        row = dict(lower=left, upper=right, midpoint=middle, gap=gap, lower_bound=bound)
        evaluations.append(row)
        if bound > 0:
            covered.append(row)
        elif right - left <= 1e-10 or middle in (left, right):
            unresolved.append([left, right])
        else:
            pending.extend([(middle, right), (left, middle)])
    covered.sort(key=lambda row: row["lower"])
    return dict(covered=covered, unresolved=unresolved, evaluations=evaluations,
                resolved=not unresolved, roundoff_margin_hartree=margin)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("scan", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--interval", nargs=2, type=float, default=[.0025, .0035])
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root, work = args.root.resolve(), args.work.resolve()
    report = json.loads(args.scan.read_text())
    for name, expected in report["inputs"].items():
        if digest(work / name) != expected:
            raise ValueError("Changed reference input: " + name)
    if digest(Path(analysis.__file__)) != report["methods"]["operators"]:
        raise ValueError("Changed operator method")
    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import read_snapshot
    from check_bloch_localizer import aii_skew, skew_sign

    for name, key in (("gaussian_states.py", "gaussian_integrals"),
                      ("check_bloch_localizer.py", "pfaffian_reference")):
        if digest(root / "build-serial" / name) != report["methods"][key]:
            raise ValueError("Changed reference method: " + name)
    snapshot = read_snapshot(work / "bismuth.topology")
    with np.load(work / "independent-ao-moments.npz", allow_pickle=False) as archive:
        h, xyz, diagnostics = analysis.operators(snapshot, archive["moments"])
    identity = np.eye(len(h))
    center = np.mean(snapshot.positions, axis=0)
    mass = h - report["spectrum"]["midpoint_hartree"] * identity
    cross = xyz[0] - center[0]*identity - 1j*(xyz[1] - center[1]*identity)
    lipschitz = float(svdvals(cross)[0])
    def localizer(kappa):
        return np.block([[mass, kappa*cross], [kappa*cross.conj().T, -mass]])

    def gap(kappa):
        values = eigvalsh(localizer(kappa), subset_by_index=[len(h)-2, len(h)+1], driver="evr")
        return float(np.min(abs(values)))

    started = time.monotonic()
    bounds = cover_gap(*args.interval, lipschitz, gap)
    midpoint = sum(args.interval) / 2
    trivial = np.diag(np.r_[np.ones(len(h)), -np.ones(len(h))]).astype(complex)
    reference_sign = skew_sign(aii_skew(trivial)[0])
    anchor_sign = skew_sign(aii_skew(localizer(midpoint))[0]) if gap(midpoint) > 1e-9 else None
    result = dict(case=work.name, interval=args.interval, lipschitz_bohr=lipschitz,
                  reference_scan_sha256=digest(args.scan), script_sha256=digest(Path(__file__)),
                  methods=report["methods"], inputs=report["inputs"], diagnostics=diagnostics,
                  anchor_kappa=midpoint, anchor_z2=None if anchor_sign is None else (1-anchor_sign*reference_sign)//2,
                  minimum_gap_lower_bound_hartree=min((r["lower_bound"] for r in bounds["covered"]), default=None),
                  bound=bounds, elapsed_seconds=time.monotonic()-started)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("bound", "methods", "inputs")}, indent=2))
    if not bounds["resolved"] or anchor_sign is None:
        raise SystemExit("Unresolved scale interval retained; no stability conclusion")


if __name__ == "__main__":
    main()
