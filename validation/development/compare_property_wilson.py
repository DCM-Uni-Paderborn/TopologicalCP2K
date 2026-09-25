"""Compare reconstructed property frames and directed links in the physical AO metric."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def overlaps(path):
    with path.open() as stream:
        next(stream)
        nb, nk, nn = map(int, next(stream).split())
        if min(nb, nk, nn) < 1:
            raise ValueError("Invalid overlap dimensions")
        result = {}
        for _ in range(nk * nn):
            key = tuple(map(int, next(stream).split()))
            values = np.array([list(map(float, next(stream).split())) for _ in range(nb * nb)])
            if len(key) != 5 or not (1 <= key[0] <= nk and 1 <= key[1] <= nk):
                raise ValueError("Invalid directed link")
            if values.shape != (nb * nb, 2) or not np.isfinite(values).all():
                raise ValueError("Invalid overlap entries")
            if key in result:
                raise ValueError("Duplicate directed link")
            result[key] = (values[:, 0] + 1j * values[:, 1]).reshape((nb, nb), order="F")
        if stream.read().strip():
            raise ValueError("Trailing overlap data")
    return result


def circular_distance(a, b):
    if a.shape != b.shape or not len(a) or not np.isfinite([a, b]).all():
        raise ValueError("Invalid circular spectrum")
    a, b = np.sort(a % 1.0), np.sort(b % 1.0)
    return min(float(np.max(np.abs((a - np.roll(b, shift) + 0.5) % 1.0 - 0.5)))
               for shift in range(len(b)))


def frame_errors(metric, left, right):
    if left.shape != right.shape or left.ndim != 3 or metric.shape != (left.shape[1], left.shape[1]):
        raise ValueError("Incompatible frame dimensions")
    if not all(np.isfinite(a).all() for a in [metric, left, right]):
        raise ValueError("Nonfinite frame or metric")
    gauge = sum(a.conj().T @ metric @ b for a, b in zip(left, right))
    eye = np.eye(gauge.shape[0])
    norm_error = max(float(np.max(np.abs(sum(a.conj().T @ metric @ a for a in frame) - eye)))
                     for frame in [left, right])
    difference = right - left @ gauge
    residual = sum(a.conj().T @ metric @ a for a in difference)
    return gauge, dict(frame_norm_error=norm_error,
                       subspace_residual=float(np.sqrt(np.max(np.abs(residual)))),
                       sewing_unitarity_error=float(np.max(np.abs(gauge.conj().T @ gauge - eye))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("reference", type=Path, help="Full-diagonalization seed, without suffix")
    parser.add_argument("candidate", type=Path, help="Reconstructed seed, without suffix")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    sys.path.insert(0, str(args.root / "build-serial"))
    from gaussian_states import cross_ao, read_snapshot

    ref = read_snapshot(args.reference.with_suffix(".topology"))
    candidate = read_snapshot(args.candidate.with_suffix(".topology"))
    if not np.array_equal(ref.kpoints, candidate.kpoints) or not np.array_equal(ref.bands, candidate.bands):
        raise ValueError("Changed property graph or band selection")
    for name in ["cell", "periodic", "positions", "kinds", "atom_sizes"]:
        if not np.array_equal(getattr(ref, name), getattr(candidate, name)):
            raise ValueError("Changed geometry or AO layout: " + name)
    if ref.coefficients.shape != candidate.coefficients.shape or ref.energies.shape != candidate.energies.shape:
        raise ValueError("Changed state dimensions")
    result = dict(reference=str(args.reference), candidate=str(args.candidate),
                  eigenvalue_error_hartree=float(np.max(np.abs(ref.energies - candidate.energies))),
                  frame_norm_error=0.0, subspace_residual=0.0, sewing_unitarity_error=0.0,
                  link_covariance_error=0.0, wilson_phase_error=0.0)
    gauges = []
    for i, k in enumerate(ref.kpoints):
        metric = cross_ao(ref, ref, k, k)
        gauge, errors = frame_errors(metric, ref.coefficients[i], candidate.coefficients[i])
        gauges.append(gauge)
        for name, value in errors.items():
            result[name] = max(result[name], value)
    left, right = overlaps(args.reference.with_suffix(".mmn")), overlaps(args.candidate.with_suffix(".mmn"))
    if left.keys() != right.keys():
        raise ValueError("Directed links or reciprocal boundary shifts differ")
    for key, matrix in left.items():
        expected = gauges[key[0] - 1].conj().T @ matrix @ gauges[key[1] - 1]
        result["link_covariance_error"] = max(result["link_covariance_error"],
                                               float(np.max(np.abs(expected - right[key]))))
    suffixes = [".topology", ".mmn"]
    result["wilson_loops"] = 0
    result["link_singular_value_error"] = 0.0
    if args.reference.with_suffix(".wilson").exists() or args.candidate.with_suffix(".wilson").exists():
        a = np.loadtxt(args.reference.with_suffix(".wilson"), ndmin=2)
        b = np.loadtxt(args.candidate.with_suffix(".wilson"), ndmin=2)
        if a.shape != b.shape or not np.array_equal(a[:, 0], b[:, 0]) or not np.isfinite([a, b]).all():
            raise ValueError("Wilson loop count differs")
        result["wilson_phase_error"] = max(circular_distance(x[2:], y[2:]) for x, y in zip(a, b))
        result["link_singular_value_error"] = float(np.max(np.abs(a[:, 1] - b[:, 1])))
        result["wilson_loops"] = len(a)
        suffixes.append(".wilson")
    result["points"] = len(ref.kpoints)
    result["links"] = len(left)
    result["tolerance"] = 2e-7
    result["accepted"] = all(result[key] <= result["tolerance"] for key in [
        "eigenvalue_error_hartree", "frame_norm_error", "subspace_residual", "sewing_unitarity_error",
        "link_covariance_error", "wilson_phase_error", "link_singular_value_error"])
    paths = [args.reference.with_suffix(suffix) for suffix in suffixes]
    paths += [args.candidate.with_suffix(suffix) for suffix in suffixes]
    paths += [Path(__file__), args.root / "build-serial/gaussian_states.py"]
    result["sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "sha256"}, indent=2), flush=True)
    if not result["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
