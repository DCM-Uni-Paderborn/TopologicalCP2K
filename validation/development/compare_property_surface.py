"""Compare converged Wilson spectra and gauge-invariant link singular values."""

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from compare_property_wilson import circular_distance, overlaps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    texts = [(directory / "output.out").read_text() for directory in [args.reference, args.candidate]]
    indices = []
    for text in texts:
        values = re.findall(r"Converged Z2 invariant:\s+(\d+)", text)
        if len(values) != 1 or "PROGRAM ENDED" not in text or "[ABORT]" in text:
            raise ValueError("Not a completed, sampling-converged surface")
        indices.append(int(values[0]))
    a, b = [np.loadtxt(d / "stanene.wilson", ndmin=2) for d in [args.reference, args.candidate]]
    if a.shape != b.shape or not np.array_equal(a[:, 0], b[:, 0]) or not np.isfinite([a, b]).all():
        raise ValueError("Different or nonfinite surfaces")
    values = [np.loadtxt(d / "stanene.eig", ndmin=2) for d in [args.reference, args.candidate]]
    if values[0].shape != values[1].shape or not np.array_equal(values[0][:, :2], values[1][:, :2]):
        raise ValueError("Different band sets")
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite energies")
    links = [overlaps(d / "stanene.mmn") for d in [args.reference, args.candidate]]
    if links[0].keys() != links[1].keys():
        raise ValueError("Different directed graphs")
    singular_error = max(float(np.max(np.abs(np.linalg.svd(m, compute_uv=False) -
                                             np.linalg.svd(links[1][key], compute_uv=False))))
                         for key, m in links[0].items())
    result = dict(indices=indices, loops=len(a), links=len(links[0]),
                  occupied_eigenvalue_error_ev=float(np.max(np.abs(values[0][:, 2] - values[1][:, 2]))),
                  phase_error=max(circular_distance(x[2:], y[2:]) for x, y in zip(a, b)),
                  loop_singular_value_error=float(np.max(np.abs(a[:, 1] - b[:, 1]))),
                  all_link_singular_value_error=singular_error,
                  scalar_diagonalizations=re.findall(r"Property scalar diagonalizations:\s+(\d+) / (\d+)", texts[1]),
                  soc_diagonalizations=re.findall(r"Property SOC diagonalizations:\s+(\d+) / (\d+)", texts[1]))
    result["tolerance"] = 2e-7
    result["accepted"] = indices[0] == indices[1] and all(result[name] <= result["tolerance"] for name in [
        "occupied_eigenvalue_error_ev", "phase_error", "loop_singular_value_error", "all_link_singular_value_error"])
    files = [d / name for d in [args.reference, args.candidate] for name in [
        "input.inp", "output.out", "stanene.wilson", "stanene.eig", "stanene.mmn"]]
    files += [Path(__file__), Path(__file__).with_name("compare_property_wilson.py")]
    result["sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "sha256"}, indent=2), flush=True)
    if not result["accepted"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
