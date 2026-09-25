"""Independently diagonalize printed AO SOC operators and check the NNKP spectrum."""

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.linalg import block_diag, eigvalsh

from analyze_bismuth_flakes import HARTREE_EV
from diagnose_bismuth_operators import matrix


def headerless(text, label, n):
    if text.count(label) != 1:
        raise ValueError(f"Expected exactly one {label}")
    rows = []
    for line in text.split(label)[1].splitlines():
        if not line.strip():
            continue
        try:
            row = list(map(float, line.split()))
        except ValueError:
            break
        rows.append(row)
    # NDIGITS 16 and before=4 produce two columns per printed page.
    result = np.full((n, n), np.nan)
    cursor = 0
    for col in range(0, n, 2):
        width = min(2, n-col)
        block = rows[cursor:cursor+n]
        if len(block) != n or any(len(row) != width for row in block):
            raise ValueError("Incomplete headerless matrix")
        result[:, col:col+width] = block
        cursor += n
    if cursor != len(rows) or not np.isfinite(result).all():
        raise ValueError("Unexpected headerless data")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--nao", type=int, default=26)
    parser.add_argument("--occupied", type=int, default=10)
    parser.add_argument("--omit-headers", action="store_true")
    args = parser.parse_args()
    work, n = args.work, args.nao
    text = (work / "output.out").read_text()
    if "PROGRAM ENDED AT" not in text or "*** SCF run converged" not in text:
        raise ValueError("Incomplete CP2K run")
    scalar_file = work / "soc-export-ao-soc-1_0.Log"
    scalar_text = scalar_file.read_text()
    reader = headerless if args.omit_headers else matrix
    # The post-SCF frozen-potential output contains the same overlap twice.
    overlap = reader("OVERLAP MATRIX" + scalar_text.split("OVERLAP MATRIX")[1], "OVERLAP MATRIX", n)
    h0 = reader(scalar_text, "KOHN-SHAM MATRIX", n)
    soc = []
    files = [scalar_file, work / "output.out", work / "soc-export.eig"]
    for axis in "XYZ":
        file = work / f"soc-export-ao-soc-SOC_{axis}-1_0.Log"
        files.append(file)
        if args.omit_headers:
            value = headerless(file.read_text(), "OVERLAP MATRIX", n)
        else:
            value = matrix(file.read_text(), f"AO SOC {axis} [hartree]; real antisymmetric component; multiply by i",
                           n, antisymmetric=True)
        if np.max(abs(value+value.T)) > 1e-12:
            raise ValueError("Non-antisymmetric AO SOC matrix")
        soc.append(value)
    x, y, z = soc
    h = np.block([[h0+1j*z, 1j*x-y], [1j*x+y, h0-1j*z]])
    ev = eigvalsh(h, block_diag(overlap, overlap))*HARTREE_EV
    exported = np.loadtxt(work / "soc-export.eig")
    if exported.shape != (args.occupied, 3) or not np.all(exported[:, 1] == 1):
        raise ValueError("Expected one Gamma point and the selected occupied bands")
    gap = ev[args.occupied]-ev[args.occupied-1]
    native_gap = float(re.findall(r"Sampled indirect gap \[eV\]:\s+(\S+)", text)[-1])
    error = float(np.max(abs(ev[:args.occupied]-exported[:, 2])))
    accepted = bool(error < 1e-7 and abs(gap-native_gap) < 1e-7)
    record = dict(accepted=accepted, native_gap_ev=native_gap, independent_gap_ev=float(gap),
                  eigenvalue_maximum_error_ev=error, method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    serialized = json.dumps(record, indent=2, allow_nan=False)
    with (work / "independent-soc-spectrum.json").open("x") as stream:
        stream.write(serialized + "\n")
    print(serialized)
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
