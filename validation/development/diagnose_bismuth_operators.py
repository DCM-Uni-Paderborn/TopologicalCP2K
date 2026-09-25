"""Compare printed finite AO operators with the full exported SOC eigensystem."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import block_diag, eigh


def matrix(text, label, n, antisymmetric=False):
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip() == label]
    if len(starts) != 1:
        raise ValueError(f"Expected exactly one {label}")
    result = np.full((n, n), np.nan)
    columns = []
    for line in lines[starts[0]+1:]:
        fields = line.split()
        if fields and all(field.isdigit() for field in fields):
            columns = [int(field)-1 for field in fields]
        elif columns and len(fields) == 4+len(columns) and fields[0].isdigit():
            row = int(fields[0])-1
            result[row, columns] = list(map(float, fields[4:]))
            if row == n-1 and columns[-1] == n-1:
                break
        elif columns and antisymmetric and len(fields) == 1+len(columns) and fields[0].isdigit():
            row = int(fields[0])-1
            result[row, columns] = list(map(float, fields[1:]))
            if row == n-1 and columns[-1] == n-1:
                break
    sign = -1 if antisymmetric else 1
    if not np.isfinite(result).all() or np.max(abs(result-sign*result.T)) > 1e-10:
        raise ValueError(f"Incomplete or incorrect transpose symmetry: {label}")
    return result


def soc_components(h):
    n = len(h)//2
    if h.shape != (2*n, 2*n) or np.max(abs(h-h.conj().T)) > 1e-10:
        raise ValueError("Expected an even Hermitian spinor matrix")
    uu, ud, du, dd = h[:n, :n], h[:n, n:], h[n:, :n], h[n:, n:]
    scalar = (uu+dd)/2
    components = np.array([(ud+du)/(2j), (du-ud)/2, (uu-dd)/(2j)])
    if max(np.max(abs(scalar.imag)), np.max(abs(components.imag)),
           np.max(abs(components+components.transpose(0, 2, 1)))) > 1e-10:
        raise ValueError("Spinor matrix is not real scalar plus imaginary antisymmetric SOC")
    return scalar.real, components.real


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.root / "build-serial"))
    from gaussian_states import read_snapshot

    run = json.loads((args.work / "run.json").read_text())
    if not run["completed"] or not run["scf_converged"] or run["returncode"]:
        raise ValueError("Only completed converged runs can be compared")
    snapshot_file = args.work / "bismuth.topology"
    snapshot_hash = hashlib.sha256(snapshot_file.read_bytes()).hexdigest()
    if snapshot_hash != run["files"][snapshot_file.name]:
        raise ValueError("Changed SOC snapshot")
    snapshot = read_snapshot(snapshot_file)
    with np.load(args.work / "independent-ao-moments.npz", allow_pickle=False) as cache:
        if str(cache["snapshot_sha256"]) != snapshot_hash:
            raise ValueError("Moment cache belongs to a different SOC snapshot")
        integrated = cache["moments"]
    n = len(integrated[0])
    s = integrated[0]
    values, vectors = eigh(s)
    inverse = (vectors / np.sqrt(values)) @ vectors.conj().T
    coefficients = snapshot.coefficients[0].reshape(2*n, 2*n)
    covariant_frame = block_diag(s, s) @ coefficients
    h = (covariant_frame * snapshot.energies[0]) @ covariant_frame.conj().T
    scalar, soc = soc_components(h)
    report = dict(nao=n, scalar_imaginary_max=float(np.max(abs(scalar.imag))), matrices=[])
    for path in sorted(args.work.glob("bismuth-ao-*.Log")):
        if hashlib.sha256(path.read_bytes()).hexdigest() != run["files"][path.name]:
            raise ValueError("Printed AO matrix changed since the run")
        text = path.read_text()
        for label, expected in [("OVERLAP MATRIX", s), ("KOHN-SHAM MATRIX", scalar),
                                *(("AO POSITION " + axis + " [bohr]; coordinate origin; covariant AO matrix",
                                   integrated[i+1]) for i, axis in enumerate("XYZ")),
                                *(("AO SOC " + axis + " [hartree]; real antisymmetric component; multiply by i",
                                   soc[i]) for i, axis in enumerate("XYZ"))]:
            if label not in text:
                continue
            if "-SOC_" in path.name and not label.startswith("AO SOC"):
                continue
            for occurrence, block in enumerate(text.split(label)[1:]):
                actual = matrix(label+block, label, n, antisymmetric=label.startswith("AO SOC"))
                difference = actual-expected
                row, col = np.unravel_index(np.argmax(abs(difference)), difference.shape)
                record = dict(file=path.name, label=label, occurrence=occurrence,
                              maximum_error=float(np.max(abs(difference))),
                              metric_error=float(np.max(abs(inverse @ difference @ inverse))),
                              largest_entry=[int(row), int(col)], native=float(actual[row, col]),
                              reference=[float(expected[row, col].real), float(expected[row, col].imag)])
                report["matrices"].append(record)
                print(json.dumps(record), flush=True)
    with (args.work / "operator-diagnostic.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
