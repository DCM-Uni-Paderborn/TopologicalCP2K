"""Compare retained neutral-window indices with the native Tacho C interface.

This is a factorization check on independently reconstructed operators, not a
second native AO-assembly or SCF calculation. No entries are screened or repaired.
"""

import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy import sparse

import analyze_bismuth_flakes as analysis


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("scan", type=Path)
    parser.add_argument("library", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    root, work = args.root.resolve(), args.work.resolve()
    report = json.loads(args.scan.read_text())
    for name, expected in report["inputs"].items():
        if digest(work / name) != expected:
            raise ValueError("Changed scan input: " + name)
    if digest(Path(analysis.__file__)) != report["methods"]["operators"]:
        raise ValueError("Changed operator reconstruction")
    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import read_snapshot
    from check_bloch_localizer import aii_skew
    from check_sparse_bloch_localizer import tacho_sign_once

    for name, key in (("gaussian_states.py", "gaussian_integrals"),
                      ("check_bloch_localizer.py", "pfaffian_reference")):
        if digest(root / "build-serial" / name) != report["methods"][key]:
            raise ValueError("Changed reference method: " + name)
    snapshot = read_snapshot(work / "bismuth.topology")
    with np.load(work / "independent-ao-moments.npz", allow_pickle=False) as archive:
        h, xyz, diagnostics = analysis.operators(snapshot, archive["moments"])
    library = ctypes.CDLL(str(args.library.resolve()))
    identity = np.eye(len(h))
    trivial = np.diag(np.r_[np.ones(len(h)), -np.ones(len(h))]).astype(complex)
    reference, _, _ = tacho_sign_once(sparse.csr_matrix(aii_skew(trivial)[0]), library)
    if reference not in (-1, 1):
        raise ValueError("Native reference Pfaffian unresolved")
    rows = []
    start = time.monotonic()
    for index, expected in enumerate(report["queries"]):
        point = expected["position_bohr"]
        mass = h - expected["energy_hartree"] * identity
        cross = expected["kappa_hartree_bohr"] * (xyz[0] - point[0]*identity - 1j*(xyz[1] - point[1]*identity))
        localizer = np.block([[mass, cross], [cross.conj().T, -mass]])
        sign, residual, fill = tacho_sign_once(sparse.csr_matrix(aii_skew(localizer)[0]), library)
        invariant = None if sign is None else (1 - sign*reference) // 2
        rows.append(dict(query=index, native_z2=invariant, reference_z2=expected["z2"],
                         residual=residual, scalar_panel_entries=fill,
                         matched=invariant is not None and invariant == expected["z2"]))
        print(json.dumps(rows[-1]), flush=True)
    accepted = all(row["matched"] for row in rows)
    result = dict(kind="native-factorization-comparison-not-native-ao-assembly", accepted=accepted,
                  scan_sha256=digest(args.scan), library_sha256=digest(args.library),
                  script_sha256=digest(Path(__file__)),
                  adapter_sha256=digest(root / "build-serial/check_sparse_bloch_localizer.py"),
                  diagnostics=diagnostics, queries=rows, elapsed_seconds=time.monotonic() - start)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    if not accepted:
        raise SystemExit("Unresolved or inconsistent native Pfaffians; negative result retained")


if __name__ == "__main__":
    main()
