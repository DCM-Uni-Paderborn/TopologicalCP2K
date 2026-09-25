"""Check diagnostic class bookkeeping without treating it as an atomic-column reduction."""
import argparse
import json
from pathlib import Path

from verify_onsite_export import read


def audit(directory):
    results = {}
    for path in sorted(directory.glob("*.little_group")):
        records, _ = read(path)
        size, = records["ATOMIC_SIGNATURE_SIZES"]
        _, nr, nc = map(int, size)
        classes = {int(j): int(rep) for j, rep in records["ATOMIC_CHARACTER_CLASS"]}
        assert len(records["ATOMIC_CHARACTER_CLASS"]) == nc
        assert set(classes) == set(range(1, nc + 1))
        matrix = {(int(i), int(j)): int(value) for i, j, value in records["ATOMIC_SIGNATURE_ENTRY"]}
        assert len(matrix) == nr * nc
        for j, rep in classes.items():
            assert 1 <= rep <= j and classes[rep] == rep
            assert all(matrix[i, j] == matrix[i, rep] for i in range(1, nr + 1))
        error, = records["ATOMIC_CHARACTER_CLASS_RESIDUAL"]
        assert 0 <= float(error[0]) < 1e-7
        results[path.name] = dict(columns=nc, classes=len(set(classes.values())),
                                 sampled_classes=len({tuple(matrix[i, j] for i in range(1, nr + 1)) for j in classes}),
                                 residual=float(error[0]))
    assert len(results) == 16
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial", type=Path)
    parser.add_argument("mpi", type=Path)
    args = parser.parse_args()
    serial, mpi = audit(args.serial), audit(args.mpi)
    assert serial.keys() == mpi.keys()
    for name in serial:
        for key in ("columns", "classes", "sampled_classes"):
            assert serial[name][key] == mpi[name][key]
    print(json.dumps(dict(serial=serial, mpi=mpi), indent=2))


if __name__ == "__main__":
    main()
