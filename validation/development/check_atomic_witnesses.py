"""Independently verify exported signature witnesses with Python's exact integers."""

import argparse
from pathlib import Path


def check(path):
    rows, entries, columns, witnesses = {}, {}, {}, {}
    witnesses["SIGNED_ATOMIC_WITNESS"] = {}
    witnesses["NONNEGATIVE_ATOMIC_WITNESS"] = {}
    sizes, result = None, None
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        key = fields[0]
        if key == "ATOMIC_SIGNATURE_SIZES":
            sizes = tuple(map(int, fields[1:]))
        elif key == "ATOMIC_SIGNATURE_ROW":
            index, point, local, dimension, target = map(int, fields[1:])
            assert index not in rows
            rows[index] = (point, dimension, target)
        elif key == "ATOMIC_SIGNATURE_COLUMN":
            index, site, local, dimension = map(int, fields[1:5])
            assert index not in columns
            columns[index] = dimension
        elif key == "ATOMIC_SIGNATURE_ENTRY":
            row, column, value = map(int, fields[1:])
            assert (row, column) not in entries
            entries[row, column] = value
        elif key in witnesses:
            column, value = map(int, fields[1:])
            assert column not in witnesses[key]
            witnesses[key][column] = value
        elif key == "ATOMIC_SIGNATURE_RESULT":
            result = tuple(map(int, fields[1:]))
    assert sizes is not None and result is not None, path
    np, m, n = sizes
    assert set(rows) == set(range(1, m + 1))
    assert set(columns) == set(range(1, n + 1))
    assert set(entries) == {(i, j) for i in rows for j in columns}
    assert all(value >= 0 for value in entries.values())
    ranks = []
    for point in range(1, np + 1):
        block = [i for i in rows if rows[i][0] == point]
        assert block
        ranks.append(sum(rows[i][1] * rows[i][2] for i in block))
        for j in columns:
            assert sum(rows[i][1] * entries[i, j] for i in block) == columns[j]
    assert len(set(ranks)) == 1
    checked = 0
    for key, vector in witnesses.items():
        if not vector:
            continue
        assert set(vector) == set(columns)
        if key == "NONNEGATIVE_ATOMIC_WITNESS":
            assert all(value >= 0 for value in vector.values())
        for i in rows:
            assert sum(entries[i, j] * vector[j] for j in columns) == rows[i][2]
        checked += 1
    if result[0] == 1:
        assert checked > 0
    if result[1] == 1:
        assert witnesses["NONNEGATIVE_ATOMIC_WITNESS"]
    if result[0] == 0:
        assert result[1] == 0 and checked == 0
    print(path.name, "size", sizes, "result", result, "verified witnesses", checked)
    return checked


parser = argparse.ArgumentParser()
parser.add_argument("directory", type=Path)
args = parser.parse_args()
files = sorted(args.directory.glob("*.little_group"))
assert files
count = sum(check(path) for path in files)
print("Verified exported matrices/witnesses:", len(files), count)
