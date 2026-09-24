"""Compare EBR dimension multisets separately at every canonical Wyckoff position.

Arguments: native reference log, pinned CSV cache, numeric reference fixture.
This does not assign individual irrep names among equal-dimensional representations.
"""
import csv
import io
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

log = Path(sys.argv[1]).read_text()
cache = Path(sys.argv[2])
fixture = iter(Path(sys.argv[3]).read_text().splitlines())
ngroups = int(next(fixture))
centering = {}
for _ in range(ngroups):
    sg, n, nf = map(int, next(fixture).split())
    operations = [next(fixture).split() for _ in range(n)]
    centering[sg] = sum(list(map(int, row[:9])) == [1, 0, 0, 0, 1, 0, 0, 0, 1] for row in operations)
    for _ in range(nf):
        next(fixture)
assert f"Independent affine-reference groups passed:         {ngroups}" in log or re.search(
    rf"Independent affine-reference groups passed:\s+{ngroups}\s", log)
families = {}
native = defaultdict(list)
for line in log.splitlines():
    if line.startswith("AFFINE_MATCH "):
        _, sg, grey, spin, label, family, operation, dimension, multiplicity, maximal = line.split()
        key = (int(sg), int(grey), int(spin), label)
        assert key not in families
        families[key] = dict(family=int(family), operation=int(operation), dimension=int(dimension),
                             multiplicity=int(multiplicity), maximal=maximal == "T")
    if line.startswith("SITE_COLUMN "):
        _, sg, grey, spin, label, local, dim, summands, maximal = line.split()
        key = (int(sg), int(grey), int(spin), label)
        assert key in families
        if maximal == "T" and int(summands) == 1:
            assert int(dim) % centering[int(sg)] == 0
            native[key].append(int(dim) // centering[int(sg)])

reference = defaultdict(list)
for sg in centering:
    for grey, directory in [(1, "elementary"), (2, "elementaryTR")]:
        table = list(csv.reader(io.StringIO((cache / f"{directory}-{sg}.csv").read_text()), delimiter="|"))
        for col in range(1, len(table[0])):
            spin = int(any("\u02e2" in row[col] or "\\bar" in row[col] for row in table[3:]))
            site = re.fullmatch(r"(\d+)([A-Za-z])\(.*\)", table[0][col])
            dimension = re.search(r"\((\d+)\)$", table[1][col])
            assert site and dimension
            key = (sg, grey, spin, site[2])
            assert key in families and families[key]["multiplicity"] == int(site[1])
            reference[key].append(int(dimension[1]))

results = []
for key, family in families.items():
    actual, expected = sorted(native[key]), sorted(reference[key])
    row = dict(sg=key[0], grey=key[1], spin=key[2], wyckoff=key[3], **family,
               native=actual, reference=expected, match=actual == expected)
    results.append(row)
    print(json.dumps(row))
summary = dict(groups=ngroups, affine_matches=len(families),
               per_site_dimension_matches=sum(row["match"] for row in results),
               mismatches=sum(not row["match"] for row in results),
               native_columns=sum(map(len, native.values())), reference_columns=sum(map(len, reference.values())))
print(json.dumps(summary))
if summary["mismatches"]:
    raise SystemExit(1)
