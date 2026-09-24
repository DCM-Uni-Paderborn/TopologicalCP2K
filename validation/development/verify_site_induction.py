"""Independent integer, path-geometry and provenance checks of Gaussian site-induction exports."""
import json
import sys
from pathlib import Path

from sympy import Matrix, Rational, eye, zeros


def verify(path):
    records = [line.split() for line in path.read_text().splitlines() if line.strip()]

    def rows(tag):
        return [row[1:] for row in records if row[0] == tag]

    sizes = rows("SITE_INDUCTION_SIZES")
    if not sizes:
        return None
    assert len(sizes) == 1
    ne, nc, composite = map(int, sizes[0])
    _, nr, na = map(int, rows("ATOMIC_SIGNATURE_SIZES")[0])
    assert na == nc
    column_records = rows("ATOMIC_SIGNATURE_COLUMN")
    assert len(column_records) == nc
    columns = {int(row[0]): (int(row[1]), int(row[2]), int(row[3]), Matrix([Rational(x) for x in row[4:]]))
               for row in column_records}
    assert set(columns) == set(range(1, nc + 1))
    inverse = {(site, local): col for col, (site, local, _, _) in columns.items()}
    assert len(inverse) == nc
    sites = {site: position for site, _, _, position in columns.values()}
    dims = Matrix([[columns[i][2] for i in range(1, nc + 1)]])
    operations = {int(row[0]): (Matrix(3, 3, list(map(int, row[2:11]))), Matrix([Rational(x) for x in row[11:]]))
                  for row in rows("ATOMIC_REFERENCE_OPERATION")}

    def periodic_zero(vector):
        return all(abs(value - round(value)) < Rational("1e-8") for value in vector)

    def sparse(tag, shape):
        result = zeros(*shape)
        used = set()
        for row in rows(tag):
            i, j, value = map(int, row)
            assert 1 <= i <= shape[0] and 1 <= j <= shape[1]
            assert (i, j) not in used and value >= 0
            used.add((i, j))
            result[i - 1, j - 1] = value
        return result

    a = sparse("ATOMIC_SIGNATURE_ENTRY", (nr, nc))
    expansion = sparse("SITE_INDUCTION_EXPANSION", (nc, nc))
    assert a*expansion == a and dims*expansion == dims
    assert sum(bool(sum(expansion[:, j]) > 1) for j in range(nc)) == composite
    parent = list(range(nc + 1))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    edges = {}
    for row in rows("SITE_INDUCTION_EDGE"):
        e, source, target, op, *shift = map(int, row[:7])
        y = Matrix([Rational(x) for x in row[7:]])
        rotation, translation = operations[op]
        assert max(map(abs, y - rotation*sites[target] - translation - Matrix(shift))) < Rational("1e-8")
        hs = [g for g, (w, tau) in operations.items() if periodic_zero(w*sites[source] + tau - sites[source])]
        ht = [g for g, (w, tau) in operations.items() if periodic_zero(w*y + tau - y)]
        assert set(hs) <= set(ht) and len(ht) % len(hs) == 0
        for g in hs:
            displacement = (operations[g][0] - eye(3))*(sites[source] - y)
            assert max(map(abs, displacement)) < Rational("1e-8")
        src = [inverse[source, j] for j in range(1, 1 + sum(site == source for site, _, _, _ in columns.values()))]
        dst = [inverse[target, j] for j in range(1, 1 + sum(site == target for site, _, _, _ in columns.values()))]
        b = zeros(len(dst), len(src))
        used = set()
        for ee, i, j, value in map(lambda row: tuple(map(int, row)), rows("SITE_INDUCTION_ENTRY")):
            if ee != e:
                continue
            assert 1 <= i <= len(dst) and 1 <= j <= len(src) and value > 0
            assert (i, j) not in used
            used.add((i, j))
            b[i - 1, j - 1] = value
        src0, dst0 = [i - 1 for i in src], [i - 1 for i in dst]
        assert a[:, dst0]*b == a[:, src0] and dims[:, dst0]*b == dims[:, src0]
        for j in range(len(src)):
            if sum(b[:, j]) == 1:
                i = list(b[:, j]).index(1)
                parent[root(src[j])] = root(dst[i])
        assert e not in edges
        edges[e] = src, dst, b
    assert set(edges) == set(range(1, ne + 1))
    labels = {int(row[0]): tuple(map(int, row[1:])) for row in rows("SITE_INDUCTION_COLUMN")}
    assert set(labels) == set(columns)
    for j, (representative, edge, local) in labels.items():
        assert root(j) == root(representative)
        assert labels[representative][0] == representative
        assert expansion[:, j - 1] == expansion[:, representative - 1]
        if edge == 0:
            assert local == 0 and sum(expansion[:, j - 1]) == 1
        else:
            src, dst, b = edges[edge]
            assert 1 <= local <= len(src) and root(src[local - 1]) == root(j)
            assert sum(b[:, local - 1]) > 1
            assert expansion[:, [i - 1 for i in dst]]*b[:, local - 1] == expansion[:, j - 1]
    result = dict(file=path.name, columns=nc, relations=ne, composite=composite, verified=True)
    print(json.dumps(result), flush=True)
    return result


results = [verify(Path(path)) for path in sys.argv[1:]]
assert any(row is not None for row in results)
print(json.dumps(dict(verified=sum(row is not None for row in results))))
