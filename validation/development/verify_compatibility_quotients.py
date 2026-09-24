"""Independent exact SymPy checks of CP2K's exported compatibility quotients."""
from pathlib import Path
import json
import sys

from sympy import Matrix, ZZ, zeros
from sympy.matrices.normalforms import smith_normal_form


def verify(path):
    records = [line.split() for line in path.read_text().splitlines()]
    records = [row for row in records if row]
    def rows(tag):
        return [list(map(int, row[1:])) for row in records if row[0] == tag]
    if not rows("COMPATIBILITY_LATTICE_STATUS"):
        return None
    assert rows("COMPATIBILITY_LATTICE_STATUS") == [[0]], path
    m, n, q = rows("COMPATIBILITY_LATTICE_SIZES")[0]
    _, atomic_rows, na = rows("ATOMIC_SIGNATURE_SIZES")[0]
    assert n == atomic_rows
    def sparse(tag, shape):
        matrix = zeros(*shape)
        seen = set()
        for i, j, value in rows(tag):
            assert (i, j) not in seen
            seen.add((i, j))
            matrix[i-1, j-1] = value
        return matrix
    c = sparse("COMPATIBILITY_CONSTRAINT", (m, n))
    k = sparse("COMPATIBILITY_KERNEL", (n, q))
    a = sparse("ATOMIC_SIGNATURE_ENTRY", (n, na))
    y = sparse("COMPATIBILITY_ATOMIC_COORDINATE", (q, na))
    u = sparse("COMPATIBILITY_QUOTIENT_MAP", (q, q))
    invariants = [value for _, value in rows("COMPATIBILITY_QUOTIENT_FACTOR")]
    b = Matrix([row[-1] for row in rows("ATOMIC_SIGNATURE_ROW")])
    assert len(invariants) == q
    assert c*k == zeros(m, q)
    assert k*y == a and c*a == zeros(m, na)
    assert c.to_DM().rank() == n-q
    kd = smith_normal_form(k, domain=ZZ)
    assert all(abs(kd[i, i]) == 1 for i in range(q)), "Kernel not saturated"
    yd = smith_normal_form(y, domain=ZZ)
    independent = [abs(yd[i, i]) if i < min(q, na) else 0 for i in range(q)]
    assert list(map(int, independent)) == invariants
    assert abs(u.det(method="domain-ge")) == 1
    transformed = u*y
    r = sum(d > 0 for d in invariants)
    normalized = zeros(r, na)
    for i, d in enumerate(invariants):
        for j in range(na):
            if d == 0:
                assert transformed[i, j] == 0
            else:
                assert transformed[i, j] % d == 0
                normalized[i, j] = transformed[i, j] // d
    nd = smith_normal_form(normalized, domain=ZZ)
    assert all(abs(nd[i, i]) == 1 for i in range(r)), "Quotient row map is not complete"
    indices = rows("COMPATIBILITY_QUOTIENT_CLASS")
    compatible = c*b == zeros(m, 1)
    if compatible:
        z, params = k.gauss_jordan_solve(b)
        assert params.rows == 0 and all(value.q == 1 for value in z)
        value = u*z
        expected = [int(value[i] % d if d else value[i]) for i, d in enumerate(invariants)]
        assert indices == [[i+1, v] for i, v in enumerate(expected)]
    else:
        assert not indices, "Incompatible target has a quotient class"
    result = dict(file=path.name, equations=m, coordinates=n, compatible_rank=q,
                  free_rank=invariants.count(0), torsion=[d for d in invariants if d > 1],
                  compatible=compatible, zero_class=all(v == 0 for _, v in indices) if compatible else None)
    print(json.dumps(result), flush=True)
    return result


results = [verify(Path(path)) for path in sys.argv[1:]]
assert any(result is not None for result in results), "No quotient records found"
