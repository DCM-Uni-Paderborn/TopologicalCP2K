"""Audit tabulated characters by their central projectors in the twisted algebra."""
import json
import sys
from pathlib import Path

import numpy as np


def complex_row(line):
    return np.array([complex(*map(float, word.strip("()").split(",")))
                     for word in line.split()])


def projector_error(product, factor, characters):
    n = len(product)
    error = 0.0
    for chi in characters:
        dim = chi[0].real
        matrix = np.zeros((n, n), dtype=complex)
        matrix[product, np.arange(n)] = (dim/n*chi.conj())[:, None]*factor
        error = max(error, np.max(np.abs(matrix-matrix.conj().T)),
                    np.max(np.abs(matrix@matrix-matrix)),
                    abs(np.trace(matrix)-dim**2))
        for g in range(n):
            left = np.empty_like(matrix)
            left[product[g], :] = factor[g, :, None]*matrix
            right = matrix[:, product[g]]*factor[g]
            error = max(error, np.max(np.abs(left-right)))
    return float(error)


lines = iter(Path(sys.argv[1]).read_text().splitlines())
counts = dict(cases=0, raw_failed=0, conjugate_failed=0, both_failed=0)
for _ in range(int(next(lines))):
    sg, spin, n, ni, point = next(lines).split()
    n, ni = int(n), int(ni)
    next(lines)
    for _ in range(n):
        next(lines)
    product = np.empty((n, n), dtype=int)
    factor = np.empty((n, n), dtype=complex)
    for j in range(n):
        product[:, j] = np.array(list(map(int, next(lines).split())))-1
        factor[:, j] = complex_row(next(lines))
    chars = []
    for j in range(ni):
        next(lines)
        chars.append(complex_row(next(lines)))
    chars = np.array(chars)
    raw = projector_error(product, factor, chars)
    conjugate = projector_error(product, factor, chars.conj())
    counts['cases'] += 1
    counts['raw_failed'] += raw > 2e-4
    counts['conjugate_failed'] += conjugate > 2e-4
    counts['both_failed'] += min(raw, conjugate) > 2e-4
    if max(raw, conjugate) > 2e-4:
        print(json.dumps(dict(sg=int(sg), spin=int(spin), point=point,
                              raw_error=raw, conjugate_error=conjugate)), flush=True)
print(json.dumps(counts))
