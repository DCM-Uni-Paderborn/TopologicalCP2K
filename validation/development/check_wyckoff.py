"""Independent affine-coordinate comparison against SPGLIB's Wyckoff.csv."""

import ast
import csv
import hashlib
import io
from itertools import product
from pathlib import Path
import sys
import tokenize

import numpy as np
import spglib


def affine(expression):
    tokens = list(tokenize.generate_tokens(io.StringIO(expression).readline))
    explicit = []
    previous = None
    for token in tokens:
        if previous == tokenize.NUMBER and token.type == tokenize.NAME:
            explicit.append((tokenize.OP, "*"))
        explicit.append((token.type, token.string))
        previous = token.type
    tree = ast.parse(tokenize.untokenize(explicit), mode="eval").body

    def evaluate(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return np.array([0.0, 0.0, 0.0, float(node.value)])
        if isinstance(node, ast.Name) and node.id in ("x", "y", "z"):
            return np.eye(4)["xyz".index(node.id)]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return evaluate(node.operand) * (-1 if isinstance(node.op, ast.USub) else 1)
        if isinstance(node, ast.BinOp):
            a, b = evaluate(node.left), evaluate(node.right)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                if not a[:3].any():
                    return a[3] * b
                if not b[:3].any():
                    return b[3] * a
            if isinstance(node.op, ast.Div) and not b[:3].any():
                return a / b[3]
        raise ValueError(f"Not an affine coordinate: {expression}")

    return evaluate(tree)


database_path, output_path = map(Path, sys.argv[1:])
reference = {}
hall = None
with database_path.open() as stream:
    for row in csv.reader(stream, delimiter=":"):
        if row[0] == "end of data":
            break
        if row[0]:
            hall = int(row[0])
            reference[hall] = []
        elif row[2]:
            coordinate = row[5].strip("()")
            transform = np.array([affine(x) for x in coordinate.split(",")])
            reference[hall].append((row[3], int(row[2]), transform))

generated = {}
specializations = {}
for line in output_path.read_text().splitlines():
    if line.startswith("SPECIALIZATIONS "):
        fields = list(map(int, line.split()[1:]))
        hall, j, count = fields[:3]
        assert len(fields[3:]) == count
        specializations[hall, j - 1] = {x - 1 for x in fields[3:]}
    if not line.startswith("WYCKOFF "):
        continue
    fields = line.split()
    hall, dimension, multiplicity = map(int, fields[1:4])
    point = np.array(list(map(float, fields[4:7])))
    normal = np.array(list(map(int, fields[7:16]))).reshape((3, 3), order="F")
    offset = np.array(list(map(float, fields[16:19])))
    generated.setdefault(hall, []).append((dimension, multiplicity, point, normal, offset))

assert set(generated) == set(reference) == set(range(1, 531))
total = 0
worst = 0.0
edges = 0
for hall, entries in reference.items():
    candidates = generated[hall]
    assert len(candidates) == len(entries), (hall, len(candidates), len(entries))
    symmetry = spglib.get_symmetry_from_database(hall)
    rotations, translations = symmetry["rotations"], symmetry["translations"]
    matches = np.zeros((len(entries), len(candidates)), dtype=bool)
    for i, (letter, multiplicity, transform) in enumerate(entries):
        direction = transform[:, :3]
        dimension = np.linalg.matrix_rank(direction)
        images = np.einsum("gij,j->gi", rotations, transform[:, 3]) + translations
        tangents = rotations @ direction
        for j, (dim, mult, point, normal, offset) in enumerate(candidates):
            if dim != dimension or mult != multiplicity:
                continue
            residuals = images @ normal.T - offset
            residuals -= np.rint(residuals)
            errors = np.maximum(np.max(abs(residuals), axis=1), np.max(abs(normal @ tangents), axis=(1, 2)))
            if errors.min() < 1e-8:
                matches[i, j] = True
                worst = max(worst, errors.min())
    assert (matches.sum(axis=0) == 1).all(), (hall, "generated coverage", matches.sum(axis=0))
    assert (matches.sum(axis=1) == 1).all(), (hall, "reference coverage", matches.sum(axis=1))
    for parent, (letter, mult, transform) in enumerate(entries):
        parent_dimension = np.linalg.matrix_rank(transform[:, :3])
        u, s, v = np.linalg.svd(transform[:, :3])
        normal = u[:, parent_dimension:].T
        child_entries = set()
        for child, (child_letter, child_mult, child_transform) in enumerate(entries):
            if np.linalg.matrix_rank(child_transform[:, :3]) >= parent_dimension:
                continue
            # Match affine closures from the independent parameterizations. Lift through the
            # unit cube, because a floating orthonormal normal is not an integer torus normal.
            directions = normal @ (rotations @ child_transform[:, :3])
            good = np.max(abs(directions), axis=(1, 2)) < 1e-8 if normal.size else np.ones(len(rotations), bool)
            images = np.einsum("gij,j->gi", rotations[good], child_transform[:, 3]) + translations[good]
            offsets = images - transform[:, 3]
            # Tabulated conventional-cell coefficients/offsets are small; [-2,2]^3 covers
            # these standard representatives (not an algorithm for arbitrary sheared cells).
            shifts = np.array(list(product(range(-2, 3), repeat=3)))
            errors = np.max(abs((offsets[:, None, :] - shifts) @ normal.T), axis=2, initial=0)
            if np.any(errors < 1e-8):
                child_entries.add(int(matches[child].argmax()))
        generated_parent = int(matches[parent].argmax())
        assert specializations[hall, generated_parent] == child_entries, (hall, letter, "specializations")
        edges += len(child_entries)
    total += len(entries)
    print(f"Hall {hall}: {len(entries)} affine families matched bijectively", flush=True)

print(f"PASS: all 530 Hall settings, {total} Wyckoff families; largest affine residual {worst:.3e}")
print(f"PASS: {edges} directed specialization relations from independent affine parameterizations")
print("SPGLIB operation database:", spglib.__version__)
print("Wyckoff.csv source: spglib/spglib a6b561fb60cdd021a1ac90852c3c2ae14405b2c8")
print("Wyckoff.csv SHA256:", hashlib.sha256(database_path.read_bytes()).hexdigest())
