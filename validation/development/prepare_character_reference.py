"""Prepare character fixtures from irreptables 3.1.0, retaining source rounding.

Usage: affine-reference.dat output.dat sources.jsonl [SG ...]
The affine fixture supplies the conventional centering lattice, not characters.
No irrep implementation is imported; only the installed reference data is read.
"""
import argparse
import hashlib
import importlib.metadata
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

assert importlib.metadata.version("irreptables") == "3.1.0"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("affine", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("sources", type=Path)
parser.add_argument("groups", type=int, nargs="*")
parser.add_argument("--native-k-sign", type=int, choices=(-1, 1), default=-1,
                    help="Use +1 only to reproduce the unconverted-convention audit")
args = parser.parse_args()
dist = importlib.metadata.distribution("irreptables")
root = Path(dist.locate_file("irreptables/data/tables"))
groups = args.groups or list(range(1, 231))
assert len(groups) == len(set(groups)) and all(1 <= sg <= 230 for sg in groups)
sigma = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
source_records = []


def rational(value):
    result = Fraction(value).limit_denominator(48)
    # Coordinates are printed to five decimal places; do not round characters.
    assert abs(float(result)-float(value)) < 5.1e-6, value
    return result


def parse_table(path):
    data = path.read_bytes()
    source_records.append(dict(package="irreptables", version="3.1.0", file=path.name,
                               sha256=hashlib.sha256(data).hexdigest(),
                               coordinate_rounding_bound=5.1e-6,
                               coordinate_max_denominator=48,
                               native_k_convention=f"{args.native_k_sign:+d} * P^T k_table; raw characters and spin lifts retained"))
    lines = iter(data.decode().splitlines())
    n = None
    for line in lines:
        if line.strip().startswith("nsym="):
            n = int(line.split("=")[1])
        if line.strip() == "symmetries=":
            break
    assert n is not None
    ops = []
    for _ in range(n):
        row = next(lines).split()
        r = np.array(list(map(int, row[:9]))).reshape(3, 3)
        t = np.array([float(rational(v)) for v in row[9:12]])
        if len(row) == 20:
            s = (np.array(list(map(float, row[12:16]))) *
                 np.exp(1j*np.pi*np.array(list(map(float, row[16:20]))))).reshape(2, 2)
        else:
            assert len(row) == 12
            s = np.eye(2, dtype=complex)
        ops.append((r, t, s))
    points = []
    for line in lines:
        if not line.strip():
            continue
        if line.strip().startswith("kpoint "):
            header, coords, indices = line.split(":")
            points.append(dict(name=header.split()[1], k=np.array([float(rational(v)) for v in coords.split()]),
                               indices=[int(v)-1 for v in indices.split()], labels=[], characters=[]))
        else:
            row = line.split()
            point = points[-1]
            count = len(point["indices"])
            label, dim = row[0], int(row[1])
            values = np.array(list(map(float, row[2:])))
            assert len(values) in (count, 2*count)
            chars = values[:count].astype(complex)
            if len(values) == 2*count:
                chars *= np.exp(1j*np.pi*values[count:])
            assert abs(chars[0]-dim) < 1e-6
            point["labels"].append(label)
            point["characters"].append(chars)
    return ops, points


def embedding(rotations, spins, spinful):
    metric = sum(r.T@r for r in rotations)/len(rotations)
    cell = np.linalg.cholesky(metric).T
    if not spinful:
        return cell, 0.0
    cart = [cell@r@np.linalg.inv(cell)*round(np.linalg.det(r)) for r in rotations]
    polar_spins = [np.linalg.svd(s)[0]@np.linalg.svd(s)[2] for s in spins]
    target = [np.array([[np.trace(a@s@b@s.conj().T).real/2 for b in sigma] for a in sigma])
              for s in polar_spins]
    equations = np.vstack([np.kron(q.T, np.eye(3))-np.kron(np.eye(3), p) for q, p in zip(cart, target)])
    _, singular, vh = np.linalg.svd(equations, full_matrices=False)
    null = vh[singular < 2e-4]
    assert len(null), ("No reference spin/space intertwiner", singular)
    best = None
    for trial in range(1, 33):
        indices = np.arange(1, len(null)+1)
        weights = np.sin(indices**2*trial*np.sqrt(2)) + np.cos(indices*trial**2*np.sqrt(3))
        raw = (weights@null).reshape(3, 3, order="F")
        u, sv, v = np.linalg.svd(raw)
        if sv[-1] < 1e-4:
            continue
        orient = u@v
        if np.linalg.det(orient) < 0:
            orient = -orient
        error = max(np.max(np.abs(orient@q@orient.T-p)) for q, p in zip(cart, target))
        if best is None or error < best[0]:
            best = (error, orient)
    assert best is not None and best[0] < 2e-4, best
    return best[1]@cell, best[0]


fixture = iter(args.affine.read_text().splitlines())
primitive = {}
for _ in range(int(next(fixture))):
    sg, n, nf = map(int, next(fixture).split())
    translations = []
    for _ in range(n):
        row = next(fixture).split()
        if list(map(int, row[:9])) == [1,0,0,0,1,0,0,0,1]:
            translations.append([round(24*float(v)) for v in row[9:]])
    basis = Matrix.hstack(24*Matrix.eye(3), *[Matrix(v) for v in translations])
    primitive[sg] = np.array(hermite_normal_form(basis)).astype(float)/24
    for _ in range(nf):
        next(fixture)


def fmt(values):
    return " ".join(f"{v:.17g}" for v in np.asarray(values).flatten(order="F"))


def cfmt(values):
    return " ".join(f"({v.real:.17g},{v.imag:.17g})" for v in np.asarray(values).flatten(order="F"))


output = []
count = 0
for sg in groups:
    p = primitive[sg]
    inv = np.linalg.inv(p)
    for spin in (0, 1):
        ops, points = parse_table(root / f"irreps-SG={sg}-{'spin' if spin else 'scal'}.dat")
        rotations = []
        for r, _, _ in ops:
            transformed = inv@r@p
            assert np.max(np.abs(transformed-np.rint(transformed))) < 1e-10
            rotations.append(np.rint(transformed).astype(int))
        tau = [inv@t for _, t, _ in ops]
        spins = [s for _, _, s in ops]
        cell, embed_error = embedding(rotations, spins, bool(spin))
        for point in points:
            indices = point["indices"]
            r = [rotations[i] for i in indices]
            t = [tau[i] for i in indices]
            s = [spins[i] for i in indices]
            # The tabulated characters realize exp(+2*pi*i*k_table*T).
            # CP2K uses exp(-2*pi*i*k*T): compare at -k_table, not by relabeling
            # columns or selectively conjugating characters until they match.
            k = args.native_k_sign*p.T@point["k"]
            n, ni = len(indices), len(point["labels"])
            assert indices[0] == 0 and ni > 0
            assert len(point["labels"]) == len(set(point["labels"]))
            products = np.empty((n, n), dtype=int)
            factors = np.empty((n, n), dtype=complex)
            rotation_indices = {tuple(op.flatten()): j for j, op in enumerate(r)}
            assert len(rotation_indices) == n
            for g in range(n):
                assert np.max(np.abs(r[g].T@k-k-np.rint(r[g].T@k-k))) < 1e-9
                for h in range(n):
                    q = rotation_indices[tuple((r[g]@r[h]).flatten())]
                    delta = t[g]+r[g]@t[h]-t[q]
                    assert np.max(np.abs(delta-np.rint(delta))) < 1e-9
                    products[g, h] = q+1
                    sign = 1
                    if spin:
                        composed = s[g]@s[h]
                        sign = 1 if np.vdot(s[q], composed).real >= 0 else -1
                        assert np.max(np.abs(composed-sign*s[q])) < 2e-4
                    factors[g, h] = sign*np.exp(-2j*np.pi*k@np.rint(delta))
            output.append(f"{sg} {spin} {n} {ni} {point['name']}")
            output.append(fmt(cell)+" "+fmt(k))
            for g in range(n):
                output.append(fmt(r[g])+" "+fmt(t[g])+" "+cfmt(s[g]))
            for h in range(n):
                output.append(" ".join(map(str, products[:, h])))
                output.append(cfmt(factors[:, h]))
            for label, chars in zip(point["labels"], point["characters"]):
                output.append(label)
                output.append(cfmt(chars))
            count += 1
        print(json.dumps(dict(sg=sg, spin=spin, points=len(points), spin_embedding_residual=embed_error)), flush=True)
args.output.write_text(str(count)+"\n"+"\n".join(output)+"\n")
args.sources.write_text("".join(json.dumps(row)+"\n" for row in source_records))
print(json.dumps(dict(cases=count, groups=len(groups))))
