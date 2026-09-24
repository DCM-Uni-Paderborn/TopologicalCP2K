"""Independent point-group character references at actual affine Wyckoff sites.

The axis correspondence is explicit, not a Bilbao site-label assignment.
Inputs: convention-aware character fixture, affine fixture, output, manifest.
Only reference data, NumPy/SciPy and exact integer centering arithmetic are used.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

POINT_GROUPS = [1, 2, 3, 6, 10, 16, 25, 47, 75, 81, 83, 89, 99, 111, 123,
                143, 147, 149, 156, 162, 168, 174, 175, 177, 183, 187, 191,
                195, 200, 207, 215, 221]
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
POINTS = np.array([[0, 0, 0], [.5, 0, 0], [0, .5, 0], [0, 0, .5],
                   [.25, .25, .25], [.173, .217, .319]])


def spin_lift(cart):
    proper = cart*round(np.linalg.det(cart))
    x, y, z, w = Rotation.from_matrix(proper).as_quat()
    return w*np.eye(2)-1j*np.einsum("i,ijk->jk", [x, y, z], SIGMA)


def induce(cell, rotations, tau, seed, ref, orientation, table, spinful):
    """Direct site-orbit trace with explicit physical spin intertwiners."""
    inverse = np.linalg.inv(cell)
    cart = cell@rotations@inverse
    lifts = np.array([spin_lift(r) for r in cart])
    v = spin_lift(orientation)
    site_lifts = v@ref["spin"]@v.conj().T
    positions, representatives = [], []
    for g, (r, t) in enumerate(zip(rotations, tau)):
        position = np.mod(r@seed+t, 1)
        if any(np.max(np.abs(position-old-np.rint(position-old))) < 1e-9 for old in positions):
            continue
        positions.append(position)
        representatives.append(g)
    terms = [[] for _ in rotations]
    for g, (r, t) in enumerate(zip(rotations, tau)):
        for position, representative in zip(positions, representatives):
            delta = r@position+t-position
            if np.max(np.abs(delta-np.rint(delta))) > 1e-9:
                continue
            a = cart[representative]
            h = a.T@cart[g]@a
            distance = np.max(np.abs(orientation@ref["cart"]@orientation.T-h), axis=(1, 2))
            index = int(np.argmin(distance))
            assert distance[index] < 1e-8
            coefficient = 1
            if spinful:
                transported = lifts[representative].conj().T@lifts[g]@lifts[representative]
                coefficient = np.vdot(site_lifts[index], transported)/2
                assert abs(abs(coefficient)-1) < 2e-4
                coefficient /= abs(coefficient)
            terms[g].append((np.rint(delta), coefficient*table[index]))
    result = []
    for k in POINTS:
        active, values = [], []
        for g, r in enumerate(rotations):
            difference = r.T@k-k
            if np.max(np.abs(difference-np.rint(difference))) > 1e-9:
                continue
            active.append(g)
            value = sum((np.exp(-2j*np.pi*k@delta)*chars for delta, chars in terms[g]),
                        np.zeros(table.shape[1], dtype=complex))
            values.append(value)
        result.append((active, values))
    return lifts, result


def complex_values(line):
    return np.array([complex(float(a), float(b)) for a, b in re.findall(r"\(([^,]+),([^\)]+)\)", line)])


def characters(path):
    rows = iter(path.read_text().splitlines())
    result = {}
    for _ in range(int(next(rows))):
        sg, spin, n, ni, point = next(rows).split()
        sg, spin, n, ni = map(int, (sg, spin, n, ni))
        values = np.fromstring(next(rows), sep=" ")
        cell, k = values[:9].reshape(3, 3, order="F"), values[9:]
        operations, lifts = [], []
        for _ in range(n):
            line = next(rows)
            operations.append(np.fromstring(line.split("(")[0], sep=" ")[:9].reshape(3, 3, order="F"))
            lifts.append(complex_values(line).reshape(2, 2, order="F"))
        products, factors = [], []
        for _ in range(n):
            products.append(np.fromstring(next(rows), sep=" ", dtype=int)-1)
            factors.append(complex_values(next(rows)))
        names, table = [], []
        for _ in range(ni):
            names.append(next(rows))
            table.append(complex_values(next(rows)))
        if sg in POINT_GROUPS and np.max(np.abs(k)) < 1e-12:
            result[sg, spin] = dict(cart=np.array([cell@r@np.linalg.inv(cell) for r in operations]),
                                    spin=np.array(lifts), characters=np.array(table).T, labels=names,
                                    product=np.array(products).T, factor=np.array(factors).T)
    assert len(result) == 64
    return result


def axis_data(rotations):
    axes = []
    for r in rotations:
        det, trace = round(np.linalg.det(r)), round(np.trace(r))
        if abs(trace) == 3:
            continue
        _, _, vh = np.linalg.svd(r-det*np.eye(3))
        axes.append(((det, trace), vh[-1]))
    return axes


def frame(a, b):
    b = b-a*np.dot(a, b)
    if np.linalg.norm(b) < 1e-8:
        return None
    b /= np.linalg.norm(b)
    return np.column_stack((a, b, np.cross(a, b)))


def orient(reference, native):
    if len(reference) != len(native):
        return None
    signature = lambda group: sorted((round(np.linalg.det(r)), round(np.trace(r))) for r in group)
    if signature(reference) != signature(native):
        return None
    axes, target = axis_data(reference), axis_data(native)
    candidates = []
    if not axes:
        candidates.append(np.eye(3))
    else:
        tag, a = axes[0]
        second = next(((t, b) for t, b in axes if abs(np.dot(a, b)) < 1-1e-8), None)
        if second is None:
            f = frame(a, np.eye(3)[np.argmin(np.abs(a))])
            for t, c in target:
                if t == tag:
                    for sign in (1, -1):
                        v = sign*c
                        candidates.append(frame(v, np.eye(3)[np.argmin(np.abs(v))])@f.T)
        else:
            tag2, b = second
            f = frame(a, b)
            for t, c in target:
                if t != tag:
                    continue
                for u, d in target:
                    if u != tag2 or abs(np.dot(c, d)) > 1-1e-8:
                        continue
                    for s in (1, -1):
                        for z in (1, -1):
                            candidates.append(frame(s*c, z*d)@f.T)
    for q in candidates:
        differences = np.max(np.abs((q@reference@q.T)[:, None]-native[None, :]), axis=(2, 3))
        mapping = np.argmin(differences, axis=1)
        if len(set(mapping)) == len(native) and np.max(np.min(differences, axis=1)) < 1e-8:
            return q
    return None


def grey_table(ref, spinful):
    table, used, columns, names = ref["characters"], set(), [], []
    for j in range(table.shape[1]):
        if j in used:
            continue
        errors = np.max(np.abs(table-table[:, j:j+1].conj()), axis=0)
        partner = int(np.argmin(errors))
        assert errors[partner] < 2e-4
        indicator = (1-2*spinful)*np.mean(np.diag(ref["factor"])*table[np.diag(ref["product"]), j])
        value = round(indicator.real)
        assert abs(indicator-value) < 2e-4 and value in (-1, 0, 1), indicator
        if value == 0:
            assert partner != j and partner not in used
            columns.append(table[:, j]+table[:, partner])
            names.append(ref["labels"][j]+"+"+ref["labels"][partner])
            used.update((j, partner))
        else:
            assert partner == j
            columns.append(table[:, j]*(1 if value == 1 else 2))
            names.append(("" if value == 1 else "2*")+ref["labels"][j])
            used.add(j)
    return np.array(columns).T, names


def fmt(values):
    return " ".join(f"{x:.17g}" for x in np.asarray(values).flatten(order="F"))


def cfmt(values):
    return " ".join(f"({x.real:.17g},{x.imag:.17g})" for x in np.asarray(values).flatten(order="F"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("characters", type=Path)
    parser.add_argument("affine", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--groups", type=int, nargs="+")
    args = parser.parse_args()
    refs = characters(args.characters)
    rows = iter(args.affine.read_text().splitlines())
    output, records, matched = [], [], {}
    theta = np.array([[0, 1], [-1, 0]])
    for _ in range(int(next(rows))):
        sg, n, nf = map(int, next(rows).split())
        raw = np.array([np.fromstring(next(rows), sep=" ") for _ in range(n)])
        rotations = raw[:, :9].reshape(n, 3, 3).transpose(0, 2, 1).astype(int)
        tau = raw[:, 9:]
        sites = [next(rows).split() for _ in range(nf)]
        if args.groups and sg not in args.groups:
            continue
        centers = [t for r, t in zip(rotations, tau) if np.array_equal(r, np.eye(3))]
        h = hermite_normal_form(Matrix.hstack(24*Matrix.eye(3), *[Matrix(np.rint(24*t).astype(int)) for t in centers]))
        p = np.array(h).astype(float)/24
        inv = np.linalg.inv(p)
        metric = np.mean([r.T@r for r in rotations], axis=0)
        cell = np.linalg.cholesky(metric).T@p
        primitive, translations, seen = [], [], set()
        for r, t in zip(rotations, tau):
            transformed = inv@r@p
            assert np.max(np.abs(transformed-np.rint(transformed))) < 1e-9
            transformed = np.rint(transformed).astype(int)
            key = tuple(transformed.flat)
            if key not in seen:
                seen.add(key)
                primitive.append(transformed)
                translations.append(np.mod(inv@t, 1))
        rotations, tau = np.array(primitive), np.array(translations)
        for row in sites:
            label, multiplicity = row[0], int(row[1])
            assert multiplicity % len(centers) == 0
            multiplicity //= len(centers)
            values = np.array(list(map(float, row[2:])))
            seed = inv@(values[:3]+values[3:].reshape(3, 3, order="F")@(np.sqrt([2, 3, 5])/7))
            displacement = rotations@seed+tau-seed
            indices = np.flatnonzero(np.max(np.abs(displacement-np.rint(displacement)), axis=1) < 1e-9)
            assert len(indices)*multiplicity == len(rotations), (sg, label)
            native = cell@rotations[indices]@np.linalg.inv(cell)
            for spin in (0, 1):
                key = (spin, tuple(np.round(native, 12).flat))
                if key not in matched:
                    selected = []
                    for group in POINT_GROUPS:
                        ref = refs[group, spin]
                        q = orient(ref["cart"], native)
                        if q is not None:
                            selected.append((group, ref, q))
                    assert len(selected) == 1, (sg, label, spin, [item[0] for item in selected])
                    matched[key] = selected[0]
                group, ref, q = matched[key]
                for grey in (0, 1):
                    chars, names = grey_table(ref, spin) if grey else (ref["characters"], ref["labels"])
                    lifts, induced = induce(cell, rotations, tau, seed, ref, q, chars, bool(spin))
                    nfull, nsite, ni = len(rotations)*(1+grey), len(indices)*(1+grey), len(names)
                    output.append(f"{sg} {label} {spin} {grey} {group} {nfull} {nsite} {ni} {multiplicity}")
                    output.append(fmt(cell)+" "+fmt(seed)+" "+fmt(q))
                    for anti in range(grey+1):
                        for r, t, u in zip(rotations, tau, lifts):
                            output.append(fmt(r)+" "+fmt(t)+" "+cfmt(u@theta if anti else u)+f" {anti}")
                    for anti in range(grey+1):
                        for r, u in zip(ref["cart"], ref["spin"]):
                            output.append(fmt(r)+" "+cfmt(u@theta if anti else u)+f" {anti}")
                    for j, name in enumerate(names):
                        output.append(name)
                        output.append(cfmt(chars[:, j]))
                    output.append(str(len(POINTS)))
                    for k, (active, values) in zip(POINTS, induced):
                        output.append(fmt(k)+f" {len(active)}")
                        for g, value in zip(active, values):
                            output.append(str(g+1)+" "+cfmt(value))
                    records.append(dict(sg=sg, site=label, spin=spin, grey=grey, reference_group=group,
                                        orbit=multiplicity, unitary_stabilizer=len(indices), columns=ni,
                                        labels=names, reciprocal_members=[[int(g)+1 for g in active] for active, _ in induced]))
            print(json.dumps(dict(sg=sg, site=label, orbit=multiplicity)), flush=True)
    args.output.write_text(str(len(records))+"\n"+"\n".join(output)+"\n")
    args.manifest.write_text(json.dumps(dict(records=records, fixture_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
        inputs={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (args.characters, args.affine)}), indent=2)+"\n")
    print(json.dumps(dict(cases=len(records), columns=sum(row["columns"] for row in records))), flush=True)


if __name__ == "__main__":
    main()
