"""Reconstruct exported atomic Bloch characters without CP2K induction caches.

Inputs are completed Gaussian .little_group files. Only NumPy is required.
Every onsite column and every exported reciprocal point is checked. This
oracle supports ordinary and grey groups, with unitary orbit representatives.
It neither guesses irrep names nor merges equal sampled band signatures.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np


def complex_array(values):
    values = np.array(values, dtype=float)
    assert len(values) % 2 == 0
    return values[::2] + 1j*values[1::2]


def read(path):
    records, points = {}, []
    complete = False
    for line in path.read_text().splitlines():
        row = line.split()
        if not row:
            continue
        tag, values = row[0], row[1:]
        if tag == "KPOINT":
            points.append(np.array(values, float))
        if tag == "END_ATOMIC_SIGNATURES":
            complete = True
        if not tag.startswith("ATOMIC_"):
            continue
        records.setdefault(tag, []).append(values)
    assert complete, f"Incomplete atomic export: {path}"
    return records, points


def reference_match(s, spinful, n, local_rows, operations, local_spin, table, references):
    from prepare_wyckoff_characters import grey_table, orient, spin_lift

    native = np.array([operations[s, g][2] for g in local_rows])
    choices = []
    for (group, spin), ref in references.items():
        if spin != spinful:
            continue
        q = orient(ref["cart"], native)
        if q is not None:
            choices.append((group, ref, q))
    assert len(choices) == 1
    group, ref, q = choices[0]
    grey = n != len(local_rows)
    if grey:
        assert n == 2*len(local_rows)
    reference, names = grey_table(ref, bool(spinful)) if grey else (ref["characters"], ref["labels"])
    assert table.shape == reference.shape
    v = spin_lift(q)
    aligned, error = [], 0.0
    for g in local_rows:
        distance = np.max(np.abs(q@ref["cart"]@q.T-operations[s, g][2]), axis=(1, 2))
        index = int(np.argmin(distance))
        assert distance[index] < 1e-8
        phase = 1
        if spinful:
            rotated = v@ref["spin"][index]@v.conj().T
            phase = np.vdot(rotated, local_spin[s, g])/2
            assert abs(abs(phase)-1) < 2e-4
            phase /= abs(phase)
            error = max(error, float(np.max(np.abs(local_spin[s, g]-phase*rotated))))
        aligned.append(phase*reference[index])
    aligned = np.array(aligned)
    dimensions = np.rint(table[0].real)
    assert np.all(dimensions > 0)
    distances = np.max(np.abs(table[:, :, None]-aligned[:, None, :]), axis=0)/dimensions[:, None]
    correspondence = np.argmin(distances, axis=1)
    assert len(set(correspondence)) == len(names)
    assert np.all(np.sum(distances < 2e-4, axis=1) == 1)
    error = max(error, float(np.max(np.min(distances, axis=1))))
    assert error < 2e-4
    return dict(site=s, point_group_reference=group, spinful=bool(spinful), grey=grey,
                reference_to_export_axes=q.tolist(), columns=[names[i] for i in correspondence],
                maximum_reference_residual=error,
                convention="declared point-group axes, not canonical Bilbao onsite/EBR labels")


def check(records, points, references=None):
    def rows(tag):
        return records[tag]

    def close(a, b, message):
        error = float(np.max(np.abs(np.asarray(a)-b)))
        assert error < 1e-8, (message, error)
        return error

    cell = np.array(rows("ATOMIC_REFERENCE_CELL")[0], float).reshape(3, 3, order="F")
    assert len(rows("ATOMIC_SIGNATURE_SIZES")) == 1
    np_, nr, nc = map(int, rows("ATOMIC_SIGNATURE_SIZES")[0])
    assert len(points) == np_
    rotation, tau, anti, spin = {}, {}, {}, {}
    for row in rows("ATOMIC_REFERENCE_OPERATION"):
        g, flag = map(int, row[:2])
        assert g not in rotation
        rotation[g] = np.array(row[2:11], int).reshape(3, 3)
        tau[g], anti[g] = np.array(row[11:], float), bool(flag)
    for row in rows("ATOMIC_REFERENCE_SPIN"):
        g = int(row[0])
        assert g not in spin
        spin[g] = complex_array(row[1:]).reshape(2, 2, order="F")
    assert rotation.keys() == spin.keys()

    def find_operation(r, t, flag=False):
        matches = [g for g in rotation if anti[g] == flag and np.array_equal(rotation[g], r)
                   and np.max(np.abs(t-tau[g]-np.rint(t-tau[g]))) < 1e-8]
        assert len(matches) == 1, matches
        return matches[0]

    columns = {}
    for row in rows("ATOMIC_SIGNATURE_COLUMN"):
        j, s, local, rank = map(int, row[:4])
        assert j not in columns
        columns[j] = (s, local, rank, np.array(row[4:], float))
    assert sorted(columns) == list(range(1, nc+1))
    sizes = {int(row[0]): list(map(int, row[1:])) for row in rows("ATOMIC_ONSITE_SIZES")}
    assert set(sizes) == {row[0] for row in columns.values()}
    operations, local_spin, factors, characters = {}, {}, {}, {}
    for row in rows("ATOMIC_ONSITE_OPERATION"):
        s, g, original, flag = map(int, row[:4])
        assert (s, g) not in operations
        operations[s, g] = (original, flag, np.array(row[4:13], float).reshape(3, 3), np.array(row[13:], float))
    for row in rows("ATOMIC_ONSITE_SPIN"):
        key = tuple(map(int, row[:2]))
        assert key not in local_spin
        local_spin[key] = complex_array(row[2:]).reshape(2, 2, order="F")
    for row in rows("ATOMIC_ONSITE_PRODUCT_ROW"):
        s, g = map(int, row[:2])
        assert len(row[2:]) == sizes[s][0]
        for h, value in enumerate(map(int, row[2:]), 1):
            assert 1 <= abs(value) <= sizes[s][0]
            key = (s, g, h)
            assert key not in factors
            factors[key] = abs(value), (1 if value > 0 else -1)
    for row in rows("ATOMIC_ONSITE_CHARACTER"):
        key = tuple(map(int, row[:3]))
        assert key not in characters
        characters[key] = complex_array(row[3:])[0]

    matrix = {(int(a), int(b)): int(c) for a, b, c in rows("ATOMIC_SIGNATURE_ENTRY")}
    assert len(matrix) == len(rows("ATOMIC_SIGNATURE_ENTRY")) == nr*nc
    signature_rows = [list(map(int, row)) for row in rows("ATOMIC_SIGNATURE_ROW")]
    assert [row[0] for row in signature_rows] == list(range(1, nr+1))
    reciprocal_characters = {(int(p), int(j), int(g)): complex(float(a), float(b))
                            for p, j, g, a, b in rows("ATOMIC_REFERENCE_CHARACTER")}
    coreps = {(int(p), int(j), int(i)): int(n) for p, j, i, n in records.get("ATOMIC_REFERENCE_COREP", [])}
    maximum, entries, site_columns = 0.0, 0, 0
    reference_matches = []
    for s, (n, nu, ni, spinful) in sizes.items():
        col = sorted([j for j, c in columns.items() if c[0] == s], key=lambda j: columns[j][1])
        assert [columns[j][1] for j in col] == list(range(1, ni+1))
        seed = np.mod(columns[col[0]][3], 1)
        onsite = {operations[s, g][0]: g for g in range(1, n+1)}
        assert len(onsite) == n
        assert len([g for g in onsite if not anti[g]]) == nu
        table = np.empty((nu, ni), complex)
        local_rows = [g for g in range(1, n+1) if not operations[s, g][1]]
        for i, g in enumerate(local_rows):
            table[i] = [characters[s, j, g] for j in range(1, ni+1)]
        if references is not None:
            reference_matches.append(reference_match(s, spinful, n, local_rows, operations, local_spin, table, references))
        for g in rotation:
            displacement = rotation[g]@seed+tau[g]-seed
            assert (g in onsite) == (np.max(np.abs(displacement-np.rint(displacement))) < 1e-8)
        for original, g in onsite.items():
            _, flag, cart, shift = operations[s, g]
            assert flag == anti[original]
            close(cart, cell@rotation[original]@np.linalg.inv(cell), "Cartesian onsite frame")
            close(shift, rotation[original]@seed+tau[original]-seed, "Onsite lattice shift")
            close(local_spin[s, g], spin[original], "Full/onsite spin lift")
            for h in range(1, n+1):
                q, factor = factors[s, g, h]
                right = local_spin[s, h].conj() if flag else local_spin[s, h]
                if spinful:
                    close(local_spin[s, g]@right, factor*local_spin[s, q], "Onsite spin product")
                else:
                    close(factor, 1, "Scalar onsite factor")
                close(operations[s, g][2]@operations[s, h][2], operations[s, q][2], "Onsite spatial product")

        orbit, representatives = [], []
        for g in rotation:
            if anti[g]:
                continue
            position = np.mod(rotation[g]@seed+tau[g], 1)
            if any(np.max(np.abs(position-q-np.rint(position-q))) < 1e-8 for q in orbit):
                continue
            orbit.append(position)
            representatives.append(g)
        close(table[0]*len(orbit), [columns[j][2] for j in col], "Band ranks")
        terms = {}
        for g in rotation:
            if anti[g]:
                continue
            terms[g] = []
            for position, t in zip(orbit, representatives):
                displacement = rotation[g]@position+tau[g]-position
                if np.max(np.abs(displacement-np.rint(displacement))) > 1e-8:
                    continue
                inverse = np.rint(np.linalg.inv(rotation[t])).astype(int)
                h = find_operation(inverse@rotation[g]@rotation[t],
                                   inverse@(tau[g]+rotation[g]@tau[t]-tau[t]))
                local = onsite[h]
                transported = spin[t].conj().T@spin[g]@spin[t]
                coefficient = np.vdot(spin[h], transported)/2
                close(transported, coefficient*spin[h], "Orbit spin intertwiner")
                close(abs(coefficient), 1, "Orbit spin phase")
                if not spinful:
                    coefficient = 1
                chars = np.array([characters[s, j, local] for j in range(1, ni+1)])
                terms[g].append((np.rint(displacement), coefficient*chars))
        for p, k in enumerate(points, 1):
            active = [g for g in rotation if not anti[g]
                      and np.max(np.abs(rotation[g].T@k-k-np.rint(rotation[g].T@k-k))) < 1e-8]
            present = {g for pp, _, g in reciprocal_characters if pp == p}
            assert set(active) == present
            local_irreps = sorted({j for pp, j, _ in reciprocal_characters if pp == p})
            chi = np.array([[reciprocal_characters[p, j, g] for j in local_irreps] for g in active])
            rs = [row for row in signature_rows if row[1] == p]
            assert [row[2] for row in rs] == list(range(1, len(rs)+1))
            counts = np.array([[matrix[row[0], j] for j in col] for row in rs])
            if any(pp == p for pp, _, _ in coreps):
                restriction = np.array([[coreps[p, j, i] for j in range(1, len(rs)+1)] for i in local_irreps])
                counts = restriction@counts
            expected = chi@counts
            for i, g in enumerate(active):
                induced = sum((np.exp(-2j*np.pi*k@shift)*chars for shift, chars in terms[g]),
                              np.zeros(ni, complex))
                maximum = max(maximum, close(induced, expected[i], "Induced Bloch character"))
                entries += ni
        site_columns += ni
    assert len(operations) == sum(v[0] for v in sizes.values())
    assert len(local_spin) == len(operations)
    assert len(factors) == sum(v[0]**2 for v in sizes.values())
    assert len(characters) == sum(v[1]*v[2] for v in sizes.values())
    return dict(sites=len(sizes), columns=site_columns, points=len(points), entries=entries,
                maximum_error=maximum, reference_matches=reference_matches)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", type=Path, nargs="+")
    parser.add_argument("--negative", action="store_true")
    parser.add_argument("--onsite-reference", type=Path)
    args = parser.parse_args()
    references, reference_sha = None, None
    if args.onsite_reference:
        from prepare_wyckoff_characters import characters
        references = characters(args.onsite_reference)
        reference_sha = hashlib.sha256(args.onsite_reference.read_bytes()).hexdigest()
    for path in args.files:
        records, points = read(path)
        result = check(records, points, references)
        rejected = []
        if args.negative:
            for tag, index in [("ATOMIC_ONSITE_CHARACTER", 3), ("ATOMIC_ONSITE_OPERATION", 13), ("ATOMIC_REFERENCE_SPIN", 1)]:
                wrong = copy.deepcopy(records)
                wrong[tag][0][index] = str(float(wrong[tag][0][index])+1)
                try:
                    check(wrong, points)
                except AssertionError:
                    rejected.append(tag)
                else:
                    raise AssertionError(f"Corrupted {tag} accepted")
            wrong = copy.deepcopy(records)
            wrong["ATOMIC_ONSITE_PRODUCT_ROW"][0][2] = str(-int(wrong["ATOMIC_ONSITE_PRODUCT_ROW"][0][2]))
            try:
                check(wrong, points)
            except AssertionError:
                rejected.append("ATOMIC_ONSITE_PRODUCT_ROW")
            else:
                raise AssertionError("Corrupted product sign accepted")
        print(json.dumps(dict(file=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                              **result, rejected_corruptions=rejected, onsite_reference_sha256=reference_sha)), flush=True)


if __name__ == "__main__":
    main()
