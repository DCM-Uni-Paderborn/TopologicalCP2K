"""Build induced-character references from labelled irreptables EBR vectors.

Uses the separately validated, convention-aware little-group character fixture
and independently archived affine Wyckoff families. No native induction code is
used. EBR labels remain metadata: equal sampled characters can be ambiguous.
"""
import argparse
import hashlib
import importlib.metadata
import json
import re
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

from prepare_wyckoff_characters import complex_values


def character_records(path):
    rows = iter(path.read_text().splitlines())
    records = defaultdict(list)
    for _ in range(int(next(rows))):
        sg, spin, n, ni, name = next(rows).split()
        sg, spin, n, ni = map(int, (sg, spin, n, ni))
        values = np.fromstring(next(rows), sep=" ")
        cell, k = values[:9].reshape(3, 3, order="F"), values[9:]
        rotations, translations, lifts = [], [], []
        for _ in range(n):
            line = next(rows)
            geometry = np.fromstring(line.split("(")[0], sep=" ")
            rotations.append(geometry[:9].reshape(3, 3, order="F").astype(int))
            translations.append(geometry[9:])
            lifts.append(complex_values(line).reshape(2, 2, order="F"))
        for _ in range(n):
            next(rows)
            next(rows)
        labels, chars = [], []
        for _ in range(ni):
            labels.append(next(rows))
            chars.append(complex_values(next(rows)))
        records[sg, spin].append(dict(name=name, cell=cell, k=k, rotations=np.array(rotations),
                                     tau=np.array(translations), spin=np.array(lifts),
                                     labels=labels, characters=np.array(chars).T))
    assert next(rows, None) is None
    return records


def affine_records(path):
    rows = iter(path.read_text().splitlines())
    records = {}
    for _ in range(int(next(rows))):
        sg, n, nf = map(int, next(rows).split())
        centers = []
        for _ in range(n):
            row = np.fromstring(next(rows), sep=" ")
            if np.array_equal(row[:9].reshape(3, 3), np.eye(3)):
                centers.append(Matrix(np.rint(24*row[9:]).astype(int)))
        primitive = np.array(hermite_normal_form(Matrix.hstack(24*Matrix.eye(3), *centers))).astype(float)/24
        families = {}
        for _ in range(nf):
            row = next(rows).split()
            values = np.array(list(map(float, row[2:])))
            families[row[0]] = (int(row[1]), values[:3], values[3:].reshape(3, 3, order="F"))
        records[sg] = (primitive, len(centers), families)
    assert next(rows, None) is None
    return records


def fmt(values):
    return " ".join(f"{v:.17g}" for v in np.asarray(values).flatten(order="F"))


def cfmt(values):
    return " ".join(f"({v.real:.17g},{v.imag:.17g})" for v in np.asarray(values).flatten(order="F"))


def supplement_type_one(root, sg, spin, full, primitive, tables, missing):
    """Read explicit missing points, checking the shared labelled characters first."""
    def coordinate(value):
        rational = Fraction(value).limit_denominator(48)
        assert abs(float(rational)-float(value)) < 5.1e-6
        return float(rational)

    selected = []
    for path in root.glob(f"irreps-SG={sg}.*-{'spin' if spin else 'scal'}.dat"):
        rows = iter(path.read_text().splitlines())
        n = None
        for line in rows:
            if line.strip().startswith("nsym="):
                n = int(line.split("=")[1])
            if line.strip() == "symmetries=":
                break
        operations = [next(rows).split() for _ in range(n)]
        if not all(row[-1] == "1" for row in operations):
            continue
        assert n == len(full["rotations"])
        selected.append((path, operations, list(rows)))
    assert len(selected) == 1, (sg, spin, "type-I table", len(selected))
    path, operations, rows = selected[0]
    inverse = np.linalg.inv(primitive)
    mapping, phases, shifts = [], [], []
    for row in operations:
        r = inverse@np.array(list(map(int, row[:9]))).reshape(3, 3)@primitive
        matches = np.flatnonzero(np.max(np.abs(full["rotations"]-r), axis=(1, 2)) < 1e-10)
        assert len(matches) == 1
        g = int(matches[0])
        delta = inverse@np.array(list(map(coordinate, row[9:12])))-full["tau"][g]
        assert np.max(np.abs(delta-np.rint(delta))) < 1e-10
        phase = 1
        if spin:
            s = (np.array(list(map(float, row[12:16])))*np.exp(1j*np.pi*np.array(list(map(float, row[16:20]))))).reshape(2, 2)
            phase = np.vdot(s, full["spin"][g])/2
            assert abs(abs(phase)-1) < 2e-4, (sg, "spin-frame")
            phase /= abs(phase)
            assert np.max(np.abs(full["spin"][g]-phase*s)) < 2e-4
        mapping.append(g)
        phases.append(phase)
        shifts.append(np.rint(delta))
    assert len(set(mapping)) == len(mapping)
    additional = []
    for line in rows:
        if not line.strip():
            continue
        if line.strip().startswith("kpoint "):
            header, coordinates, indices = line.split(":")
            indices = np.array(list(map(int, indices.split())))-1
            active = np.array(mapping)[indices]
            k = -primitive.T@np.array(list(map(coordinate, coordinates.split())))
            phase = np.array(phases)[indices]*np.exp(2j*np.pi*np.array(shifts)[indices]@k)
            additional.append(dict(name=header.split()[1], cell=full["cell"], k=k,
                rotations=full["rotations"][active], tau=full["tau"][active], spin=full["spin"][active],
                labels=[], characters=[]))
        else:
            row = line.split()
            values = np.array(list(map(float, row[2:])))
            count = len(indices)
            assert len(values) in (count, 2*count)
            value = values[:count].astype(complex)
            if len(values) == 2*count:
                value *= np.exp(1j*np.pi*values[count:])
            additional[-1]["labels"].append(row[0])
            additional[-1]["characters"].append(phase*value)
    largest = 0.0
    result, common = [], 0
    for table in additional:
        table["characters"] = np.array(table["characters"]).T
        previous = [t for t in tables if t["name"] == table["name"]]
        if previous:
            old, = previous
            assert np.max(np.abs(old["k"]-table["k"])) < 1e-9
            assert set(old["labels"]) == set(table["labels"])
            indices = [int(np.flatnonzero(np.all(table["rotations"] == r, axis=(1, 2)))[0]) for r in old["rotations"]]
            columns = [table["labels"].index(label) for label in old["labels"]]
            error = np.max(np.abs(old["characters"]-table["characters"][np.ix_(indices, columns)])/old["characters"][0].real)
            largest = max(largest, float(error))
            assert error < 2e-4, (sg, spin, table["name"], error)
            common += 1
        elif set(table["labels"]) & missing:
            assert set(table["labels"]) <= missing
            result.append(table)
    assert {l for t in result for l in t["labels"]} == missing
    assert common == len(tables)
    return result, dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       shared_points=common, shared_character_residual=largest,
                       added_points=[t["name"] for t in result])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("characters", type=Path)
    parser.add_argument("affine", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("groups", type=int, nargs="*")
    args = parser.parse_args()
    assert importlib.metadata.version("irreptables") == "3.1.0"
    root = Path(importlib.metadata.distribution("irreptables").locate_file("irreptables/data/ebrs"))
    chars = character_records(args.characters)
    affine = affine_records(args.affine)
    groups = args.groups or list(range(1, 231))
    assert len(set(groups)) == len(groups) and all(1 <= sg <= 230 for sg in groups)
    output, sources, cases, supplements = [], [], [], []
    for sg in groups:
        path = root/f"{sg}_ebrs.json"
        data = path.read_bytes()
        source = json.loads(data)
        sources.append(dict(file=path.name, sha256=hashlib.sha256(data).hexdigest()))
        p, centering, families = affine[sg]
        for spin in (0, 1):
            tables = chars[sg, spin]
            full, = [r for r in tables if np.max(np.abs(r["k"])) < 1e-12]
            basis = source["double" if spin else "single"]["basis"]
            ebrs = source["double" if spin else "single"]["ebrs"]
            basis_labels = basis["irrep_labels"]
            assert len(set(basis_labels)) == len(basis_labels)
            dimensions = dict(zip(basis_labels, basis["degeneracies"]))
            missing = set(basis_labels)-{l for t in tables for l in t["labels"]}
            if missing:
                extra, provenance = supplement_type_one(root.parent/"correptables", sg, spin, full, p, tables, missing)
                tables = tables+extra
                supplements.append(dict(sg=sg, spin=spin, **provenance))
            covered = []
            for table in tables:
                assert np.max(np.abs(table["cell"]-full["cell"])) < 1e-12
                for label, dim in zip(table["labels"], table["characters"][0].real):
                    assert label in dimensions and abs(dim-dimensions[label]) < 1e-5, (sg, spin, label)
                    covered.append(label)
            assert set(covered) == set(basis_labels) and len(covered) == len(basis_labels), (sg, spin, covered, basis_labels)
            sites = defaultdict(list)
            for ebr in ebrs:
                multiplicity, label = re.match(r"(\d+)([a-z]+)\(", ebr["wyckoff_position"]).groups()
                assert families[label][0] == int(multiplicity)
                sites[label].append(ebr)
            for site, references in sites.items():
                mult, origin, directions = families[site]
                assert mult % centering == 0
                seed = np.linalg.solve(p, origin+directions@(np.sqrt([2, 3, 5])/7))
                n, nr, nq = len(full["rotations"]), len(references), len(tables)
                vectors = np.array([r["vector"] for r in references]).T
                assert vectors.shape == (len(basis_labels), nr)
                assert np.all(np.isfinite(vectors)) and np.all(vectors >= 0)
                assert np.array_equal(vectors, np.rint(vectors))
                names = [r["ebr_name"] for r in references]
                assert len(set(names)) == nr
                output.append(f"{sg} {spin} {n} {nr} {nq} {mult//centering} {site}")
                output.append(fmt(full["cell"])+" "+fmt(seed))
                for r, t, s in zip(full["rotations"], full["tau"], full["spin"]):
                    output.append(fmt(r)+" "+fmt(t)+" "+cfmt(s))
                for name in names:
                    output.append(name.replace("\u2191", "_up_"))
                band_dims = np.array([int(re.search(r"\((\d+)\)$", name)[1]) for name in names])
                for table in tables:
                    subset = [basis_labels.index(label) for label in table["labels"]]
                    expected = table["characters"]@vectors[subset]
                    assert np.max(np.abs(expected[0]-band_dims)) < 1e-4
                    output.append(fmt(table["k"])+f" {len(table['rotations'])}")
                    for r, t, s, row in zip(table["rotations"], table["tau"], table["spin"], expected):
                        matches = np.flatnonzero(np.all(full["rotations"] == r, axis=(1, 2)))
                        assert len(matches) == 1
                        g = int(matches[0])
                        assert np.max(np.abs(full["tau"][g]-t)) < 1e-12
                        assert np.max(np.abs(full["spin"][g]-s)) < 1e-12
                        output.append(f"{g+1} "+cfmt(row))
                cases.append(dict(sg=sg, spin=spin, site=site, columns=nr, points=nq,
                                  entries=sum(len(t["rotations"])*nr for t in tables), labels=names))
                print(json.dumps({k: v for k, v in cases[-1].items() if k != "labels"}), flush=True)
    args.output.write_text(str(len(cases))+"\n"+"\n".join(output)+"\n")
    args.manifest.write_text(json.dumps(dict(package="irreptables", version="3.1.0", sources=sources,
        characters_sha256=hashlib.sha256(args.characters.read_bytes()).hexdigest(),
        affine_sha256=hashlib.sha256(args.affine.read_bytes()).hexdigest(),
        convention="ordinary groups; k_native=-P^T*k_table; raw table characters and lifts",
        supplements=supplements, cases=cases), indent=2)+"\n")


if __name__ == "__main__":
    main()
