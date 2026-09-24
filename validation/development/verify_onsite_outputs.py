"""Validate all onsite label bijections and compare serial, MPI and debug records."""
import hashlib
import json
import re
import sys
from pathlib import Path


def records(filename, rank=None):
    cases, labels, characters, largest, finished, screw, magnetic = {}, {}, {}, 0.0, False, False, False
    for line in Path(filename).read_text().splitlines():
        if rank is not None:
            match = re.fullmatch(r"\[\d+,(\d+)\]<stdout>: ?(.*)", line)
            if not match or int(match[1]) != rank:
                continue
            line = match[2]
        finished |= "Site-induced atomic band tests passed." in line
        screw |= "P6_3 two-site C3 labels, scalar/spin/grey screw-phase induction: passed." in line
        magnetic |= "Magnetic C2*T stabilizer and site-exchanging half translation: passed." in line
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "ONSITE_MATCH":
            key = tuple(map(int, fields[1:4]))
            assert key not in cases, (filename, rank, key)
            cases[key] = (int(fields[4]), int(fields[5]))
            largest = max(largest, float(fields[6]))
        elif fields[0] == "ONSITE_REFERENCE":
            key = tuple(map(int, fields[1:4]))
            labels.setdefault(key, []).append((int(fields[4]), fields[5]))
        elif fields[0] == "ONSITE_CHARACTER":
            key = (*map(int, fields[1:5]), fields[5])
            assert key not in characters
            characters[key] = complex(float(fields[6]), float(fields[7]))
    assert finished and screw and magnetic
    assert set(cases) == {(sg, spin, variant) for sg in range(1, 231)
                         for spin in (0, 1) for variant in (1, 2, 3)}
    assert set(labels) == set(cases)
    assert sum(size for size, order in cases.values()) == 6384 and largest < 2e-4
    rows, expected_characters = [], set()
    for key, (size, order) in sorted(cases.items()):
        columns = labels[key]
        assert len(columns) == size
        assert {column for column, _ in columns} == set(range(1, size + 1))
        assert len({label for _, label in columns}) == size
        expected_characters.update((*key, operation, label) for operation in range(1, order + 1)
                                   for _, label in columns)
        rows.extend(f"{key} {column} {label}" for column, label in sorted(columns))
    assert expected_characters == characters.keys()
    coordinate_error = 0.0
    for sg in range(1, 231):
        for spin in (0, 1):
            names = [{label for _, label in labels[sg, spin, variant]} for variant in (1, 2, 3)]
            assert names[0] == names[1] == names[2]
            for variant in (2, 3):
                for operation in range(1, cases[sg, spin, 1][1] + 1):
                    for label in names[0]:
                        difference = abs(characters[sg, spin, variant, operation, label]
                                         - characters[sg, spin, 1, operation, label])
                        coordinate_error = max(coordinate_error, difference)
    assert coordinate_error < 1e-9, (filename, rank, coordinate_error)
    return rows, characters, largest, coordinate_error


serial, characters, largest, coordinate_error = records(sys.argv[1])
reorderings = {}
character_error = 0.0
for filename, rank in ((sys.argv[2], 0), (sys.argv[2], 1), (sys.argv[3], None)):
    other, actual, error, variant_error = records(filename, rank)
    assert actual.keys() == characters.keys()
    deviation = max(abs(actual[k] - characters[k]) for k in characters)
    assert deviation < 1e-9, (filename, rank, deviation)
    character_error = max(character_error, deviation)
    reorderings[f"{filename}:{rank}"] = sum(a != b for a, b in zip(other, serial))
    largest = max(largest, error)
    coordinate_error = max(coordinate_error, variant_error)
print(json.dumps(dict(reference_gamma_tables=460, coordinate_variants=3,
                      cases=1380, label_records=len(serial), largest_residual=largest,
                      character_records=len(characters), largest_cross_build_character_error=character_error,
                      largest_cross_coordinate_character_error=coordinate_error,
                      native_label_order_differences=reorderings,
                      builds=["ssmp", "psmp-rank-0", "psmp-rank-1", "instrumented"],
                      key_sha256=hashlib.sha256(("\n".join(map(str, sorted(characters)))+"\n").encode()).hexdigest()), indent=2))
