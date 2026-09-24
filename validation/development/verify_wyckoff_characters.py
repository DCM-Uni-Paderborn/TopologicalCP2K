"""Audit complete actual-site character records against fixture coverage and each build."""
import hashlib
import io
import json
import re
import sys
import tarfile
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def text_stream(specification):
    """Accept a local file or archive.tar.gz::relative/member without extraction."""
    if "::" not in specification:
        with Path(specification).open() as stream:
            yield stream
    else:
        filename, member = specification.split("::", 1)
        with tarfile.open(filename, "r:gz") as archive:
            info = archive.getmember(member)
            assert info.isfile()
            with io.TextIOWrapper(archive.extractfile(info)) as stream:
                yield stream


def lines(filename):
    with text_stream(filename) as stream:
        for line in stream:
            yield line.rstrip("\r\n")


def key_digest(records):
    digest = hashlib.sha256()
    for key in sorted(records):
        digest.update((str(key)+"\n").encode())
    return digest.hexdigest()


def read_records(filename, expected, rank=None):
    cases, chars, bloch, induction = {}, {}, {}, {}
    finished = False
    for line in lines(filename):
        if rank is not None:
            match = re.fullmatch(r"\[\d+,(\d+)\]<stdout>: ?(.*)", line)
            if match is None or int(match[1]) != rank:
                continue
            line = match[2]
        finished |= "Site-induced atomic band tests passed." in line
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "WYCKOFF_ONSITE":
            key = (int(fields[1]), fields[2], *map(int, fields[3:6]))
            assert key not in cases
            cases[key] = (*map(int, fields[6:10]), float(fields[10]))
        elif fields[0] == "WYCKOFF_CHARACTER":
            key = (int(fields[1]), fields[2], *map(int, fields[3:8]), fields[8])
            assert key not in chars
            chars[key] = complex(float(fields[9]), float(fields[10]))
        elif fields[0] == "WYCKOFF_BLOCH":
            key = (int(fields[1]), fields[2], *map(int, fields[3:9]), fields[9])
            assert key not in bloch
            bloch[key] = complex(float(fields[10]), float(fields[11]))
        elif fields[0] == "WYCKOFF_INDUCTION":
            key = (int(fields[1]), fields[2], *map(int, fields[3:8]))
            assert key not in induction
            induction[key] = (int(fields[8]), int(fields[9]), float(fields[10]), float(fields[11]))
    assert finished and cases.keys() == expected.keys(), (filename, rank, len(cases), len(expected))
    largest = 0.0
    for key, value in cases.items():
        assert value[:4] == expected[key], (key, value, expected[key])
        largest = max(largest, value[4])
    assert largest < 2e-4
    labels = {}
    for key in chars:
        case, pg, operation, label = key[:5], key[5], key[6], key[7]
        assert case in cases and pg == cases[case][0] and 1 <= operation <= cases[case][2]
        labels.setdefault(case, set()).add(label)
    assert labels.keys() == cases.keys()
    for key, names in labels.items():
        pg, orbit, order, ni, error = cases[key]
        assert len(names) == ni
        assert names == set(by_site[key[:4]]["labels"]), (key, names, by_site[key[:4]]["labels"])
        assert all((*key, pg, operation, name) in chars
                   for name in names for operation in range(1, order+1))
    coordinate_error = 0.0
    for key, value in chars.items():
        if key[4] != 2:
            continue
        other = (*key[:4], 1, *key[5:])
        coordinate_error = max(coordinate_error, abs(value-chars[other]))
    assert coordinate_error < 1e-9, (filename, rank, coordinate_error)
    assert induction.keys() == expected_induction.keys()
    for key, (size, ni, scaled, absolute) in induction.items():
        assert (size, ni) == expected_induction[key]
        assert scaled < 2e-4
    assert bloch.keys() == expected_bloch
    bloch_coordinate_error = max(abs(value-bloch[(*key[:4], 1, *key[5:])])
                                 for key, value in bloch.items() if key[4] == 2)
    assert bloch_coordinate_error < 1e-9, (filename, rank, bloch_coordinate_error)
    return (cases, chars, largest, coordinate_error, bloch, bloch_coordinate_error,
            max(row[2] for row in induction.values()), max(row[3] for row in induction.values()))


with text_stream(sys.argv[1]) as source:
    manifest = json.load(source)
rows = manifest["records"]
assert len(rows) == 6924
sites = {(row["sg"], row["site"]) for row in rows}
assert len(sites) == 1731 and {sg for sg, site in sites} == set(range(1, 231))
expected = {(row["sg"], row["site"], row["spin"], row["grey"], variant):
            (row["reference_group"], row["orbit"], row["unitary_stabilizer"], row["columns"])
            for row in rows for variant in (1, 2)}
assert len(expected) == 13848
by_site = {(row["sg"], row["site"], row["spin"], row["grey"]): row for row in rows}
expected_induction, expected_bloch = {}, set()
for row in rows:
    for variant in (1, 2):
        for query, members in enumerate(row["reciprocal_members"], 1):
            key = (row["sg"], row["site"], row["spin"], row["grey"], variant, row["reference_group"], query)
            expected_induction[key] = (len(members), row["columns"])
            expected_bloch.update((*key, operation, label) for operation in members for label in row["labels"])
(cases, reference, largest, coordinate_error, bloch_reference, bloch_coordinate_error,
 induction_error, absolute_error) = read_records(sys.argv[2], expected)
cross_error, bloch_cross_error = 0.0, 0.0
for filename, rank in ((sys.argv[3], 0), (sys.argv[3], 1), (sys.argv[4], None)):
    _, values, error, change, bloch, bloch_change, scaled, absolute = read_records(filename, expected, rank)
    assert values.keys() == reference.keys()
    deviation = max(abs(values[key]-reference[key]) for key in reference)
    assert deviation < 1e-9, (filename, rank, deviation)
    cross_error = max(cross_error, deviation)
    coordinate_error = max(coordinate_error, change)
    largest = max(largest, error)
    deviation = max(abs(value-bloch_reference[key]) for key, value in bloch.items())
    assert deviation < 1e-9, (filename, rank, deviation)
    bloch_cross_error = max(bloch_cross_error, deviation)
    bloch_coordinate_error = max(bloch_coordinate_error, bloch_change)
    induction_error = max(induction_error, scaled)
    absolute_error = max(absolute_error, absolute)
print(json.dumps(dict(actual_wyckoff_families=len(sites), variants=4, coordinate_variants=2,
                      complete_table_matches=len(cases), onsite_columns=sum(row["columns"] for row in rows),
                      character_entries=len(reference), largest_reference_residual=largest,
                      largest_cross_build_character_difference=cross_error,
                      largest_cross_coordinate_character_difference=coordinate_error,
                      induced_tables=len(expected_induction), induced_character_entries=len(bloch_reference),
                      largest_induced_reference_dimension_scaled_error=induction_error,
                      largest_induced_reference_absolute_error=absolute_error,
                      largest_induced_cross_build_error=bloch_cross_error,
                      largest_induced_cross_coordinate_error=bloch_coordinate_error,
                      fixture_sha256=manifest["fixture_sha256"],
                      character_key_sha256=key_digest(reference),
                      induced_character_key_sha256=key_digest(bloch_reference)), indent=2))
