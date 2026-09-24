"""Compare native graph quotients with Tables I-IV, IX-XII of arXiv:1703.00911v2.

Arguments: pdftotext -layout output, native --all-settings log, libsymspg path.
This is a development oracle, not a CP2K dependency or a completeness proof.
"""
import ctypes as ct
import json
import re
import sys
from pathlib import Path


class SpacegroupType(ct.Structure):
    _fields_ = [("number", ct.c_int), ("international_short", ct.c_char * 11),
                ("international_full", ct.c_char * 20), ("international", ct.c_char * 32),
                ("schoenflies", ct.c_char * 7), ("hall_number", ct.c_int),
                ("hall_symbol", ct.c_char * 17), ("choice", ct.c_char * 6),
                ("pointgroup_international", ct.c_char * 6),
                ("pointgroup_schoenflies", ct.c_char * 4),
                ("arithmetic_crystal_class_number", ct.c_int),
                ("arithmetic_crystal_class_symbol", ct.c_char * 7)]


def read_table(text, roman, quotient=False):
    section = text.split(f"TABLE {roman}.", 1)[1]
    section = section.split("Space groups", 1)[1]
    section = section.split("XBS :" if quotient else "d:", 1)[0]
    result = {sg: [] for sg in range(1, 231)} if quotient else {}
    seen = set()
    value = None
    for line in section.splitlines():
        if not line.strip():
            continue
        if quotient:
            match = re.match(r"\s*(Z\d+(?:\s*\N{MULTIPLICATION SIGN}\s*Z\d+)*)\s+(\d.*)$", line)
        else:
            match = re.match(r"\s*(\d+)\s+(\d[\d,\s]*)$", line)
        if match:
            value = list(map(int, re.findall(r"\d+", match[1]))) if quotient else int(match[1])
            numbers = match[2]
        else:
            assert value is not None and re.fullmatch(r"[\d,\s]+", line), line
            numbers = line
        for sg in map(int, re.findall(r"\d+", numbers)):
            assert 1 <= sg <= 230 and sg not in seen, (roman, sg)
            seen.add(sg)
            result[sg] = value
    assert quotient or len(seen) == 230, (roman, len(seen))
    return result


text = Path(sys.argv[1]).read_text()
references = {}
for grey, spin, rank_table, quotient_table in [(1, 1, "I", "III"), (1, 0, "II", "IV"),
                                              (0, 1, "IX", "XI"), (0, 0, "X", "XII")]:
    references[grey, spin] = (read_table(text, rank_table), read_table(text, quotient_table, True))
lib = ct.CDLL(sys.argv[3])
lib.spg_get_spacegroup_type.argtypes = [ct.c_int]
lib.spg_get_spacegroup_type.restype = SpacegroupType
pending = []
results = []
for line in Path(sys.argv[2]).read_text().splitlines():
    if "Graph quotient spin/coordinates/rank/free/cyclic:" in line:
        spin, coordinates, rank, free, cyclic = map(int, line.split(":")[1].split())
        pending.append(dict(spin=spin, rank=rank, free=free, cyclic=cyclic))
    elif "Graph cyclic orders:" in line:
        pending[-1]["orders"] = list(map(int, line.split(":")[1].split()))
    elif "Reciprocal Hall, grey, strata:" in line:
        hall, grey, _ = map(int, line.split(":")[1].split())
        sg = lib.spg_get_spacegroup_type(hall)
        assert sg.hall_number == hall and 1 <= sg.number <= 230
        assert len(pending) >= 2
        for row in pending[-2:]:
            ranks, groups = references[grey, row["spin"]]
            row.update(hall=hall, sg=sg.number, symbol=sg.international_short.decode(), grey=grey,
                       reference_rank=ranks[sg.number], reference_orders=groups[sg.number])
            assert row["cyclic"] == len(row["orders"])
            row["match"] = row["rank"] == row["reference_rank"] and row["free"] == 0 and row["orders"] == row["reference_orders"]
            results.append(row)
        pending = []
assert len(results) == 2120 and len({(x["hall"], x["grey"], x["spin"]) for x in results}) == 2120
assert "Reciprocal-stratum tests passed." in Path(sys.argv[2]).read_text()
for row in results:
    print(json.dumps(row))
print(json.dumps(dict(total=len(results), matched=sum(x["match"] for x in results),
                      mismatched=sum(not x["match"] for x in results))))
raise SystemExit(0 if all(row["match"] for row in results) else 1)
