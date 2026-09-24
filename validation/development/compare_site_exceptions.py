"""Compare real-space composite detection with Cano et al. PRB 97, 035139 (2018).

The reference sets are the space-group coverage of Tables II-IV, not a
row-by-row canonical-irrep identification or a complete EBR database.
Arguments: native site-unit --all-settings log, libsymspg path.
"""

import ctypes as ct
import json
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


scalar = {163, 165, 167, 228, 230, 223, 211, 208, 210, 188, 190, 192, 193,
          207, 222, 124, 140, 229, 226, 215, 217, 224, 131, 132, 139}
scalar_grey = {229, 226, 215, 217, 224, 131, 132, 139, 140, 223,
               84, 87, 135, 136, 112, 116, 120, 121, 126, 130, 133, 138, 142,
               218, 230, 222, 219, 228}
spinor = {224, 227, 225, 223, 211, 208, 210, 228, 188, 190, 192, 193,
          163, 165, 167, 230, 194, 207, 222, 124, 140, 183, 51, 63, 67, 74,
          138, 99, 107, 115, 137, 195, 197, 201, 209, 218, 177, 214, 112,
          116, 120, 121, 126, 130, 133, 142, 111, 132, 134, 49, 66, 69, 72,
          128, 135, 89, 97}
reference = {(1, 0): scalar, (1, 1): spinor, (2, 0): scalar_grey, (2, 1): set()}
lib = ct.CDLL(sys.argv[2])
lib.spg_get_spacegroup_type.argtypes = [ct.c_int]
lib.spg_get_spacegroup_type.restype = SpacegroupType
text = Path(sys.argv[1]).read_text()
assert "Site-specialization induction tests passed." in text
rows = []
seen = set()
sg_rows = {}
for line in text.splitlines():
    if "Maximal site Hall/grey/spin/composite/remaining:" not in line:
        continue
    hall, grey, spin, composite, remaining = map(int, line.split(":")[1].split())
    key = (hall, grey, spin)
    assert key not in seen
    seen.add(key)
    group = lib.spg_get_spacegroup_type(hall)
    assert group.hall_number == hall and 1 <= group.number <= 230
    expected = group.number in reference[grey, spin]
    row = dict(hall=hall, sg=group.number, symbol=group.international_short.decode(),
               grey=grey, spin=spin, composite=composite, remaining=remaining,
               expected_exception=expected, match=(composite > 0) == expected)
    rows.append(row)
    sg_rows.setdefault((group.number, grey, spin), row)
assert len(rows) == 2120 and len(sg_rows) == 920
for row in rows:
    print(json.dumps(row))
totals = {f"grey={g},spin={s}": sum(x["remaining"] for x in sg_rows.values()
                                    if x["grey"] == g and x["spin"] == s)
          for g, s in reference}
print(json.dumps(dict(total=len(rows), matched=sum(x["match"] for x in rows),
                      mismatched=sum(not x["match"] for x in rows), remaining_maximal_columns=totals)))
if not all(x["match"] for x in rows):
    raise SystemExit(1)
