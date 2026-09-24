"""Compare complete scalar/spin reference correspondences across local builds."""
import hashlib
import json
import re
import sys
from pathlib import Path


def records(filename, rank=None):
    rows, cases, names, largest, completed = [], {}, {}, 0.0, False
    for line in Path(filename).read_text().splitlines():
        if rank is not None:
            match = re.fullmatch(r"\[\d+,(\d+)\]<stdout>: ?(.*)", line)
            if not match or int(match[1]) != rank:
                continue
            line = match[2]
        completed |= 'Character correspondence tests passed.' in line
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'CHARACTER_MATCH':
            key = tuple(fields[1:4])
            assert key not in cases
            cases[key] = int(fields[4])
            largest = max(largest, float(fields[5]))
            rows.append(' '.join(fields[:-1]))
        elif fields[0] == 'IRREP_REFERENCE':
            rows.append(' '.join(fields))
            names.setdefault(tuple(fields[1:4]), []).append(fields)
    assert completed and len(cases) == 2700 and sum(cases.values()) == 8907
    assert {(int(key[0]), int(key[1])) for key in cases} == {
        (sg, spin) for sg in range(1, 231) for spin in (0, 1)}
    assert len(rows) == 11607 and largest < 2e-4
    for key, size in cases.items():
        named = names.get(key, [])
        assert len(named) == size
        assert {int(row[4]) for row in named} == set(range(1, size+1))
        assert len({row[5] for row in named}) == size
    return sorted(rows), largest


serial, largest = records(sys.argv[1])
for name, rank in ((sys.argv[2], 0), (sys.argv[2], 1), (sys.argv[3], None)):
    other, error = records(name, rank)
    assert other == serial, (name, rank)
    largest = max(largest, error)
print(json.dumps(dict(groups=230, variants=2, cases=2700, irreps=8907,
                      builds=['ssmp', 'psmp-rank-0', 'psmp-rank-1', 'instrumented'],
                      records=len(serial), largest_printed_residual=largest,
                      sha256=hashlib.sha256(('\n'.join(serial)+'\n').encode()).hexdigest()), indent=2))
