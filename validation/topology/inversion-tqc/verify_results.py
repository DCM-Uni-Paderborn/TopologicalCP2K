"""Check frozen CP2K inversion tests with the version-matched test matchers."""

import argparse
import dataclasses
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("cp2k_source", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parent
sys.path.insert(0, str(args.cp2k_source / "tests"))
from matchers import run_matcher

specs = tomllib.loads((root / "regression" / "TEST_FILES.toml").read_text())
report = {"source_commit": "d97839c4bc47c7d921137641d3a06086886a6056", "regression": {}}
for name, checks in specs.items():
    path = root / "regression" / (name + ".out")
    output = path.read_text()
    assert output.count("PROGRAM ENDED AT") == 1 and "[ABORT]" not in output, name
    assert re.search(r"source code revision number:\s+d97839c\s", output), (name, "source revision")
    matches = []
    for check in checks:
        result = run_matcher(output, **check)
        assert result.status == "OK", (name, check, result)
        matches.append({"matcher": check["matcher"], **dataclasses.asdict(result)})
    report["regression"][name] = matches

pattern = re.compile(r"TQC\| TRIM\s+([\d. -]+)even/odd states\s+(\d+)\s+(\d+)\s+residuals/gap\s+(\S+)\s+(\S+)\s+(\S+)")
expected = {"neon-tqc.inp": [3, 1, 1, 3, 1, 3, 3, 1],
            "neon-tqc-plane.inp": [3, 3, 3, 3],
            "stanene-tqc.inp": [1, 2, 2, 2]}
report["parities"] = {}
for name, odd in expected.items():
    output = (root / "regression" / (name + ".out")).read_text()
    rows = pattern.findall(output)
    assert len(rows) == len(odd), name
    assert [int(row[2]) // 2 for row in rows] == odd, name
    assert all(int(row[1]) + int(row[2]) == 8 and int(row[2]) % 2 == 0 for row in rows)
    residual = max(float(row[3]) for row in rows)
    commutator = max(float(row[4]) for row in rows)
    minimum_gap = min(float(row[5]) for row in rows)
    assert residual < 1e-6 and commutator < 1e-6 and minimum_gap > 1e-6
    report["parities"][name] = {"odd_kramers_pairs": odd, "max_residual": residual,
                              "max_energy_commutator_hartree": commutator,
                              "minimum_trim_gap_hartree": minimum_gap}

serial = (root / "stanene-serial.out").read_text()
assert "PROGRAM ENDED AT" in serial and "[ABORT]" not in serial
assert re.search(r"source code revision number:\s+d97839c\s", serial)
assert "TQC| Fu-Kane parity index: 1" in serial
assert [int(row[2]) // 2 for row in pattern.findall(serial)] == expected["stanene-tqc.inp"]
assert "Native characters, parity, inversion indicators and EBR tests passed." in (root / "unit.txt").read_text()
negative = (root / "bad-center.out").read_text()
assert "[ABORT]" in negative and "inversion does not preserve geometry and atomic kinds" in negative
assert "TQC| Fu-Kane parity index" not in negative
report["unit_candidate_signatures"] = {"2D": 81, "3D": 6561}
report["serial_stanene"] = "PASS"
report["wrong_inversion_center"] = "REJECTED"
report["source_patch_sha256"] = hashlib.sha256((root / "source.patch").read_bytes()).hexdigest()
print(json.dumps(report, indent=2))
