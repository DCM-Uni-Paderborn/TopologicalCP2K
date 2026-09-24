"""Compare deterministic affine-reference records across serial/MPI/debug logs.

Arguments: serial log, tagged two-rank MPI log, instrumented no-database log.
Runtime diagnostics are retained in raw logs, but are not scientific records.
"""
import hashlib
import json
import re
import sys
from pathlib import Path

prefixes = ("AFFINE_MATCH ", "SITE_COLUMN ", "Affine convention reference passed:")


def records(text):
    return [" ".join(line.split()) for line in text.splitlines()
            if line.strip().startswith(prefixes)]


serial = Path(sys.argv[1]).read_text()
mpi = Path(sys.argv[2]).read_text()
debug = Path(sys.argv[3]).read_text()
output = {"serial": records(serial), "debug": records(debug)}
assert "Independent affine-reference groups passed:" in serial
assert "Independent affine-reference groups passed:" in debug
for rank in (0, 1):
    lines = []
    for line in mpi.splitlines():
        match = re.match(rf"\[1,{rank}\]<stdout>:(.*)", line)
        if match:
            lines.append(match[1])
    assert any("Independent affine-reference groups passed:" in line for line in lines)
    output[f"mpi{rank}"] = records("\n".join(lines))
assert all(rows == output["serial"] for rows in output.values())
assert sum(row.startswith("AFFINE_MATCH ") for row in output["serial"]) == 6924
assert sum(row.startswith("Affine convention") for row in output["serial"]) == 230
result = {name: dict(records=len(rows), sha256=hashlib.sha256("\n".join(rows).encode()).hexdigest())
          for name, rows in output.items()}
print(json.dumps(dict(identical=True, outputs=result), indent=2))
