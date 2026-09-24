"""Compare maximal-site generator counts with the Bilbao EBR tables archived by Crystalline.jl.

This is an independent count check, not a row-by-row canonical-irrep certificate.
Arguments: exception-comparison JSONL, download cache directory.
"""
import csv
import hashlib
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

revision = "9cfb644a1a4cc1c7baec457b321885459bb98bd1"
base = f"https://raw.githubusercontent.com/thchr/Crystalline.jl/{revision}/data/bandreps/3d"
cache = Path(sys.argv[2])
cache.mkdir(parents=True, exist_ok=True)
native = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
native = [row for row in native if "hall" in row]
assert len(native) == 2120
results = []
totals = [0, 0, 0, 0]
for grey, directory in [(1, "elementary"), (2, "elementaryTR")]:
    for sg in range(1, 231):
        filename = cache / f"{directory}-{sg}.csv"
        url = f"{base}/{directory}/maxpaths/{sg}.csv"
        if not filename.exists():
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url, timeout=30) as response:
                        data = response.read()
                    break
                except OSError:
                    if attempt == 2:
                        raise
                    time.sleep(2)
            filename.write_bytes(data)
        data = filename.read_bytes()
        table = list(csv.reader(io.StringIO(data.decode()), delimiter="|"))
        assert len(table) > 3 and table[1][0] == "Band-Rep."
        n = len(table[0])
        assert all(len(row) == n for row in table)
        spinful = [any("\u02e2" in row[col] or "\\bar" in row[col] for row in table[3:])
                   for col in range(1, n)]
        for spin in (0, 1):
            count = sum(flag == bool(spin) for flag in spinful)
            totals[2*(grey - 1) + spin] += count
            for row in native:
                if (row["sg"], row["grey"], row["spin"]) != (sg, grey, spin):
                    continue
                result = dict(hall=row["hall"], sg=sg, grey=grey, spin=spin,
                              native=row["remaining"], reference=count, match=row["remaining"] == count,
                              revision=revision, sha256=hashlib.sha256(data).hexdigest(), url=url)
                results.append(result)
                print(json.dumps(result), flush=True)
assert len(results) == 2120
print(json.dumps(dict(total=len(results), matched=sum(row["match"] for row in results),
                      mismatched=sum(not row["match"] for row in results), reference_totals=totals)))
if not all(row["match"] for row in results):
    raise SystemExit(1)
