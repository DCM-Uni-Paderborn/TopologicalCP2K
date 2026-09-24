"""Require identical site-induction result records in optimized, instrumented and two-rank runs."""
import json
import re
import sys
from pathlib import Path

prefixes = ("Site induction", "Site-induction", "Site-specialization", "Maximal site", "Transitive site")


def select(text):
    return [line.strip() for line in text.splitlines() if line.strip().startswith(prefixes)]


serial, debug, mpi = [Path(path).read_text() for path in sys.argv[1:]]
reference = select(serial)
assert reference[-1] == "Site-specialization induction tests passed."
assert len(reference) == 4774 and select(debug) == reference
for rank in (0, 1):
    prefix = f"[1,{rank}]<stdout>:"
    records = select("\n".join(re.sub(r"^" + re.escape(prefix) + " ?", "", line)
                               for line in mpi.splitlines() if line.startswith(prefix)))
    assert records == reference, rank
print(json.dumps(dict(records=len(reference), debug_identical=True, mpi_identical=[0, 1])))
