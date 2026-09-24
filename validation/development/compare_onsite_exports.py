"""Compare named onsite provenance between the serial and MPI export audits."""
import json
import sys
from pathlib import Path

import numpy as np


def read(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines()]
    result = {Path(row["file"]).name: row for row in rows}
    assert len(result) == len(rows) == 16
    return result


serial, mpi = map(read, sys.argv[1:])
assert serial.keys() == mpi.keys()
axis_error = 0.0
for name, left in serial.items():
    right = mpi[name]
    for key in ("sites", "columns", "points", "entries", "rejected_corruptions", "onsite_reference_sha256"):
        assert left[key] == right[key], (name, key)
    assert len(left["rejected_corruptions"]) == 4
    assert len(left["reference_matches"]) == left["sites"]
    assert len(right["reference_matches"]) == right["sites"]
    for a, b in zip(left["reference_matches"], right["reference_matches"]):
        for key in ("site", "point_group_reference", "spinful", "grey", "columns", "convention"):
            assert a[key] == b[key], (name, key)
        axis_error = max(axis_error, float(np.max(np.abs(np.array(a["reference_to_export_axes"])-b["reference_to_export_axes"]))))
assert axis_error < 1e-9
records = list(serial.values()) + list(mpi.values())
print(json.dumps(dict(cases_per_build=16, builds=2, sites_per_build=sum(r["sites"] for r in serial.values()),
                      columns_per_build=sum(r["columns"] for r in serial.values()),
                      reciprocal_queries_per_build=sum(r["points"] for r in serial.values()),
                      induced_entries_per_build=sum(r["entries"] for r in serial.values()),
                      largest_induction_error=max(r["maximum_error"] for r in records),
                      largest_reference_error=max(s["maximum_reference_residual"] for r in records for s in r["reference_matches"]),
                      largest_cross_build_axis_error=axis_error,
                      corruptions_rejected=sum(len(r["rejected_corruptions"]) for r in records)), indent=2))
