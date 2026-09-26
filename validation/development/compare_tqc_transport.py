"""Compare ordered native TQC records without relaxing integer certificates.

Usage: python compare_tqc_transport.py REFERENCE_DIR CANDIDATE_DIR
Optional --extra permits one explicitly named new candidate file.
Floating fields use an absolute 1e-10 tolerance; all other tokens are exact.
This checks transport consistency, not scientific completeness of the graph.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re


REAL = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[EeDd][+-]?\d+)?\Z")


def compare(left, right):
    a, b = left.splitlines(), right.splitlines()
    if len(a) != len(b):
        raise ValueError("Record counts differ")
    maximum = 0.0
    counts = Counter()
    numeric_fields = 0
    for lineno, (x, y) in enumerate(zip(a, b), 1):
        xt, yt = x.split(), y.split()
        if len(xt) != len(yt):
            raise ValueError(f"Changed record shape at line {lineno}")
        if xt:
            counts[xt[0]] += 1
        for u, v in zip(xt, yt):
            if u.lower() in {"nan", "inf", "infinity"} or v.lower() in {"nan", "inf", "infinity"}:
                raise ValueError(f"Nonfinite token at line {lineno}")
            if REAL.fullmatch(u) and REAL.fullmatch(v):
                error = abs(float(u.replace("D", "E")) - float(v.replace("D", "E")))
                if not math.isfinite(error) or error > 1e-10:
                    raise ValueError(f"Floating mismatch at line {lineno}: {u}, {v}")
                maximum = max(maximum, error)
                numeric_fields += 1
            elif u != v:
                raise ValueError(f"Exact-token mismatch at line {lineno}: {u}, {v}")
    return dict(records=len(a), types=dict(counts), floating_fields=numeric_fields,
                maximum_absolute_difference=maximum)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--extra", action="append", default=[])
    args = parser.parse_args()
    a = {p.name: p for p in args.reference.glob("*.little_group")}
    b = {p.name: p for p in args.candidate.glob("*.little_group")}
    if not a or set(b) - set(a) != set(args.extra) or set(a) - set(b):
        raise ValueError("Missing or unexpected native result files")
    results = []
    for name in sorted(a):
        item = compare(a[name].read_text(), b[name].read_text())
        item.update(name=name, reference_sha256=hashlib.sha256(a[name].read_bytes()).hexdigest(),
                    candidate_sha256=hashlib.sha256(b[name].read_bytes()).hexdigest())
        results.append(item)
    print(json.dumps(dict(files=len(results), results=results, accepted=True), indent=2))


if __name__ == "__main__":
    main()
