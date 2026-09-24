"""Audit exact EBR reference coverage and preserve ambiguous character matches."""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path


def read_log(path, manifest):
    expected = {(c["sg"], c["site"], c["spin"]): c for c in manifest["cases"]}
    rows, candidates, completed = defaultdict(dict), defaultdict(lambda: defaultdict(list)), set()
    for line in path.read_text().splitlines():
        match = re.match(r"\[\d+,(\d+)\]<stdout>:(.*)", line)
        rank, text = (int(match[1]), match[2]) if match else (0, line)
        fields = text.split()
        if not fields:
            continue
        if fields[0] == "EBR_REFERENCE":
            sg, site, spin, run, native, count, points, unique, ambiguous, error = fields[1:]
            key = (int(sg), site, int(spin), int(run))
            assert key not in rows[rank]
            row = dict(native=int(native), count=int(count), points=int(points), unique=int(unique),
                       ambiguous=int(ambiguous), error=float(error))
            reference = expected[key[:3]]
            assert row["count"] == reference["columns"] and row["points"] == reference["points"]
            assert row["unique"]+row["ambiguous"] == row["count"]
            assert 0 <= row["error"] < 2e-4
            rows[rank][key] = row
        elif fields[0] == "EBR_CANDIDATE":
            sg, site, spin, run, native, label, error = fields[1:]
            key = (int(sg), site, int(spin), int(run))
            assert 0 <= float(error) < 2e-4
            candidates[rank][key].append((int(native), label, float(error)))
        elif text.strip() == "Site-induced atomic band tests passed.":
            completed.add(rank)
    assert set(rows) == completed
    target = {(*k, run) for k in expected for run in (1, 2)}
    summaries = []
    baseline = None
    for rank, records in sorted(rows.items()):
        assert set(records) == target and set(candidates[rank]) == target
        signatures = {}
        for key, row in records.items():
            by_label = defaultdict(list)
            for native, label, error in candidates[rank][key]:
                assert 1 <= native <= row["native"]
                by_label[label].append(native)
            assert set(by_label) == {name.replace("\u2191", "_up_") for name in expected[key[:3]]["labels"]}
            assert all(len(v) == len(set(v)) for v in by_label.values())
            assert sum(len(v) == 1 for v in by_label.values()) == row["unique"]
            signatures[key] = {label: len(v) for label, v in by_label.items()}
        for key in expected:
            assert signatures[(*key, 1)] == signatures[(*key, 2)]
        if baseline is None:
            baseline = signatures
        else:
            assert signatures == baseline
        normal = [v for k, v in records.items() if k[-1] == 1]
        summaries.append(dict(rank=rank, groups=len({k[0] for k in expected}), site_spin_cases=len(expected),
            columns=sum(r["count"] for r in normal), unique=sum(r["unique"] for r in normal),
            ambiguous=sum(r["ambiguous"] for r in normal), conventions=2,
            native_columns=sum(r["native"] for r in normal),
            largest_residual=max(r["error"] for r in records.values())))
    assert summaries
    return summaries, baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    results, reference = [], None
    for path in args.logs:
        summaries, signatures = read_log(path, manifest)
        if reference is None:
            reference = signatures
        else:
            assert signatures == reference
        results.append(dict(log=str(path), ranks=summaries))
    print(json.dumps(dict(results=results,
        queries_per_convention=sum(r["points"] for r in manifest["cases"]),
        character_entries_per_convention=sum(r["entries"] for r in manifest["cases"]),
        supplemental_tables=len(manifest["supplements"]),
        supplemental_points=sum(len(r["added_points"]) for r in manifest["supplements"]),
        largest_shared_table_residual=max(r["shared_character_residual"] for r in manifest["supplements"])), indent=2))


if __name__ == "__main__":
    main()
