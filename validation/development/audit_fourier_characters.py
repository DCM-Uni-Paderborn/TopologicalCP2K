"""Audit per-rank invariant-torus comparisons without interpreting bundle equivalence."""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path


def read_log(path):
    rows = defaultdict(dict)
    totals = {}
    completed = set()
    for line in path.read_text().splitlines():
        match = re.match(r"\[\d+,(\d+)\]<stdout>:(.*)", line)
        rank, text = (int(match[1]), match[2]) if match else (0, line)
        fields = text.split()
        if not fields:
            continue
        if fields[0] == "EBR_FOURIER":
            sg, site, spin, run, label, count, separated = fields[1:]
            key = (int(sg), site, int(spin), int(run), label)
            assert key not in rows[rank]
            assert int(count) > 1 and 0 <= int(separated) < int(count)
            rows[rank][key] = (int(count), int(separated))
        elif text.strip().startswith("Fourier character bands/queries/largest residual:"):
            assert rank not in totals
            values = text.split(":", 1)[1].split()
            totals[rank] = dict(bands=int(values[0]), queries=int(values[1]), residual=float(values[2]))
            assert 0 <= totals[rank]["residual"] < 1e-8
        elif text.strip() == "Site-induced atomic band tests passed.":
            completed.add(rank)
    assert completed == set(rows) == set(totals) and completed
    summaries, baseline = [], None
    for rank, records in sorted(rows.items()):
        normal = {k[:3] + k[4:]: v for k, v in records.items() if k[3] == 1}
        changed = {k[:3] + k[4:]: v for k, v in records.items() if k[3] == 2}
        assert len(normal) == 378 and normal == changed and len(records) == 756
        if baseline is not None:
            assert normal == baseline
        baseline = normal
        summaries.append(dict(rank=rank, ambiguous_references=378,
                              references_with_distinct_character_functions=sum(v[1] > 0 for v in normal.values()),
                              references_with_equal_character_functions=sum(v[1] == 0 for v in normal.values()),
                              **totals[rank]))
    return summaries, baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()
    results, baseline = [], None
    for path in args.logs:
        summary, records = read_log(path)
        if baseline is not None:
            assert records == baseline
        baseline = records
        results.append(dict(log=str(path), ranks=summary))
    print(json.dumps(dict(results=results, interpretation="Unitary characters only; no Bloch-bundle equivalence claim"), indent=2))


if __name__ == "__main__":
    main()
