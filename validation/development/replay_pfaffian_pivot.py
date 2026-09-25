"""Verify pivot-fix evidence and optionally rerun selected material queries."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from replay_stanene_convergence import digest, unpack


def exact_pfaffian(a):
    if not a:
        return 1
    return sum((-1)**(j + 1) * a[0][j] * exact_pfaffian([
        [a[k][l] for l in range(len(a)) if l not in (0, j)]
        for k in range(len(a)) if k not in (0, j)]) for j in range(1, len(a)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--recompute", nargs="*", default=[])
    args = parser.parse_args()
    base = args.bundle.resolve()
    index = json.loads((base / "index.json").read_text())
    original = base.parent / "stanene-convergence"
    assert digest(original / "index.json") == index["original_index_sha256"]
    prior_index = json.loads((original / "index.json").read_text())
    labels = {case["label"] for case in prior_index["cases"]}
    if args.recompute and (not args.library or set(args.recompute) - labels - {"all"}):
        parser.error("Recomputation requires --library and known case labels (or all)")
    with tempfile.TemporaryDirectory(prefix="pfaffian-pivot-replay-") as temporary:
        work = Path(temporary)
        current = work / "current"
        manifest = unpack(base / index["archive"]["file"], current, index["archive"]["sha256"])
        original_methods = next(a for a in prior_index["archives"] if a["file"] == "methods-results.tar.gz")
        unpack(original / original_methods["file"], work / "original", original_methods["sha256"])
        fixture = json.loads((current / "evidence/pivot-fixture-v2.json").read_text())
        assert exact_pfaffian([[int(x) for x in row] for row in fixture["matrix"]]) == 60
        before = json.loads((current / "evidence/fixture-v2-baseline.json").read_text())
        after = json.loads((current / "evidence/fixture-v2-typed.json").read_text())
        assert len(before["results"]) == len(after["results"]) == 6
        assert all(row["sign"] is None for row in before["results"])
        for row in after["results"]:
            order = row["ordering"]
            assert sorted(order) == list(range(8)) and row["scale"] in (1.0, 1e-120, 1e120)
            permuted = [[int(fixture["matrix"][i][j]) for j in order] for i in order]
            exact = exact_pfaffian(permuted)
            assert row["sign"] == (1 if exact > 0 else -1) and row["residual"] < 1e-10
        for name in ("baseline-local-blas-tzvp", "patched-homebrew-tzvp", "max-pivot-local-tzvp"):
            control = json.loads((current / "evidence" / f"{name}.json").read_text())
            assert len(control["queries"]) == 1
            row = control["queries"][0]
            assert row["eta"] == 0.75 and row["nu"] == 1 and row["factor_residual"] < 1e-10
        summary, recovered, unresolved, newly_unresolved, gaps = [], [], [], [], []
        for case in prior_index["cases"]:
            label = case["label"]
            old = json.loads((work / "original/reports" / f"{label}.json").read_text())
            new = json.loads((current / "reports" / f"{label}.json").read_text())
            assert new["inputs"]["sparse-final-scan.json"] == digest(work / "original/reports" / f"{label}.json")
            for name in ("stanene.topology", "stanene.mmn"):
                assert old["inputs"][name] == new["inputs"][name]
            assert new["reference_script_sha256"] == digest(current / "scripts/check_sparse_bloch_localizer.py")
            for left, right in zip(old["queries"], new["queries"], strict=True):
                assert left["eta"] == right["eta"] and left["energy"] == right["energy"]
                gaps.append(abs(left["gap"] - right["gap"]))
                assert gaps[-1] < 1e-10 and right["eigen_residual"] < 1e-8
                if left["nu"] is not None and right["nu"] is not None:
                    assert left["nu"] == right["nu"]
                if right["nu"] is None:
                    assert all(not a["resolved"] for a in right["factor_attempts"])
                    unresolved.append((label, right["eta"]))
                    if left["nu"] is not None:
                        newly_unresolved.append((label, right["eta"], left["nu"]))
                else:
                    assert right["nu"] in (0, 1) and right["factor_residual"] <= 1e-10
                    assert right["factor_attempts"][-1]["resolved"]
                    if left["nu"] is None:
                        recovered.append((label, right["eta"], right["nu"]))
                summary.append((label, right["eta"], right["nu"]))
            if "all" in args.recompute or label in args.recompute:
                if case["new_export"]:
                    archive = next(a for a in prior_index["archives"] if a["file"] == f"{label}.tar.gz")
                    case_path = work / label
                    unpack(original / archive["file"], case_path, archive["sha256"])
                else:
                    previous = work / "previous"
                    if not previous.exists():
                        unpack(original / prior_index["previous_reference_archive"], previous,
                               prior_index["previous_reference_sha256"])
                    case_path = previous / case["directory"]
                command = [sys.executable, str(current / "scripts/check_sparse_bloch_localizer.py"),
                           str(case_path), "--library", str(args.library.resolve()),
                           "--energy", str(case["energy"]), "--eta",
                           *[str(q["eta"]) for q in new["queries"]], "--output", "recomputed-pivot.json"]
                if label == "mesh3":
                    command += ["--dense-check"]
                subprocess.run(command, check=True)
                recomputed = json.loads((case_path / "recomputed-pivot.json").read_text())
                for left, right in zip(new["queries"], recomputed["queries"], strict=True):
                    assert abs(left["gap"] - right["gap"]) < 1e-10 and left["nu"] == right["nu"]
        assert len(summary) == 33
        result = dict(verified_files=len(manifest["files"]), queries=len(summary),
                      recovered=recovered, unresolved=unresolved, newly_unresolved=newly_unresolved,
                      maximum_gap_change=max(gaps),
                      recomputed=args.recompute)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
