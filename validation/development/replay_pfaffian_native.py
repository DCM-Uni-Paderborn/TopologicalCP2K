"""Check native evidence against the independently verified prototype archive."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from replay_stanene_convergence import digest, unpack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    args = parser.parse_args()
    base = args.bundle.resolve()
    index = json.loads((base / "index.json").read_text())
    previous = base.parent / "pfaffian-delayed"
    assert digest(previous / "index.json") == index["previous_index_sha256"]
    previous_index = json.loads((previous / "index.json").read_text())
    with tempfile.TemporaryDirectory(prefix="pfaffian-native-replay-") as temporary:
        work = Path(temporary)
        manifest = unpack(base / index["archive"]["file"], work / "native", index["archive"]["sha256"])
        unpack(previous / previous_index["archive"]["file"], work / "previous", previous_index["archive"]["sha256"])
        subprocess.run([sys.executable, str(Path(__file__).with_name("replay_pfaffian_delayed.py")),
                        str(previous), str(work / "previous-replay.json")], check=True, stdout=subprocess.PIPE)
        baseline = json.loads((work / "previous-replay.json").read_text())
        assert baseline["material_queries"] == baseline["resolved"] == 33
        native, old = work / "native/build-mpi", work / "previous/build-serial/pfaffian-stability"
        controls = []
        for path in sorted((old / "delayed-material-recheck").glob("mesh*.json")):
            before = json.loads(path.read_text())
            after = json.loads((native / "delayed-native-materials-verified" / path.name).read_text())
            assert before["inputs"] == after["inputs"]
            assert after["library_sha256"] == index["material_library_sha256"]
            assert before["script_sha256"] == after["script_sha256"]
            for a, b in zip(before["queries"], after["queries"], strict=True):
                for key in ("eta", "energy", "flatten_scale", "nu"):
                    assert a[key] == b[key], (path.name, key)
                assert abs(a["gap"] - b["gap"]) < 1e-10
                assert b["eigen_residual"] < 1e-8 and b["factor_residual"] <= 1e-10
                assert len(b["factor_attempts"]) == 1
                controls.append(dict(case=path.stem, eta=b["eta"], nu=b["nu"], gap=b["gap"], residual=b["factor_residual"]))
        assert len(controls) == 33
        synthetic = []
        for name, count, resolved, singular, script in (
            ("independent", 600, 468, 132, "validate_adaptive.py"),
            ("stress", 636, 336, 300, "validate_delayed.py")):
            report = json.loads((native / f"delayed-native-verified-{name}.json").read_text())
            expected = dict(queries=count, resolved=resolved, singular_queries=singular,
                            nonsingular_unresolved=0, wrong_accepted_signs=0, singular_accepted=0)
            assert report["summary"] == expected and len(report["results"]) == count
            for row in report["results"]:
                if row["singular"]:
                    assert row["sign"] is None
                else:
                    assert row["sign"] == row["expected"] and row["residual"] <= 1e-10
            log = (native / f"delayed-native-verified-{name}.log").read_text()
            assert "runtime error:" not in log and "Assertion" not in log
            if args.library:
                output = work / (name + ".json")
                subprocess.run([sys.executable, str(work / "native/build-serial/pfaffian-stability" / script),
                                str(args.library.resolve()), str(output)], check=True, stdout=subprocess.PIPE)
                assert json.loads(output.read_text())["summary"] == expected
            synthetic.append(expected)
        before = json.loads((native / "delayed-native-pattern-chain-before.json").read_text())
        after = json.loads((native / "delayed-native-pattern-chain-after.json").read_text())
        assert len(before) == len(after) == 930
        failures = sum(row["status"] != 0 for row in before)
        assert failures == 435
        for a, b in zip(before, after, strict=True):
            assert (a["i"], a["j"]) == (b["i"], b["j"])
            assert b["status"] == 0 and b["sign"] == 1 and b["residual"] < 1e-10
        for ranks in (2, 4):
            log = (native / f"delayed-native-complete-unit-mpi{ranks}.log").read_text()
            assert "Sparse Pfaffian polynomial/permutation/scale/component tests passed" in log
            assert "ERROR STOP" not in log
        assert "Status: OK" in (native / "delayed-native-final-pretty.log").read_text()
        assert "Summary: correct: 74 / 74" in (native / "delayed-native-verified-regtests.log").read_text()
        assert "Summary: correct: 45 / 45" in (work / "native/build-serial/delayed-native-regtests.log").read_text()
        result = dict(verified_files=len(manifest["files"]), material_queries=len(controls), resolved=len(controls),
                      max_material_residual=max(row["residual"] for row in controls), synthetic=synthetic,
                      roundoff_queries=len(after), roundoff_failures_before=failures,
                      dense_comparisons_in_verified_reference=baseline["dense_comparisons"],
                      max_dense_gap_error_in_reference=baseline["max_dense_gap_error"], controls=controls,
                      recomputed_synthetic=bool(args.library))
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
        print(json.dumps({key: value for key, value in result.items() if key != "controls"}, indent=2))


if __name__ == "__main__":
    main()
