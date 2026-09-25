"""Verify nullspace-probe evidence; optionally repeat all 600 native checks."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

from replay_stanene_convergence import digest, unpack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--library", type=Path)
    args = parser.parse_args()
    base = args.bundle.resolve()
    index = json.loads((base / "index.json").read_text())
    previous = base.parent / "pfaffian-pivot/index.json"
    assert digest(previous) == index["previous_index_sha256"]
    previous_index = json.loads(previous.read_text())
    material_index_path = base.parent / "stanene-convergence/index.json"
    assert digest(material_index_path) == previous_index["original_index_sha256"]
    material_index = json.loads(material_index_path.read_text())
    reference_archive = base.parent / "stanene-flattening.tar.gz"
    assert digest(reference_archive) == material_index["previous_reference_sha256"]
    with tempfile.TemporaryDirectory(prefix="pfaffian-probes-replay-") as temporary:
        work = Path(temporary)
        manifest = unpack(base / index["archive"]["file"], work, index["archive"]["sha256"])
        evidence = work / "build-serial/pfaffian-stability"
        before = json.loads((evidence / "production-probes-before.json").read_text())
        after = json.loads((evidence / "production-probes-after.json").read_text())
        prototype = json.loads((evidence / "adaptive-checked-independent.json").read_text())
        assert before["summary"]["singular_accepted"] == 72
        assert after["summary"] == dict(queries=600, resolved=450, singular_queries=132,
                                       nonsingular_unresolved=18, wrong_accepted_signs=0, singular_accepted=0)
        assert prototype["summary"] == dict(queries=600, resolved=468, singular_queries=132,
                                           nonsingular_unresolved=0, wrong_accepted_signs=0, singular_accepted=0)
        for old, new, trial in zip(before["results"], after["results"], prototype["results"], strict=True):
            for key in ("case", "scale", "swap", "singular", "expected"):
                assert old[key] == new[key] == trial[key]
            if not new["singular"]:
                assert old["sign"] == new["sign"]
                assert trial["sign"] == trial["expected"]
            else:
                assert new["sign"] is None and trial["sign"] is None
        fixture = json.loads((evidence / "singular-probe-fixture.json").read_text())
        n, seed = fixture["order"], fixture["seed"]
        assert (n, seed) == (24, 3)
        lower, blocks = np.eye(n), np.zeros((n, n))
        for j in range(n):
            for i in range(j + 1, min(n, j + 5)):
                lower[i, j] = ((17 * (i + 1) + 13 * (j + 1) + seed) % 5 - 2) / 8
        for i in range(0, n - 2, 2):
            blocks[i, i + 1] = (-1.0)**(i // 2)
            blocks[i + 1, i] = -blocks[i, i + 1]
        matrix = lower @ blocks @ lower.T
        assert np.array_equal(matrix, fixture["matrix"])
        assert np.array_equal(64 * matrix, np.rint(64 * matrix))
        assert np.linalg.matrix_rank(matrix) == 22 and np.all(np.any(matrix != 0, axis=1))
        native = json.loads((evidence / "independent-probes-mesh3.json").read_text())
        assert len(native["queries"]) == 3
        for row in native["queries"]:
            assert row["nu"] == row["dense_nu"] and row["dense_gap_error"] < 1e-10
            assert row["factor_residual"] < 1e-10
        material = json.loads((evidence / "adaptive-checked-mesh6.json").read_text())
        query, = material["queries"]
        assert query["nu"] == 1 and query["eta"] == 0.5 and query["factor_residual"] < 1e-10
        assert query["pfaffian_fill"] == 10951808
        assert abs(query["gap"] - 0.03623404693181014) < 1e-10
        with tarfile.open(reference_archive) as archive:
            prefix = "build-mpi/stanene-flattening/mesh6/"
            dense = json.load(archive.extractfile(prefix + "flattened-scan.json"))
            row, = [r for r in dense["queries"] if r["eta"] == query["eta"]]
            assert row["nu"] == query["nu"] and abs(row["gap"] - query["gap"]) < 1e-10
            for name, expected in material["inputs"].items():
                with archive.extractfile(prefix + name) as stream:
                    assert hashlib.file_digest(stream, "sha256").hexdigest() == expected
        assert material["library_sha256"] == index["prototype_library_sha256"]
        assert native["library_sha256"] == index["production_library_sha256"]
        for name, total in (("independent-probes-regtests.log", 45), ("independent-probes-sparse-regtests.log", 29)):
            assert f"Summary: correct: {total} / {total}" in (evidence / name).read_text()
        if args.library:
            subprocess.run([sys.executable, str(evidence / "validate_adaptive.py"),
                            str(args.library.resolve()), str(work / "recomputed.json")], check=True)
            recomputed = json.loads((work / "recomputed.json").read_text())
            assert recomputed["summary"] == after["summary"]
            for old, new in zip(after["results"], recomputed["results"], strict=True):
                assert (old["case"], old["scale"], old["swap"], old["sign"]) == (
                    new["case"], new["scale"], new["swap"], new["sign"])
        result = dict(verified_files=len(manifest["files"]), production=after["summary"],
                      prototype=prototype["summary"], native_material_queries=3,
                      native_regtest_checks=74, recomputed=bool(args.library))
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
