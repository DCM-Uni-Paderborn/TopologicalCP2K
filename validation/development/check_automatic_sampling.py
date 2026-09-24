"""Local differential and rejection checks for automatic reciprocal sampling."""
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np

root = Path(__file__).resolve().parents[2]
work = Path(__file__).resolve().parent
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2", OMP_STACKSIZE="128M",
           CP2K_DATA_DIR=str(root / "data"),
           DYLD_INSERT_LIBRARIES=str(root / "build-serial/openblas-thread-safe/lib/libopenblas.dylib"))
binary = root / "build-serial/bin/cp2k.ssmp"
results = []

def run(name, source=None, failure=None):
    if source is not None:
        (work / (name + ".inp")).write_text(source)
    with (work / (name + ".launch.log")).open("w") as output:
        process = subprocess.run([binary, "-i", name + ".inp", "-o", name + ".out"],
                                 cwd=work, env=env, stdout=output, stderr=subprocess.STDOUT, timeout=180)
    text = (work / (name + ".out")).read_text()
    if failure:
        normalized = " ".join(re.sub(r"[|*]", " ", text).split())
        assert process.returncode != 0 and failure in normalized, (name, process.returncode, text[-2000:])
        print(name, "expected rejection:", failure, flush=True)
        results.append({"case": name, "expected_rejection": failure})
    else:
        assert process.returncode == 0 and "PROGRAM ENDED AT" in text, (name, text[-2000:])
    return text

def expanded(name):
    source = (root / "tests/QS/regtest-little-group" / (name + ".inp")).read_text()
    return source.replace("@INCLUDE automatic-common.inc", (work / "automatic-common.inc").read_text())

def blocks(path, start, end):
    result = {}
    current = None
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == start:
            current = tuple(fields[1:])
            result[current] = []
        elif fields[0] == end:
            current = None
        elif current is not None:
            result[current].append(fields)
    return result

for name in ("automatic-helium", "automatic-gapw", "automatic-soc", "automatic-shear"):
    text = run(name)
    assert not (work / (name + ".mmn")).exists(), "Native analysis wrote a misleading MMN"
    detail = work / (name + ".little_group")
    records = [line.split() for line in detail.read_text().splitlines()]
    points = [list(map(float, row[2:])) for row in records if row and row[0] == "AUTO_POINT"]
    edges = defaultdict(list)
    for row in records:
        if row and row[0] == "AUTO_EDGE":
            source, *target = map(int, row[2:])
            edges[source].append(target)
    degree = max(map(len, edges.values()))
    shear = 4.0 if name == "automatic-shear" else 0.0
    lattice = np.array([[4., 0., 0.], [0., 5., 0.], [shear, 0., 8.]])
    nnkp = []
    for label, matrix in (("real_lattice", lattice), ("recip_lattice", 2*np.pi*np.linalg.inv(lattice).T)):
        nnkp += ["begin " + label, *[" ".join(map(str, row)) for row in matrix], "end " + label]
    nnkp += ["begin kpoints", str(len(points)), *[" ".join(map(str, row)) for row in points], "end kpoints"]
    nnkp += ["begin nnkpts", str(degree)]
    for point in range(1, len(points) + 1):
        padded = list(edges[point])
        for target in range(1, len(points) + 1):
            candidate = [target, 0, 0, 0]
            if len(padded) == degree:
                break
            if candidate not in padded:
                padded.append(candidate)
        nnkp.extend(" ".join(map(str, [point, *target])) for target in padded)
    nnkp += ["end nnkpts"]
    replay = name + "-replay"
    (work / (replay + ".nnkp")).write_text("\n".join(nnkp) + "\n")
    source = expanded(name).replace("@SET PROJECT " + name, "@SET PROJECT " + replay)
    source = source.replace("KPOINTS_SOURCE SYMMETRY", "KPOINTS_SOURCE NNKP\n        LITTLE_GROUP_COMPATIBILITY T\n        NNKP_FILE " + replay + ".nnkp")
    source = source.replace("REQUIRE_GLOBAL_GAP T", "REQUIRE_GLOBAL_GAP F")
    run(replay, source)
    eigenvalues = np.loadtxt(work / (name + ".eig"))
    reference = np.loadtxt(work / (replay + ".eig"))
    assert np.array_equal(eigenvalues[:, :2], reference[:, :2])
    error = float(np.max(np.abs(eigenvalues[:, 2] - reference[:, 2])))
    assert error < 1e-9, (name, error)
    left = blocks(detail, "KPOINT", "END_KPOINT")
    right = blocks(work / (replay + ".little_group"), "KPOINT", "END_KPOINT")
    assert left.keys() == right.keys()
    character_error = 0.0
    for kpoint in left:
        for tag in ("OPERATION", "SIZES", "IRREP", "COREP"):
            assert [r for r in left[kpoint] if r[0] == tag] == [r for r in right[kpoint] if r[0] == tag]
        a = np.array([list(map(float, r[1:])) for r in left[kpoint] if r[0] == "BAND_CHARACTER"])
        b = np.array([list(map(float, r[1:])) for r in right[kpoint] if r[0] == "BAND_CHARACTER"])
        character_error = max(character_error, float(np.max(np.abs(a-b))))
    assert character_error < 1e-9
    native_edges = blocks(detail, "EDGE", "END_EDGE")
    replay_edges = blocks(work / (replay + ".little_group"), "EDGE", "END_EDGE")
    assert native_edges.keys() <= replay_edges.keys()
    assert all(value == replay_edges[key] for key, value in native_edges.items())
    atomic = detail.read_text().split("ATOMIC_SIGNATURE_SIZES ", 1)[1]
    atomic_replay = (work / (replay + ".little_group")).read_text().split("ATOMIC_SIGNATURE_SIZES ", 1)[1]
    assert atomic == atomic_replay
    result = dict(case=name, points=len(points), edges=sum(map(len, edges.values())), eigenvalue_error_ev=error,
                  character_error=character_error, exact_compatibility_and_atomic_witnesses=True)
    print(json.dumps(result), flush=True)
    results.append(result)

base = expanded("automatic-helium")
for name, old, new, rejection in (
    ("budget", "KPOINTS_SOURCE SYMMETRY", "KPOINTS_SOURCE SYMMETRY\n LITTLE_GROUP_MAX_POINTS 2", "exceeds LITTLE_GROUP_MAX_POINTS"),
    ("path-count", "KPOINTS_SOURCE SYMMETRY", "KPOINTS_SOURCE SYMMETRY\n LITTLE_GROUP_PATH_POINTS 1", "Invalid automatic little-group sampling limits"),
    ("wilson", "KPOINTS_SOURCE SYMMETRY", "KPOINTS_SOURCE SYMMETRY\n WILSON_LOOP T", "SYMMETRY segments are not Wilson loops"),
    ("gap", "@SET EXCLUDED 3", "@SET EXCLUDED 2 3", "Band symmetry needs isolated selected bands"),
):
    run(name, base.replace("@SET PROJECT automatic-helium", "@SET PROJECT " + name).replace(old, new), rejection)
for name, old, new in (
    ("refined", "KPOINTS_SOURCE SYMMETRY", "KPOINTS_SOURCE SYMMETRY\n LITTLE_GROUP_PATH_POINTS 5"),
    ("unitary", "TIME_REVERSAL T", "TIME_REVERSAL F"),
):
    text = run(name, base.replace("@SET PROJECT automatic-helium", "@SET PROJECT " + name).replace(old, new))
    assert "TQC| Incompatible segments: 0" in text
    assert "TQC| Sampled nonnegative atomic membership (1/0/-1): 1" in text
    if name == "unitary":
        assert "TQC| Antiunitary compatibility segments: 0" in text
    summary = [line.strip() for line in text.splitlines() if "Automatic" in line or "sampled" in line and "gap" in line]
    print(name, summary, flush=True)
    results.append(dict(case=name, summary=summary))
(work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
