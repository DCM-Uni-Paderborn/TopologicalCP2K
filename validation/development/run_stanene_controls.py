"""Separate basis, grid and cell-height effects in full-band SOC localizers.

Uses the independently validated complete-band research reference, not a new
production CP2K localizer implementation. Each analysis runs in a fresh process
so sparse factors are released before the next case. Existing evidence is never
overwritten and is reused only after checking its inputs and provenance.
"""

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, data):
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def cases():
    controls = []
    for basis, cutoff in (("DZVP", 200), ("DZVP", 400), ("DZVP", 600),
                          ("TZVP", 400), ("TZVP", 600)):
        controls.append(dict(label=f"{basis.lower()}-{cutoff}-scf8", basis=basis,
                             cutoff=cutoff, scf_mesh=8))
    for scf in (12, 18):
        controls.append(dict(label=f"tzvp-600-scf{scf}", basis="TZVP",
                             cutoff=600, scf_mesh=scf))
    controls += [dict(label="tzvp-600-rel60", basis="TZVP", cutoff=600,
                      scf_mesh=12, rel_cutoff=60),
                 dict(label="tzvp-600-height25", basis="TZVP", cutoff=600,
                      scf_mesh=12, height=25),
                 dict(label="dzvp-200-torus21", basis="DZVP", cutoff=200,
                      scf_mesh=8, size=21),
                 dict(label="tzvp-600-wilson", basis="TZVP", cutoff=600,
                      scf_mesh=12, mode="reference")]
    return [dict(dict(size=12, mode="bloch", height=20, rel_cutoff=40), **c)
            for c in controls]


def inputs(case, root):
    sys.path.insert(0, str(root / "build-serial"))
    from check_stanene_localizer import input_text, mesh_text

    text = input_text(SimpleNamespace(**case))
    # Keep the Cartesian buckling at 0.852 A when varying vacuum height.
    text = text.replace("C 0 0 20", f"C 0 0 {case['height']}")
    text = text.replace("0.5213", f"{0.5 + 0.426 / case['height']:.16g}")
    text = text.replace("0.4787", f"{0.5 - 0.426 / case['height']:.16g}")
    text = text.replace("REL_CUTOFF 40", f"REL_CUTOFF {case['rel_cutoff']}")
    result = {"input.inp": text}
    if case["mode"] == "bloch":
        import numpy as np

        mesh = mesh_text(case["size"])
        old_cell = "0 0 20\n"
        old_reciprocal = f"0 0 {2 * np.pi / 20:.16g}\n"
        assert mesh.count(old_cell) == mesh.count(old_reciprocal) == 1
        mesh = mesh.replace(old_cell, f"0 0 {case['height']}\n")
        mesh = mesh.replace(old_reciprocal, f"0 0 {2 * np.pi / case['height']:.16g}\n")
        result["mesh.nnkp"] = mesh
    return result


def provenance(root):
    paths = ["build-mpi/bin/cp2k.psmp", "build-mpi/src/libcp2k.2026.2.dylib",
             "data/BASIS_MOLOPT_UZH", "data/GTH_SOC_POTENTIALS"]
    paths += ["build-serial/" + name for name in ("check_stanene_localizer.py",
              "check_sparse_bloch_localizer.py", "check_bloch_localizer.py", "gaussian_states.py")]
    return dict(source_commit=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        sources={name: digest(root / name) for name in paths},
        runner_sha256=digest(Path(__file__)),
        environment={name: os.environ.get(name) for name in
                     ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "DYLD_LIBRARY_PATH")})


def export(case, work, root, stamp):
    expected = inputs(case, root)
    if (work / "run.json").exists():
        record = json.loads((work / "run.json").read_text())
        assert record["case"] == case and record["provenance"] == stamp
        assert record["returncode"] == 0 and record["scf_converged"] and record["completed"]
        for name, content in expected.items():
            assert (work / name).read_text() == content
        for name, value in record["outputs"].items():
            assert digest(work / name) == value
        return
    work.mkdir(parents=True, exist_ok=False)
    for name, content in expected.items():
        with (work / name).open("x") as stream:
            stream.write(content)
    command = ["mpiexec", "-n", "2", str(root / "build-mpi/bin/cp2k.psmp"),
               "-i", "input.inp", "-o", "output.out"]
    started = time.monotonic()
    env = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OMP_NUM_THREADS="2", OMP_STACKSIZE="64M")
    with (work / "stdout.log").open("x") as log:
        result = subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
    text = (work / "output.out").read_text()
    outputs = ["output.out", "stdout.log"]
    if result.returncode == 0:
        outputs += (["stanene.topology", "stanene.mmn", "stanene.eig"]
                    if case["mode"] == "bloch" else ["stanene.wilson"])
    record = dict(case=case, provenance=stamp, command=command, returncode=result.returncode,
                  elapsed_seconds=time.monotonic() - started,
                  scf_converged="*** SCF run converged" in text,
                  completed="PROGRAM ENDED AT" in text,
                  total_energy_hartree=re.findall(r"ENERGY\|.*?([-\d.]+)\s*$", text, re.M),
                  outputs={name: digest(work / name) for name in outputs})
    save(work / "run.json", record)
    if result.returncode or not record["scf_converged"] or not record["completed"]:
        raise RuntimeError(f"CP2K did not converge/finish: {work}")


def analyze_one(work, root, stamp):
    import numpy as np

    sys.path.insert(0, str(root / "build-serial"))
    from gaussian_states import read_snapshot
    from check_sparse_bloch_localizer import analyze

    snapshot = read_snapshot(work / "stanene.topology")
    valence = float(np.max(snapshot.energies[:, 7]))
    conduction = float(np.min(snapshot.energies[:, 8]))
    if not np.isfinite([valence, conduction]).all() or conduction - valence <= 1e-8:
        raise RuntimeError("No resolved sampled global electronic gap")
    energy = (valence + conduction) / 2
    report = analyze(work, energy, [0.75, 1., 1.25], 1.,
                     ctypes.CDLL(str(root / "build-mpi/src/libcp2k.2026.2.dylib")))
    report.update(inputs={name: digest(work / name) for name in ("stanene.topology", "stanene.mmn")},
                  provenance=stamp, sampled_indirect_gap_hartree=conduction-valence,
                  maximum_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                  * (1 if sys.platform == "darwin" else 1024))
    save(work / "analysis.json", report)
    if any(q["nu"] is None for q in report["queries"]):
        raise RuntimeError("Unresolved Pfaffian; retain the failed result")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--analyze-one", action="store_true")
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    stamp = provenance(root)
    if args.analyze_one:
        analyze_one(output, root, stamp)
        return
    selected = [c for c in cases() if args.cases is None or c["label"] in args.cases]
    if not selected or (args.cases and len(selected) != len(set(args.cases))):
        raise ValueError("Unknown or missing case")
    for case in selected:
        if shutil.disk_usage(output.parent).free < 2 * 1024**3:
            raise RuntimeError("Less than 2 GiB free; no new calculation started")
        work = output / case["label"]
        print(json.dumps({"start_case": case}), flush=True)
        export(case, work, root, stamp)
        if case["mode"] == "bloch":
            if (work / "analysis.json").exists():
                report = json.loads((work / "analysis.json").read_text())
                assert report["provenance"] == stamp
                for name, value in report["inputs"].items():
                    assert digest(work / name) == value
                assert all(q["nu"] is not None for q in report["queries"])
            else:
                with (work / "analysis.log").open("x") as log:
                    subprocess.run([sys.executable, __file__, str(root), str(work), "--analyze-one"],
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
            report = json.loads((work / "analysis.json").read_text())
            print(json.dumps({"finished": case["label"],
                              "electronic_gap_hartree": report["sampled_indirect_gap_hartree"],
                              "peak_gib": report["maximum_resident_bytes"] / 1024**3,
                              "gaps": [q["gap"] for q in report["queries"]],
                              "indices": [q["nu"] for q in report["queries"]]}), flush=True)
        else:
            print(json.dumps({"finished": case["label"]}), flush=True)


if __name__ == "__main__":
    main()
