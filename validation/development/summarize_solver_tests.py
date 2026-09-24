"""Recheck retained numerical records; does not launch CP2K or fit references."""

import hashlib
import json
import re
import tarfile
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "solver-test-records.tar.gz"


def values(text, section, label):
    prefix = f"{section}| {label}:"
    return np.array([
        [float(value) for value in line.split(prefix, 1)[1].split()]
        for line in text.splitlines() if prefix in line
    ])


def scalar(text, section, label):
    result = values(text, section, label)
    assert result.shape == (1, 1), (section, label, result.shape)
    return float(result[0, 0])


def solver_record(text, label):
    marker = f"{label} order/time(dense,sparse)/gap/nnz/factor/work/count/MB"
    rows = [line.split(marker, 1)[1].split() for line in text.splitlines()
            if line.startswith(marker)]
    assert len(rows) == 1 and len(rows[0]) == 11, (label, rows)
    order, dense, sparse, gap, lower, upper, nnz, fill, work, count, memory = rows[0]
    gap, lower, upper = map(float, (gap, lower, upper))
    assert lower - 1e-9 <= gap <= upper + 1e-9, (label, gap, lower, upper)
    return dict(realified_order=int(order), dense_s=float(dense), sparse_s=float(sparse),
                gap=gap, lower=lower, upper=upper, input_entries=int(nnz),
                factor_entries=int(fill), workspace_entries=int(work),
                factorization_count=int(count), solver_memory_sum_mb=float(memory),
                factor_input_ratio=int(fill) / int(nnz))


def overlap_matrix(text, n):
    block = text.split(" OVERLAP MATRIX")[-1]
    result = np.full((n, n), np.nan)
    columns = []
    for line in block.splitlines():
        fields = line.split()
        if fields and all(field.isdigit() for field in fields):
            columns = [int(field) - 1 for field in fields]
        elif columns and len(fields) == 4 + len(columns) and fields[0].isdigit():
            row = int(fields[0]) - 1
            result[row, columns] = [float(field) for field in fields[4:]]
            if row == n - 1 and columns[-1] == n - 1:
                break
    assert np.isfinite(result).all()
    assert np.max(np.abs(result - result.T)) < 1e-14
    np.linalg.cholesky(result)
    return result


def coefficients(text):
    columns = []
    for block in text.split("# State ")[1:]:
        rows = [line.split() for line in block.splitlines()[1:]
                if line and not line.startswith("#")]
        columns.append([complex(float(row[1]), float(row[2])) for row in rows])
    result = np.array(columns).T
    assert np.isfinite(result).all()
    return result


def main():
    with tarfile.open(ARCHIVE) as archive:
        records = {}
        hashes = {}
        for member in archive.getmembers():
            if member.isfile():
                raw = archive.extractfile(member).read()
                records[member.name] = raw.decode("utf-8")
                hashes[member.name] = hashlib.sha256(raw).hexdigest()

    qwz = []
    for length in (10, 20, 30):
        for ranks in (1, 2, 4):
            text = records[f"build-mpi/localizer-final-{length}-{ranks}.log"]
            assert "ERROR STOP" not in text and "status(dense,sparse)" not in text
            row = solver_record(text, "QWZ")
            other = solver_record(text, "QWZ nonorthogonal")
            assert row["realified_order"] == 8 * length**2
            assert abs(row["gap"] - other["gap"]) < 1e-9
            qwz.append(dict(length=length, ranks=ranks, **row))

    pfaffian_text = records["build-mpi/tacho-unit-final-large.log"]
    assert "Sparse Pfaffian polynomial/permutation/scale/component tests passed" in pfaffian_text
    assert "ERROR STOP" not in pfaffian_text and "status(dense,sparse)" not in pfaffian_text
    pfaffian = {}
    for label in ("AII spin-mixed", "AII trivial QWZ pair", "AII complex nonorthogonal"):
        row = solver_record(pfaffian_text, label)
        line = next(line for line in pfaffian_text.splitlines()
                    if line.startswith(label + " Pfaffian residual/fill"))
        residual, fill = line.split("Pfaffian residual/fill")[1].split()
        row.update(solve_probe_residual=float(residual), symbolic_skew_blocks=int(fill))
        pfaffian[label] = row

    near_text = records["build-mpi/periodic-sparse-unit.log"]
    near = {label: solver_record(near_text, label) for label in
            ("gap not pivot size", "gap not pivot size nonorthogonal", "near closing", "unresolved", "singular")}

    quadratic = {}
    section = "QUADRATIC_PSEUDOSPECTRUM"
    for method in ("gpw", "gapw"):
        output = [records[f"build-mpi/quadratic-bismuth-{method}-{solver}.out"]
                  for solver in ("dense", "iterative")]
        assert all("PROGRAM ENDED AT" in text for text in output)
        n = int(scalar(output[0], section, "AO/spinor dimension"))
        metric = np.kron(np.eye(2), overlap_matrix(output[0], n // 2))
        gaps = [values(text, section, "Gap [hartree]").ravel() for text in output]
        eigenres = [values(text, section, "Eigenpair residual [hartree^2]").ravel()
                    for text in output]
        states = [coefficients(records[f"build-mpi/quadratic-bismuth-{method.upper()}-{solver}-1_0.states"])
                  for solver in ("DENSE", "ITERATIVE")]
        norm_errors = [float(np.max(np.abs(c.conj().T @ metric @ c - np.eye(4)))) for c in states]
        assert all(c.shape == (n, 4) for c in states)
        singular = np.linalg.svd(states[0].conj().T @ metric @ states[1], compute_uv=False)
        gap_error = float(np.max(np.abs(gaps[0] - gaps[1])))
        decomposition_errors = []
        for text, gap in zip(output, gaps):
            energy = values(text, section, "Energy residual [hartree]").ravel()
            position = values(text, section, "Position residuals [bohr]")
            kappa = scalar(text, section, "Kappa [hartree/bohr]")
            error = float(np.max(np.abs(gap**2 - energy**2 - kappa**2 * np.sum(position**2, axis=1))))
            assert error < 1e-12
            decomposition_errors.append(error)
        assert gap_error < 1e-14
        assert max(norm_errors) < 1e-10
        assert np.max(np.abs(singular - 1)) < 1e-10
        assert max(eigenres[1]) < 1e-9
        quadratic[method] = dict(dimension=n, dense_gaps_ha=gaps[0].tolist(),
                                iterative_gaps_ha=gaps[1].tolist(), max_gap_error_ha=gap_error,
                                max_eigenpair_residual_ha2=[float(max(v)) for v in eigenres],
                                metric_norm_errors=norm_errors, subspace_singular_values=singular.tolist(),
                                max_subspace_error=float(np.max(np.abs(singular - 1))),
                                residual_decomposition_errors_ha2=decomposition_errors,
                                iterative_steps=int(scalar(output[1], section, "Ritz iterations")))

    torus = []
    for case in ("kpoints", "supercell"):
        text = records[f"build-serial/torus-converged-{case}.out"]
        assert "PROGRAM ENDED AT" in text
        energy = float(re.search(r"ENERGY\| Total FORCE_EVAL.*?\[hartree\]\s+([-0-9.]+)", text)[1])
        torus.append(dict(case=case, scf_energy_ha=energy,
                          gap_ha=scalar(text, "SPECTRAL_LOCALIZER", "Gap [hartree]"),
                          z2=int(scalar(text, "SPECTRAL_LOCALIZER", "Z2 index")),
                          dimension=int(scalar(text, "SPECTRAL_LOCALIZER", "Localizer dimension"))))
    assert all(row["z2"] == 0 and row["dimension"] == 312 for row in torus)
    torus_gap_error = abs(torus[0]["gap_ha"] - torus[1]["gap_ha"])
    assert torus_gap_error < 5e-10
    result = dict(archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
                  source_sha256=hashes, qwz=qwz, pfaffian=pfaffian, near_gap=near,
                  quadratic=quadratic, periodic_torus=torus,
                  periodic_gap_error_ha=torus_gap_error,
                  periodic_energy_error_per_cell_ha=abs(torus[0]["scf_energy_ha"] - torus[1]["scf_energy_ha"] / 6))
    destination = ROOT / "solver-test-summary.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "source_sha256"}, indent=2))


if __name__ == "__main__":
    main()
