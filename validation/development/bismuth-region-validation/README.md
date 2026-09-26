# Joint finite-Bi parameter regions

This bundle extends the fixed-position scale intervals to connected boxes
in energy, scale and both in-plane positions. Each Hamiltonian and its
physical orthonormal AO metric remain frozen. It does not infer a bulk
classification or recompute the SCF potential at each query.

All boxes use E = neutral midpoint +/- 0.02 times the full neutral gap,
x/y = atomic centroid +/- 0.05 Angstrom, and the following scale ranges:

| Case | Kappa (hartree/bohr) | Index | Evaluations/leaves | Minimum lower bound (hartree) |
| --- | --- | ---: | ---: | ---: |
| 16 Bi, DZVP | 0.003 to 0.0035 | 1 | 9/5 | 2.5285740685e-4 |
| 16 Bi, TZVP | 0.003 to 0.0035 | 1 | 49/25 | 6.0315449051e-6 |
| 30 Bi, DZVP | 0.003 to 0.0035 | 0 | 15/8 | 2.5985717088e-5 |
| 30 Bi, DZVP | 0.00025 to 0.0004 | 1 | 21/11 | 4.9742893333e-6 |
| 30 Bi, DZVP | 0.0015 to 0.0025 | 1 | 19/10 | 2.2941118272e-4 |

The covered boxes have no unresolved leaves. The minimum leaf lower bound
is not the exact minimum spectral gap. The 1e-10 hartree margin accounts
for numerical evaluation only; this is not interval arithmetic.
All 16 corners plus the anchor of every box are checked with a complete
divide-and-conquer spectrum and Pfaffian. These checks also compare the
middle-four eigenvalue gap, spectral +/- pairing and actual perturbation
norm against the bound. The cover, not the finite corner sampling, connects
the index throughout each box.

Ten tests exercise independent Hermitian-matrix perturbations, exact scalar
limits, analytic gapped eigenvalues, origin translation, partition coverage,
budget/gap-closing failures, invalid inputs and exact discrete replay fields.

`methods-results.tar.gz` contains methods, reports, the underlying scans and
test log. `index.json` names checksum-pinned native snapshot archives in
the adjacent `bismuth-size-basis-validation` directory. Those retained
native inputs/outputs and complete spinor exports are not duplicated here.
No native DFT calculation or new native AO assembly is claimed by replay.

From the manuscript repository root, the full numerical replay is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  validation/development/archive_bismuth_regions.py . . \
  validation/development/bismuth-region-validation --replay \
  --output .build/bismuth-regions-replayed.json
```

The command verifies checksums and all report fields except elapsed time,
and refuses to overwrite the chosen output. Temporary extraction uses
fresh timestamps in the checkout's ignored `.build` directory. The
retained `replay.json` checks 72 archive members and reproduces all five
complete region reports, including the 85 independent corner/anchor
audits, with zero numerical difference. It does not rerun native SCF.
