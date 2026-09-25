# Silicon translation validation

This record supports the periodic Gaussian-translation section and SI S7.9/Table S12.
It is an algebraic/material control, not a momentum-, cutoff-, basis-, or defect-converged
unfolding study.

## Completed comparisons

- Two-atom diamond Si, PBE/GPW, UZH TZVP/GTH-q4, a Gamma-centered 3 x 1 x 1 mesh.
- 102 AOs in the full three-cell property space; five wave-vector and seven energy
  queries, three states per query.
- Ordinary CP2K band energies predict the quadratic eigenvalues within
  1.21e-10 Ha^2 (serial dense) and 1.17e-10 Ha^2 (two-rank iterative).
- Dense/iterative gaps differ by at most 6.80e-12 Ha.
- An independently converged six-atom Gamma supercell agrees in its gaps within
  5.90e-12 Ha. A NumPy normal Gram from ordinary Gamma H/S and a cyclic cell
  permutation reproduces all 105 eigenvalues within 1.42e-12 Ha^2.
- The independently recomputed supercell Ritz residual is at most 9.83e-11 Ha^2;
  the largest S-orthogonality error is 7.21e-14.
- Fourier blocks of the ordinary supercell H/S reproduce native primitive band
  energies within 1.84e-10 Ha, at the precision of the printed band output.

## Screening controls

EPS_DEFAULT=1e-12 inherits EPS_PPNL=1e-8. The primitive matrix has a relative
Hermiticity defect of 4.11425e-8 and correctly fails the unchanged 1e-10 guard.
Tightening EPS_GVG_RSPACE, then EPS_PGF_ORB, does not remove it. With both at
1e-12 and EPS_PPNL=1e-14, the defect is below 4.94e-14. No torus-phase changes,
matrix averaging, tolerance relaxation, or reference-value updates were made.

Temporary dense dumps before the original guard abort were used to inspect the
failed matrices. Those raw matrices and their input/output records are archived;
this diagnostic gather is not part of the committed implementation. The production
change only exposes the distributed relative residual in output and on failure.
Its unit checks cover Hermitian, non-Hermitian and nonfinite matrices, including a
four-rank column distribution. Focused suites pass 30/30 serial and 54/54 MPI tests.

## Reproduce

The archive contains 212 checksummed files, including source snapshots, input
generation/analysis scripts, original logs, printed H/S, bands and complex states.
The companion JSON records the archive SHA-256. No build products or large restart
checkpoints are included. Successful runs reused earlier converged densities;
the input generator also supports fresh atomic guesses.

Run the independent replay with Python and NumPy:

```sh
python3 validation/development/replay_silicon_translation.py \
  validation/development/silicon-translation-records.tar.gz
```

The replay verifies all hashes, reconstructs the band and Gamma-matrix references,
compares serial/MPI gaps, and checks that the coarse-screening controls failed.
The input generator in the archive can also run fresh calculations in a matching
CP2K checkout with `build-serial/bin/cp2k.ssmp` and `build-mpi/bin/cp2k.psmp`:

```sh
python3 build-serial/check_silicon_translation.py si-serial --pristine-only
python3 build-serial/check_silicon_translation.py si-mpi --pristine-only \
  --solver ITERATIVE --ranks 2
```

Defective-Si/shifted-ghost cases in the generator are future test scaffolding and
are not included as completed results in this record.
