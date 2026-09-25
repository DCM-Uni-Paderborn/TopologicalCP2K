# Displaced silicon and frozen-potential validation

## Completed

- One Si atom displaced by 0.1 angstrom in the six-atom Gamma cell used by the
  earlier pristine silicon check. The 102-AO TZVP model, PBE, 400 Ry cutoff,
  tight projector/orbital screening and 35 energy/momentum queries are unchanged.
- Two MPI ranks, two OpenMP threads/rank, three quadratic eigenpairs/query,
  40-dimensional iterative trial space and 1e-10 Hartree-squared Ritz tolerance.
- A separate calculation with translated ghost bases provides an independent
  overlap reference for the full translation matrix. Coincident ghost/real centers
  require `CONNECTIVITY OFF`; no ghost electronic state is used as a reference.
- Independent generalized diagonalization reproduces all 105 low eigenvalues
  within 3.7277e-14 Hartree squared. The true Ritz residual is at most 9.7696e-11.
  Metric normalization, complex means, full/projected translation residuals,
  energy residuals and projection leakage are independently contracted from
  the original SCF matrices and native complex states.
- The distorted system has a nonzero Hamiltonian/translation commutator
  (0.351205 Hartree) and compression leakage operator norm (0.146770).
  Selected-state leakage is 0.00364918--0.06824526; the largest gap change from
  pristine Si is 0.01325324 Hartree. These are implementation controls, not
  converged unfolding intensities or defect predictions.

## Frozen-density correction

`create_kp_from_gamma` previously left the newly initialized kinetic-energy
density and GAPW_XC grids at their atomic guesses. He2/TZVP controls compare all
twelve Gamma band eigenvalues against the original printed SCF H/S matrices:

| Method | Maximum error before (Ha) | Maximum error after (Ha) |
|---|---:|---:|
| GPW/PBE | 1.55e-10 | 1.55e-10 |
| GPW/TPSS | 4.17e-3 | 1.82e-10 |
| GAPW/TPSS | 8.54e-3 | 1.63e-10 |
| GAPW_XC/PBE | 1.05e-2 | 1.83e-10 |
| GAPW_XC/TPSS | 1.56e-2 | 1.44e-10 |

After-errors cover serial and two-rank MPI runs and are limited by the printed
eight-decimal-eV band precision. Seven new regressions include explicit-Gamma
k-point counterparts that bypass the conversion. The new eigenvalue matcher
rejects all four affected before-controls. The focused suite passes 35/35 serial
and 58/58 MPI assertions; `make_pretty.sh --no-cache` passes all 13 changed files.
No existing reference or numerical acceptance threshold was relaxed.

## Evidence and replay

The 155-file `silicon-defect-frozen-density-records.tar.gz` archive contains inputs,
raw matrices and states, successful and failed density controls, test output,
source snapshots, scripts, basis/potential files and SHA-256 checksums. Its SHA-256:

`ee47ea30ace7122343dcf91db6d2dccecfb96aadb2f71a2b136e673ad15825ac`

Recompute without CP2K using Python 3.14 and NumPy:

```text
python3 replay_defect_density.py silicon-defect-frozen-density-records.tar.gz
```

The successful replay verifies all file hashes, reconstructs the matrix references,
checks the known before-failures and verifies the post-fix regression summaries.
It does not invoke an external Wannier or band-unfolding program.

## Remaining scope

These tests do not validate hybrid/nonlocal-exchange Gamma image densities or
implicit-solvent state transfer. Collocated XC gradients retain the existing
initialization guard. GAPW_XC band transfer is tested, but the quadratic adapter
still permits GPW/GAPW only. Larger defective-cell and basis/cutoff/momentum studies,
realistic converged SOC localizer classifications and multi-node scaling remain
distinct tasks. The code correction is locally committed; no CP2K push is implied.
