# Periodized Gaussian translation checks

Local verification on 25 September 2026. This is a correctness record for
full-AO periodic translation analysis, not a converged material-unfolding or
parallel-scaling benchmark.

## Reference and Scope

The primitive Gamma-centered 3x1x1 mesh and explicit three-cell Gamma supercell
are converged independently. A separate supercell containing shifted ghost
bases supplies the original/ghost cross block of the ordinary overlap matrix.
It does not use the new translation matrix assembly. The original Gamma H and S
and the ghost-overlap T define an independent NumPy normal Gram. The Gaussian
integral kernel remains shared; the earlier finite archive contains an analytic
single-Gaussian check of that convention.

Cells, input parameters and quantitative results are in the main manuscript's
periodic Gaussian translation subsection and the corresponding SI subsection.
All five ghost-reference runs use three low states in an 18-AO torus. A separate
lattice-translation control checks the complete 18-state spectrum with T=S P,
where P cyclically permutes cell copies. Zero translation, commensurate wrapping,
rigid rotation, opposite momenta and reduced-SCF controls supplement this check.
Regression fixtures additionally cover a one-dimensional He chain and
two-dimensionally periodic Ne with 68 physical-spin SOC coefficients.

The initial skew-GAPW MPI calculation is retained in
`build-mpi/periodic-translation-skew` as a **pre-fix failure**, not a passing case.
Its onsite hard/soft density coefficients were not transferred by the existing
Gamma-to-k-point setup. The corrected calculations are named `gapw-checked`.
No independent spectrum or tolerance was adjusted to hide the failure.

## Replay

`torus-translation-records.tar.gz` contains 635 files (5,770,165 compressed
bytes). SHA-256:

```text
a7ed4c91a4674dae4351366c9fe944acd1311ccbc340e1f46463d51e17fbcbee
```

```sh
python3 validation/development/replay_torus_translation.py \
  validation/development/torus-translation-records.tar.gz
```

Requires Python and NumPy, not CP2K. The script checks every archived file hash,
recomputes independent spectra and true residuals from the original AO matrices
and complex states, and compares all new dense/iterative regression gaps.
The archive includes source snapshots, basis/potential files, input fixtures,
matrix/state output, formatter/build/test logs and a SHA-256 manifest. It omits
binaries, wavefunction restarts and large build products. Its manifest records
the local implementation checkpoint; this does not signify an upstream release.
Earlier finite and model-translation archives are retained unchanged.

The successful replay is `torus-translation-replay.log`. It also reproduces the
pre-fix GAPW failure deliberately, checking its spectrum error and true residual.
All 18 new regression gaps agree between serial dense and two-rank iterative
calculations within 2.8e-14 hartree; the largest iterative residual is 9.2e-11
hartree squared. The formatted/rebuilt broader quadratic, localizer, Kubo and
Wannier90 subset passes 163/163 serial and 210/210 MPI assertions. Those driver
counts include numerical unit-test executables as individual tests, not every
internal assertion. No pre-existing numerical reference was changed.

To regenerate a new independent run in a checkout with the requisite builds:

```sh
python3 build-serial/check_periodic_translation.py fresh-cubic-serial
python3 build-serial/check_periodic_translation.py fresh-skew-gapw-mpi \
  --skew --method GAPW --ranks 4 --solver ITERATIVE
python3 build-serial/check_translation_covariance.py fresh-covariance \
  --reference fresh-cubic-serial
```

The existing local environment uses GNU Fortran 16.1, OpenMPI 5.0.9,
thread-safe OpenBLAS 0.3.33, DBCSR 2.10.0 and MUMPS 5.9.1. Runs use two OpenMP
threads and one BLAS thread per process. The periodic assembly is a full-torus
reference implementation. It is not a linear-scaling or multi-node claim.

Arbitrary translated Gaussians can leave the AO span. The exact continuum
normal Gram retains this leakage; the energy term still uses the projected AO
Hamiltonian. The gap is not a spectral weight or a topological invariant.
Defect/disorder studies, basis and physical finite-size convergence, translation
selection and boundary suppression remain separate scientific work.
