# Finite Gaussian translation checks

Validated locally on 25 September 2026. This record concerns finite Gaussian
Hamiltonian/translation analysis, not periodic material band unfolding.

## Evidence

`gaussian-translation-records.tar.gz` contains 501 files (1,717,359 compressed
bytes), including source snapshots, input fixtures, AO matrix output, complex
states, successful driver logs, input rejection logs, and a per-file SHA-256
manifest. Archive SHA-256:

```text
e91be13a10c0c8751f1aaca3d09fc6335e878cffcd2f8a61ec1ffc229c65228e
```

The manifest records the local implementation checkpoint. This is not an upstream
release or a CP2K push. Source snapshots reflect the final formatted implementation;
binary banners may identify its preceding checkpoint with local modifications.
Earlier translation-model and periodic-torus archives are retained separately.

Replay using Python with NumPy, without running CP2K:

```sh
python3 validation/development/replay_gaussian_translation.py \
  validation/development/gaussian-translation-records.tar.gz
```

The replay validates every file checksum before re-evaluating spectra, complex
means, metric orthogonality, true eigenpair residuals and projection leakage from
the raw AO matrices and states. `gaussian-translation-replay.log` is the successful
local replay. No wavefunction restarts, binaries or large build products are in
the archive.

## Independent Gaussian Reference

For each system, a separate calculation places ghost copies of its Gaussian basis
at R_B-a. The original/ghost cross block of the ordinary overlap matrix supplies
T independently of the new translation matrix assembly. The ghost calculation's
Hamiltonian is not used. Both paths share CP2K's existing Gaussian integral kernel;
the single-Gaussian overlap additionally has a closed-form reference.

Cases: one normalized He s Gaussian (exponent 0.7 bohr^-2), TZVP He/GPW (6 AOs),
TZVP H2/GPW (12 AOs), and TZVP Ne/GAPW (17 AOs, including d functions). The last
three use BASIS_MOLOPT_UZH and POTENTIAL_UZH. Query E=-0.1 hartree, translation
scale 0.3 hartree, a=(0.63,-0.37,0.21) angstrom, k=(0.31,-0.22,0.18) bohr^-1.
Translation convention: U(a) psi(r)=psi(r+a), with phase exp(i k.a).

Base, ghost, opposite-momentum and zero-translation runs give 16 serial dense and
16 two-rank iterative calculations. Largest squared-gap reference discrepancy:
6.0771e-12 hartree^2. Largest independently reconstructed eigenpair residual:
1.1202e-11 hartree^2. Maximal leakage norms range from 0.47708 to 0.87416, so they
cannot be neglected. Zero-translation squared residuals are at roundoff; taking
their square root gives at most 2.1074e-8.

The exact continuum translation normal Gram is
(1+|z|^2)S - conjugate(z)T - z T^dagger. Its difference from the projected normal
Gram is S - T^dagger S^-1 T. The energy term still uses the projected AO
Hamiltonian. No polar unitarization or boundary weighting is applied.

To regenerate these calculations with suitable builds, use the archived scripts:

```sh
python3 build-serial/gaussian-translation-validation/check_gaussian_translation.py \
  --root . --binary build-serial/bin/cp2k.ssmp --output fresh-translation-serial
python3 build-serial/gaussian-translation-validation/check_gaussian_translation.py \
  --root . --binary build-mpi/bin/cp2k.psmp --solver ITERATIVE --ranks 2 \
  --output fresh-translation-mpi
python3 build-serial/gaussian-translation-validation/check_translation_extensions.py \
  --root . --output fresh-translation-extensions
```

## Further Checks

- Dense versus four-rank iterative all-electron H2/GAPW gaps differ by at most
  7.50e-14 hartree; post-SCF GTH-SOC Ne gaps by 1.619e-10 hartree. These comparisons
  include independently converged SCF potentials, not only identical-matrix solves.
- Opposite momenta agree within 6.0e-15 hartree in these additional fixtures.
- Missing translations, POSITION/KAPPA conflicts, incompatible formulation and
  MP_GRID usage are rejected with the expected diagnostics. One preliminary
  checker expected different wording for the mesh error; the corrected replay
  verifies the actual intended diagnostic, not a change of numerical tolerance.
- A compressed three-state unitary has an exact two-state analytic reference and
  discarded norm |c_1|^2; its eigenvalues agree within 2.3e-16 in model units.
- The numerical kernel passes bounds checks, undefined-behavior checks and traps
  for invalid operations, division by zero and overflow. External LAPACK calls
  retain the existing narrowly scoped trap handling. This is not a wholly
  instrumented CP2K/MPI build.
- A false ~1e-6 identity-translation residual in an initial MPI run traced to
  cached metric-image drift after Ritz restarts. Final states are now normalized
  and verified with fresh S/Q applications. No physical reference values were
  adjusted to hide it.
- Final formatted/rebuilt drivers pass 57/57 serial and 96/96 MPI assertions,
  covering the affected quadratic/localizer subset and available MUMPS/Tacho
  paths. New GPW/GAPW inputs use two translations and two momenta, with shared
  dense/iterative gap and leakage references. Earlier references are unchanged.

The environment is the existing GNU Fortran 16.1/OpenMPI 5.0.9 build with
thread-safe OpenBLAS 0.3.33, DBCSR 2.10.0 and MUMPS 5.9.1. Runs use two OpenMP
threads per process and one BLAS thread. These are correctness checks on one host,
not performance or multi-node scaling measurements. Periodic unfolding,
primitive/supercell material comparisons, boundary suppression and physical
basis/finite-size convergence remain separate work.
