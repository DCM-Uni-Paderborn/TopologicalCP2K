# Periodic quadratic pseudospectrum validation

Recorded 25 September 2026. These are small numerical integration tests,
not converged localization predictions, a topological classifier, or band unfolding.

## Reproduce the comparisons

From the manuscript root, with Python 3 and NumPy:

```sh
python3 validation/development/summarize_periodic_quadratic.py \
  validation/development/periodic-quadratic-records.tar.gz \
  validation/development/periodic-quadratic-summary.json
```

This verifies all 171 member checksums and reads the raw CP2K output without
extracting files. It independently recomputes the comparisons and the decomposition
of each squared gap into energy and scaled chord residuals. It also checks the four
expected input rejections, final regression totals, and successful make_pretty.
The archive contains 2,987,919 uncompressed bytes, compressed to 668,105 bytes;
its SHA-256 is
`138a4bbd49dc16264daf1a5fcb3c3a865c95eb12ae54108fd6ebbfa2ca5d1ed8`.
The manifest records the tested source commit and hashes of the changed modules.
It is a source/evidence snapshot, not a complete standalone CP2K build distribution.

## Results and scope

- Seven matched serial-dense/two-rank-iterative inputs give fourteen gaps;
  maximum difference 4.997e-15 hartree, energy-residual difference 1.931e-13 hartree.
- The maximum iterative eigenpair residual is 6.666e-11 hartree squared, in Ne/SOC.
  Degenerate coefficient gauges are not inferred from eigenvalue agreement.
- Full/reduced SCF, query wrapping, equivalent MP/MACDONALD/GENERAL meshes and
  rigid skew-cell rotation agree within 1.1e-15 hartree. Axis permutations agree
  within 2.0e-15 hartree. These do not imply primitive-basis invariance of the
  reciprocal-coordinate chord metric.
- Primitive He and independently converged Gamma-supercell gaps differ by
  9.296e-10 hartree; energies per primitive cell differ by 9.170e-10 hartree.
  Explicit and implicit complex-Gamma setup paths differ by 2.151e-12 hartree.
- Across all successful runs the gap/residual decomposition agrees within
  4.49e-15 hartree squared.
- Official affected localizer/quadratic directories: 52/52 serial and 87/87 MPI
  assertions. The serial build intentionally skips MUMPS/Tacho-only directories.
  Existing finite-system reference values are unchanged.

The archived inputs specify UZH MOLOPT bases, GTH potentials, LDA/PADE,
300 Ry single grids, 1e-10 SCF tolerance, 0.1 hartree/bohr scaling, and two states.
They include periodicities X, Y, Z, XY, XZ, YZ and XYZ across the base/variant sets.
GPW and GAPW tests use pseudopotentials; the periodic set does not add a new
all-electron GAPW benchmark. SOC is post-SCF GTH on a restricted scalar potential.
The operator assembly uses the complete property mesh, even when SCF is reduced.
Open directions are not localized by the periodic formulation.

## Build and test details

GNU Fortran 16.1, OpenMPI 5.0.9, OpenBLAS, local MUMPS and DBCSR builds;
one BLAS thread, two OpenMP threads per process. The official driver uses
`--maxtasks 1`, `--mpiranks 1` for SSMP or `--mpiranks 2` for PSMP, and
`--ompthreads 2`. The final logs retain the exact selected directories and outcomes.
The serial smoke runner and variant generator are included in the archive.
The initial implicit-Gamma supercell run is retained as independent evidence;
the registered supercell regression uses explicit complex Gamma to avoid the
expensive Gamma-to-k-point conversion. No CP2K runtime optimization is inferred.

The seven-observable kernel and its native unit program were additionally
compiled with `-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fsanitize=undefined -fno-sanitize-recover=all -D__HAS_IEEE_EXCEPTIONS` and passed.
The initial custom build omitted `__HAS_IEEE_EXCEPTIONS` and trapped while running
the native test; the passing instrumented build enables the solver's optional
IEEE exception handling around LAPACK. The release build does not enable these
traps. Both logs are retained. This was a test-harness build
configuration correction, not an adjustment to physical references or tolerances.

The analysis is not an eigenvector-subspace comparison for the periodic DFT
cases, not a basis/cutoff/mesh convergence study, not a mixed open/periodic
position operator, and not evidence for multi-node scalability.
