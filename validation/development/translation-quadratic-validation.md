# Complex translation queries: numerical validation

Date: 25 September 2026. These are model/kernel tests, not Gaussian DFT
unfolding results. The main text and SI report model energy units with t=1;
quadratic eigenvalues and Ritz residuals have units t squared.

## Evidence

`translation-quadratic-records.tar.gz` contains source fixtures, final build
logs, native and distributed output, formatting results, focused instrumented
output, and the affected local regression runs. Its manifest identifies the
source revision and records a size and SHA-256 for every archived file.
Executables, libraries, restart wavefunctions and regenerable build trees
are excluded. `translation-quadratic-summary.json` is derived from this archive,
not a manually transcribed reference dataset.
The archive contains 252 checksummed files and occupies 663,850 bytes.
Its SHA-256 is
`78d7837b3f28b8805c4d20eb48d5bb1f50239ead7a060fbcae87a0403b61525d`.

| Check | Observed maximum or outcome |
| --- | --- |
| 18-site trimerized ring versus six independent Bloch blocks | 1,170 eigenvalues per 65-query scan; 1.598721e-14 |
| Explicit LU/generalized-eigenvalue normal-product reference | 3.552714e-15 |
| Nonsingular complex basis change, spectrum | 5.329071e-15 |
| Complex expectations after basis change, scaled model units | 2.248202e-14 |
| Deliberately incorrect independent Hermitian squares | 8.604947e-2 spectral discrepancy |
| Distributed forward/adjoint products versus dense | 4.440892e-16 / 2.482534e-16 |
| Complex metric solves | 4.577567e-16 |
| Distributed iterative versus dense eigenvalues | 8.770762e-15 |
| Largest iterative eigenpair residual | 9.833999e-11 |
| Final local regression drivers | 53/53 serial, 88/88 two-rank MPI |

The periodic scalar chain also checks a known Bloch phase and energy residual.
The open nonnormal translation fixture compares complete spectra, the two-state
projector, arbitrary complex query shifts, and expectations. An intentionally
incorrect matrix-free adjoint is rejected. The archived Fortran fixtures fully
specify each Hamiltonian, translation, metric, query and solver setting.

One-, two- and four-rank MPI runs exercise row- and column-oriented DBCSR process
grids. Both reusable dense factors and MUMPS metric factors are checked. The
focused MPI runs use one OpenMP and one OpenBLAS thread. Official drivers use
two OpenMP threads per process and one OpenBLAS thread. Repeated native-unit
output on two MPI processes is not counted as independent physical sampling.
No performance or strong-scaling conclusion follows from these small fixtures.

The numerical module and native unit were additionally rebuilt with gfortran
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined
-fno-sanitize-recover=all`. The successful instrumented output is archived.
The surrounding prebuilt CP2K library was not fully rebuilt with these flags.
As in the production numerical wrapper, IEEE traps are disabled temporarily
inside external LAPACK calls. This also applies to the unit reference's LU
solve, whose OpenBLAS complex-division implementation otherwise triggers a
trap. No numerical reference, input Hamiltonian or convergence tolerance was
changed to obtain the passing result.

## Reproduce

Re-extract all quoted maxima and verify every archived checksum without CP2K:

```sh
python3 summarize_translation_quadratic.py translation-quadratic-records.tar.gz
```

The printed JSON should match `translation-quadratic-summary.json`.
The summary contains individual logs and distinguishes the MUMPS and dense
metric modes; its repeated records must not be summed into new sample counts.
Creating a fresh archive from the retained checkout uses:

```sh
python3 archive_translation_quadratic.py /path/to/cp2k new-records.tar.gz
```

The script refuses to overwrite an archive and requires the included source
files to match the committed tree. With the recorded source and suitable
serial/MPI CMake builds, run the native executables:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M build-serial/bin/quadratic_pseudospectrum_unittest.ssmp
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M build-serial/bin/quadratic_dbcsr_unittest.ssmp
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_STACKSIZE=64M mpiexec -n 4 build-mpi/bin/quadratic_dbcsr_unittest.psmp
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_STACKSIZE=64M mpiexec -n 4 build-mpi/bin/quadratic_dbcsr_unittest.psmp --column-grid
```

Repeat the latter two commands with one and two ranks. The MUMPS-enabled build
runs both metric modes; without MUMPS the dense-metric mode still runs.
The exact restricted regression directory lists and settings are retained
in both final driver logs. The existing finite and periodic Gaussian
references were not changed by this extension.

## Interpretation

For arbitrary complex covariant operators and queries, the implemented normal
products use `(A-zS)^dagger S^-1 (A-zS)`. Replacing the adjoint with the forward
operator or discarding the commutator of the Hermitian components changes the
problem. Basis covariance and nonnormal negative controls test these points.

Gaussian continuum unfolding is still separate work. A projected unitary
translation generally leaves a nonnegative complement
`S - T^dagger S^-1 T`. The projected quadratic norm alone omits leakage outside
the Gaussian span. No translated Gaussian integrals, boundary smoothing,
converged material dispersion or full continuum residual is validated here.

The conceptual unfolding reference is Bairnsfather, Kaufmann, Loring and
Cerjan, [arXiv:2605.05423](https://arxiv.org/abs/2605.05423). The nonorthogonal
normal-product implementation and tests were independently written.
