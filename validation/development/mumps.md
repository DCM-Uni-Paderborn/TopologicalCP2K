# MUMPS

The optional MPI MUMPS backend supplies sparse symmetric-indefinite inertia for finite-system
spectral localizers. Enable it with `-DCP2K_USE_MUMPS=ON` and provide `MUMPS_ROOT` (headers and
libraries), or `MUMPS_INCLUDE_DIR` and `MUMPS_LIBRARIES`. Static builds must include their ordering
dependencies in `MUMPS_LIBRARIES`, for example `dmumps;mumps_common;pord`. MUMPS must use the same
MPI implementation, Fortran ABI, and integer/BLAS convention as CP2K. There is no automatic download
and no new dependency for default CP2K builds.

`SPECTRAL_LOCALIZER/SOLVER MUMPS` constructs the operator in DBCSR. Complex Hermitian matrices are
realified as `[Re(A),-Im(A);Im(A),Re(A)]`; each eigenvalue is duplicated. Distributed lower-triangle
coordinates are passed to MUMPS without gathering a dense matrix. Separate symbolic analyses are
reused for the metric stage and for all indefinite gap-bracketing shifts at one query point. The
metric uses `SYM=1` (internally non-pivoted LDL, with an explicit negative-pivot check); the
localizer uses `SYM=2` with dynamic pivoting. Neither stage uses static pivot perturbation or
low-rank compression. `ICNTL(13)=1` is essential: otherwise the ScaLAPACK root node does not
contribute to the reported negative inertia. The remaining factorization remains distributed. See
the [MUMPS FAQ](https://mumps-solver.org/index.php?page=faq).

When numerical fill exceeds the symbolic workspace estimate (MUMPS errors -8, -9, -17 or -20), the
adapter retries the same numerical factorization with increased `ICNTL(14)`, at most five times per
call, with the relaxation capped at 640 percent. It does not change matrix values or pivot
tolerances. Other failures are not retried, and native error codes are reported on failure.
Workspace retries count toward the reported factorization total.

The optional solver currently implements the class A half-signature, including complex SOC
Hamiltonians. It is **not** a Pfaffian-sign solver. Optional [Tacho](tacho.md) supplies class AII/Z2
Pfaffian signs while reusing MUMPS for the metric and gap checks. The returned protection gap is a
numerical lower/upper bracket for the generalized pencil, not a pivot magnitude. The metric ratio is
a conservative lower bound using a row-norm upper bound for the largest overlap eigenvalue.
Ill-conditioned or unresolved cases are rejected conservatively.

MUMPS is distributed under CeCILL-C; see the
[upstream licensing and references](https://mumps-solver.org/index.php?page=dwnld). This integration
does not vendor MUMPS sources. Sparse fill-in and the root-front size can still limit memory and
scaling, and must be measured for the intended systems.

## Reusable metric solutions

The [quadratic pseudospectrum](../methods/properties/quadratic_pseudospectrum.md) uses the same
adapter for positive AO-overlap factors. `sparse_metric_start` checks the metric and retains its
factorization; `sparse_metric_solve` executes MUMPS `JOB=3` for multiple right hand sides without
refactorization. Realification also supports complex AO coefficients. The factors are reused for all
matrix-free quadratic-operator applications and query points. Right hand sides are centralized on
communicator rank zero and solutions broadcast to the bounded replicated Ritz space. This avoids a
dense overlap inverse but does not remove sparse fill-in or vector-communication costs. The MPI unit
test includes repeated multi-RHS solves and verifies that the factorization counter does not
increase during solutions.

## Validation and indicative timings

Local tests on 2026-09-24 used an Apple M4 Max (14 cores, 36 GiB), GNU Fortran 16.1, OpenMPI 5.0.9,
OpenBLAS 0.3.33 and MUMPS 5.9.1 with PORD ordering. Both OpenMP and BLAS used one thread in the
following kernel measurements. Each row is a single final run, not a statistical performance
guarantee. Reproduce it with:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 \
  build-mpi/bin/spectral_localizer_sparse_unittest.psmp 30
```

The argument is the side length of an open Qi-Wu-Zhang model. The complex localizer has order
`4*length**2`; realification doubles it. Dense times include source-rank generalized diagonalization
and the LDL cross-check. Sparse times include symbolic analysis, metric checks and all shifted
factorizations needed for the gap, but not test-fixture/COO construction.

| Mesh    | Realified order | MPI ranks | Dense [s] | Sparse [s] | Factor/input entries | Solver memory sum [MB] |
| ------- | --------------: | --------: | --------: | ---------: | -------------------: | ---------------------: |
| 10 x 10 |             800 |         1 |    0.0332 |     0.0585 |                 15.0 |                      1 |
| 10 x 10 |             800 |         2 |    0.0333 |     0.0592 |                 15.2 |                      2 |
| 10 x 10 |             800 |         4 |    0.0343 |     0.0568 |                 15.1 |                      4 |
| 20 x 20 |            3200 |         1 |     1.896 |      0.513 |                 21.9 |                      6 |
| 20 x 20 |            3200 |         2 |     1.889 |      0.488 |                 22.3 |                      9 |
| 20 x 20 |            3200 |         4 |     1.899 |      0.337 |                 22.5 |                     17 |
| 30 x 30 |            7200 |         1 |    22.945 |      2.401 |                 36.0 |                     19 |
| 30 x 30 |            7200 |         2 |    23.120 |      2.207 |                 37.2 |                     32 |
| 30 x 30 |            7200 |         4 |    23.235 |      1.370 |                 36.4 |                     45 |

All configurations give index -1, with the dense generalized gap inside the sparse bracket within
the unit-test tolerance. The test also repeats the comparison after a nonorthogonal basis change.
The larger four-rank case exercises workspace recovery. Stored factor/input ratios report the
maximum factor entries across the shifted factorizations divided by the input pattern size; they
illustrate fill-in, not linear memory scaling.

A single dense complex matrix of order 3600 occupies 207.4 MB before eigensolver workspace. The
MUMPS memory column counts only solver allocations, excluding application inputs and DBCSR. The test
deliberately retains dense reference matrices, so its process RSS must not be presented as the
memory consumption of the production sparse path. Dense is faster for the smallest case; MPI speedup
is limited and trades time for additional memory. These are model-kernel timings, not DFT-material
benchmarks or a comparison against parallel dense eigensolvers.

Additional local checks cover zero-diagonal 2-by-2 pivots, near/exact gap closings, complex and
invalid metrics, and a kernel build with bounds checks and floating-point traps. The molecular
directories pass 12 dense checks (serial and MPI) and 7 optional sparse checks (MPI), including GPW,
GAPW and GTH SOC. Missing-MUMPS and Z2 requests with `SOLVER MUMPS` are explicitly rejected. The
entire upstream CP2K test suite and multi-node scaling have not been validated by these checks.

## Backend selection

The integration was tested with MPI MUMPS 5.9.1. PEXSI is already optional in CP2K, but its current
CP2K adapter is a real DFT-density driver, not a general inertia service. The upstream
[expert interface](https://pexsi.readthedocs.io/en/latest/tutorial.html) is a possible future
backend; its distributed CSC conversion, ordering and supported matrix types need independent
validation. It was not available in the local test installation.

[oneMKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/pardiso-iparm-parameter.html)
also reports inertia for symmetric-indefinite systems. It is an attractive option where MKL is
already linked, but perturbed pivots must not silently be accepted as an exact topological index.
oneMKL PARDISO and the separately distributed PANUA PARDISO are different products; neither is
linked by this integration. The Apple ARM/OpenBLAS test configuration cannot validate an MKL
backend.

Neither MUMPS nor ordinary PARDISO inertia supplies the orientation-sensitive Pfaffian sign needed
for class AII. The
[Tacho skew-LDL extension](https://github.com/iyamazaki/Trilinos/tree/tacho-sk/packages/shylu/shylu_node/tacho)
is integrated separately through optional `CP2K_USE_TACHO`. Revision
`02ef047cb4e1e16959714ce50462bb09511d0e9a` has a BSD-2-Clause Tacho package and a native C++
`SkewLDL`/`pfaffian()` interface. A small C ABI shim and Fortran `ISO_C_BINDING` avoid MATLAB
entirely. See [the adapter documentation](tacho.md) for Kokkos/Boost dependencies, permutation
conventions, error propagation and scaling limits.

There is an important readiness issue: at that revision, `Driver::pfaffian()` returns zero on its
small-problem branch instead of computing a sign. The adapter bypasses this branch and tests
singular matrices, permutation parity, near-zero pivots and the same orientation calibration as the
dense oracle. Tacho is selected explicitly, never silently. The
[MATLAB research wrapper](https://github.com/acerjan/sparse_sign_pfaffian) is not a CP2K dependency.
