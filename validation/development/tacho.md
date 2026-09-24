# Tacho Sparse Pfaffians

`CP2K_USE_TACHO=ON` adds optional native sparse real skew-LDL Pfaffian signs for
`FORCE_EVAL/PROPERTIES/SPECTRAL_LOCALIZER/SOLVER TACHO` with `INVARIANT Z2`. It also requires
`CP2K_USE_MUMPS=ON`: MUMPS independently checks the positive AO metric and brackets the generalized
localizer gap. Neither a determinant nor ordinary inertia can replace the orientation-sensitive
Pfaffian sign.

The required library is the **skew-LDL extension**, not the standard Trilinos release:
[iyamazaki/Trilinos, branch tacho-sk](https://github.com/iyamazaki/Trilinos/tree/tacho-sk). The
tested source revision is `02ef047cb4e1e16959714ce50462bb09511d0e9a`. Tacho is BSD-2-Clause
licensed. The native adapter uses Boost Graph's header-only general-graph matching (Boost Software
License) and the serial Kokkos execution space. No MATLAB, SuperLU_DIST or proprietary matching
library is required. No external sources are copied into CP2K.

Build Trilinos with `Trilinos_ENABLE_ShyLU_NodeTacho=ON`, `Kokkos_ENABLE_SERIAL=ON`,
`Tacho_ENABLE_INT_INT=ON`, BLAS/LAPACK enabled, and MATLAB/tests/examples disabled. Set
`ShyLU_NodeTacho_DIR` to the installed `lib/cmake/ShyLU_NodeTacho` and provide Boost's CMake
configuration through `CMAKE_PREFIX_PATH` as needed. Use the same C++ ABI and BLAS integer width as
CP2K. Tacho's skew API is experimental; pin the tested revision when reproducing results.

## Algorithm and Safeguards

DBCSR assembly and MUMPS gap checks remain distributed. Only sparse coordinate data are gathered to
the source rank for the native Pfaffian factorization. The fixed AII transformation
`Q = (I - i C)/sqrt(2)` has two entries per row and is applied by sparse permutations/additions,
without dense localizer or metric matrices. The imaginary residual and real skew structure are
checked; the discarded part must also be small relative to the independently bounded gap in the AO
metric.

Connected components are factored separately. A largest-entry-first greedy matching is augmented
with general-graph Edmonds matching to supply nonzero 2x2 pivots; its permutation parity is
retained. Matched pairs stay together during fill-reducing ordering. Numerical pivoting is enabled
and pivot perturbations are disabled. The adapter always bypasses Tacho's incomplete dense
small-problem branch, and evaluates isolated 2x2 components analytically. Tacho multiplies lower
skew-pivot entries; one minus sign per pivot converts this to the conventional upper-entry Pfaffian.
Both matching and ordering parity are included.

The atomic reference orientation is `(-1)**nao` in CP2K's fixed auxiliary/spin/AO order. A positive,
time-reversal-compatible metric can be continuously deformed to the identity without closing this
reference gap, so its sign needs no additional factorization. The physical localizer is not
orthogonalized or projected onto occupied states.

Four deterministic sparse-solve probes diagnose inaccurate factors. Their residual is reported, but
is not a rigorous backward-error certificate for the entire factorization. Small pivots, nonfinite
factors, failed probes or an unresolved spectral gap do not yield an index. Solver exceptions are
propagated collectively. There is no silent dense fallback.

## Limits

This is a **single-rank sparse Pfaffian factorization**, not an MPI-distributed Pfaffian solver.
Factor fill-in, Boost matching and the source-rank coordinate storage limit the largest usable
problem. The reported MUMPS memory excludes Tacho, matching and coordinate buffers. The Pfaffian
symbolic upper 2x2-block count is a separate diagnostic, not scalar factor storage or peak memory.
The inspected Tacho factor-export function is unimplemented and is not used. Periodic systems and
k-points remain outside the finite-system localizer implementation.

The mathematical and independent-model validation is in `spectral_localizer_sparse_unittest`; the
optional GTH-SOC end-to-end inputs are in `QS/regtest-spectral-localizer-tacho`.

Local validation covers 338 dense/sparse Pfaffian comparisons over sizes 2 through 24, odd
permutations, disconnected components, numerical singularities with perfect graph matchings, and
scales from `1e-120` to `1e120`. Additional checks cover invalid coordinates, non-skew input,
additive duplicates and unresolved pivots. Spin-mixed QWZ/Kramers-pair models exercise both Z2
indices, complex nonorthogonal metrics, roundoff and broken time reversal, and gap closings. The
largest checked complex localizer has order 800 (realified order 1600).

The dense/MUMPS/Tacho molecular directories pass 23 checks with four MPI ranks and two OpenMP
threads; the dense directory also passes 12 checks in a build without the optional libraries. Neon
and bismuth GPW SOC gaps agree with the dense oracle. A separate bismuth GAPW SOC smoke test also
agrees, and the numerical Fortran modules pass bounds checks and floating-point traps with two MPI
ranks. These are solver and integration tests, not size-converged topological-material benchmarks or
a claim of multi-node Pfaffian scaling.
