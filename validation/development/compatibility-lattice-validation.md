# Joint integer compatibility lattice

Verified 24 September 2026. This is a local CP2K research implementation,
not a new published CP2K pull request or a complete space-group classifier.

## Construction and scope

All supplied directed edge restrictions are assembled in the same endpoint
unitary-irrep/corepresentation coordinates as the generated atomic matrix A.
The equations C*n=0 include endpoint corepresentation restrictions and reciprocal
target images. A self-edge adds both endpoint terms to the same block, so a
screw/glide sector permutation is not discarded as an identity connection.
Homogeneous rows are normalized by gcd/sign and deduplicated exactly.

The checked integer layer computes a saturated kernel K of C, solves A=K*Y,
and factors Y in Smith form. It verifies C*K=0 and K*Y=A. Positive divisors
give cyclic factors (1 is trivial); zeros denote free integer coordinates.
For b=K*z, the Smith row map U_Y gives its quotient class. A target violating
compatibility has no class. Zero class is equivalent to signed atomic membership,
not nonnegative membership or Bloch-bundle triviality.

The additional analysis uses existing ATOMIC_SIGNATURES with compatibility
enabled, including the automatic SYMMETRY graph. It needs no new keyword or
runtime dependency. Matrix storage estimates include edge and integer-certificate
work but exclude external-library workspace. Budget/arithmetic failures are
unresolved, not topological obstructions. The integer construction is replicated
on MPI ranks, not a distributed factorization or scaling benchmark.

Missing star-arm identifications, open-cell relations or global constraints
can leave free factors. Even finite quotients do not prove completeness of
omitted congruence relations, EBRs or reciprocal coverage. Sampled isolation
does not prove a bulk gap. These limitations remain explicit in output and text.

## Independent and physical checks

- Inversion parity columns at four/eight TRIM with equal-rank constraints give
  Z2 and Z2^3 x Z4. All 81/6,561 two-pair patterns have zero class exactly when
  the separately tested signed atomic membership holds.
- Nonsaturated constraint rows, a Z plus Z2 quotient, zero-dimensional kernels,
  incompatible targets/generators, reversed edges and exhausted budgets are tested.
- The complete 530 Hall-setting ordinary/grey and scalar/spinful graph sweep
  checks all generated atomic columns in the joint kernel and in zero quotient
  class, in addition to 2,481,572 individual edge/column comparisons.
- Serial and two-rank MPI official drivers each pass 165/165 assertions:
  129 in the native/little-group group and 36 existing exporter/k-point checks.
  The 23 added assertions do not change existing physical reference values.
- The new integer/compatibility modules, reciprocal module and its unit program
  are additionally compiled with bounds, floating-point traps and undefined-behavior
  instrumentation. Other dependencies use the existing library.

The full serial and instrumented sweep logs are byte-identical. After removing
MPI rank prefixes, both rank-specific geometry/quotient summaries agree with
serial output. Each sweep has 2,120 database-backed quotients and four additional
analytic graph cases. The instrumented database-free subset also passes.

SymPy 1.14.0 independently verifies eleven Gaussian quotient exports per serial/MPI
run. It computes rational rank(C), Smith(K) (all kernel invariant factors one),
and Smith(Y), and checks C*K=0, K*Y=A and C*A=0. The quotient row transform must
be unimodular, have zero free rows on Y, and have primitive normalized nonzero
rows after division by its divisors. Target coordinates and residues are
recomputed with arbitrary-size integers. Serial/MPI summaries agree exactly.

| Fixture | Equations | Coordinates | Kernel rank | Free rank | Nontrivial finite factors | Target |
|---|---:|---:|---:|---:|---|---|
| Automatic He/GPW, He/GAPW, sheared He | 109 | 91 | 6 | 0 | 2 | compatible, class zero |
| Automatic Ne/SOC | 81 | 67 | 5 | 0 | 2, 2, 4 | compatible, class zero |
| Explicit He screw / GAPW screw | 8 | 9 | 3 | 0 | none | compatible, class zero |
| Explicit Ne/SOC screw | 4 | 5 | 2 | 0 | none | compatible, class zero |
| Partial He / sheared He graph | 1 | 8 | 8 | 4 | none | compatible, class zero |
| Refined He Wilson graph | 44 | 75 | 42 | 23 | none | compatible, class zero |
| Single He screw branch | 1 | 2 | 1 | 0 | none | incompatible, no class |

The automatic SOC example's Z2^2 x Z4 must not be confused with the analytic
eight-TRIM inversion example's Z2^3 x Z4. The physical cases are small atomic
fixtures; these tests do not establish a new nontrivial material phase.

## Reproduction and records

Build the usual cp2k-bin, topology_integer_unittest and topology_reciprocal_unittest
CMake targets. Use topology_reciprocal_unittest --all-settings for the complete
native sweep. The regression-directory restriction is
`UNIT/topology_(integer|atomic|wyckoff|band|reciprocal)_unittest|QS/regtest-little-group|QS/regtest-topology$|QS/regtest-kp-1$`.
SSMP uses one rank/two OpenMP threads, PSMP two ranks/two threads. Local Gaussian
tests set one BLAS thread, OMP_STACKSIZE=128M and use the retained thread-safe
OpenBLAS preload. GNU Fortran 16.1.0 and SPGLIB 2.7.0 are used on Apple ARM.
Elapsed times under concurrent local jobs are not performance/scaling results.

The instrumented compile uses
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined -fno-sanitize-recover=all`.
The database-free unit variant disables the unit program's SPGLIB macro and
links against the existing library; it is not a full no-SPGLIB build.

Retained files with prefix `compatibility-lattice-` include compressed serial/MPI
regtest and all-setting logs, instrumented tests, formatting checks and Gaussian
certificate archives. `verify_compatibility_quotients.py` requires SymPy 1.14.0
only for independent validation. After extracting either certificate archive:

```sh
python3 verify_compatibility_quotients.py *.little_group
```

Each record contains the physical target, atomic columns and character/provenance
rows as well as the new sparse certificate entries. No CP2K source, executable,
SCF restart or large scratch build is copied into this manuscript repository.
