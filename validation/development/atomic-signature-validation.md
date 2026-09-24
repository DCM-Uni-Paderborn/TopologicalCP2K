# Atomic signature matrix validation

Verified 24 September 2026. This builds on `wyckoff-validation.md` and
`atomic-band-validation.md`; those earlier records describe their individual
construction stages, not the combined matrix layer documented here.

## Scope

`topology_atomic_signatures` builds an integer atomic reference matrix from all
generated site families and their onsite irreps/corepresentations. Each column
keeps its site and local-representation identity at every supplied k point.
Rows are unitary irreps, or full corepresentations where antiunitary symmetry is
present. Provenance, row dimensions and total band dimensions are retained.
Nonmaximal sites are included; equal sampled columns are not merged.

The six fractional reciprocal points are Gamma, (1/2,1/2,1/2), (1/2,0,0),
(1/3,1/3,0), (0,0,1/4), and (0.173,0.217,0.319). These are an explicit test
sample, not all reciprocal symmetry strata for every space group. Conventional
translation sectors are preserved, not silently unfolded into a primitive cell.

## Complete checks

Every complete serial, MPI and debug sweep reports:

- 530 Hall settings, each with ordinary/grey and scalar/spinful variants;
- 2,120 atomic signature matrices;
- 31,079 provenance-preserving columns;
- 707,038 integer matrix entries;
- 217,553 directed column/segment compatibility checks.

Every k block gives the same physical band rank for every column. Antiunitary
partner/doubling conditions are enforced. The seven connections per column
start at Gamma and include the zero segment, the five other test points and a
Gamma image at (0,0,1). Independent row permutations recover integer linear
combinations of the columns. The basis changes from unitary irreps to
corepresentations as the little group changes, without losing physical rank.

Analytic inversion tests recover all sixteen centered atomic signatures and
the generic-site signature on all eight TRIM, both scalar and spinful-grey.
Changing integer operation representatives and their order gives the expected
Bloch rephasing. Duplicate/missing operation rows, odd Kramers multiplicities
and empty reciprocal sets are rejected. Character decomposition assumes the
same supplied coset/spin gauge; it does not prove physical subspace closure.

## Evidence

- `atomic-signature-serial-all-settings.log.gz`: complete optimized serial sweep.
- `atomic-signature-mpi-all-settings.log.gz`: two ranks, each running the kernel.
- `atomic-signature-debug-all-settings.log.gz`: bounds checks and invalid/zero/overflow traps.
- `atomic-signature-serial-regtests.log`, `atomic-signature-mpi-regtests.log`:
  50/50 focused driver assertions in each build.
- `atomic-signature-without-database.log`: analytic tests without the test's
  SPGLIB preprocessor option, linked against the existing CP2K library.
- `signature_reallocation.f90`: isolated local compiler-behavior reproducer.

The regular CI sample uses four Hall settings, two reciprocal points and all
four ordinary/grey/spin variants. It includes centering, screws and complex C3
characters. The full matrix sweep is requested with `--all-settings` on the
`topology_atomic_unittest` executable. MPI execution here validates replicated
native algebra; it is not a distributed catalogue or multi-node scaling test.

## Local runtime and compiler qualification

Apple ARM runs use one OpenBLAS thread and the retained thread-safe OpenBLAS
installation. MPI uses two ranks and two OpenMP threads per rank. The focused
driver sets OMP_STACKSIZE=128M. Debug compilation of the new kernel and unit test
uses `-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.

With local GNU Fortran 16.1.0 at -O3, a minimal allocatable integer vector assigned
from MATMUL retained its old extent when the matrix slice changed from one row
to two. The reproducer contains no CP2K or BLAS calls. The test now allocates
the expected-vector extents explicitly at each k point and deallocates them
after use. The reference matrices and their physical values were not modified.
No upstream compiler bug report or fix is asserted by this observation.

## Remaining requirements

This is an atomic signature generator for explicit k sets, not an EBR
elementarity classifier or a new production input mode. A general workflow
still needs reciprocal-stratum coverage, integer-lattice and nonnegative
solvability, elementarity/equivalence analysis, and matching of the reference
convention to actual Gaussian-band data. Neither equality of sampled signatures
nor a nonnegative atomic signature decomposition proves a trivial Bloch bundle.
