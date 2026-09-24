# Reciprocal fixed families and sampling graph

Verified 24 September 2026. The CP2K implementation is committed locally;
this record does not imply publication of that research branch.

## Implemented construction

`topology_reciprocal` applies the already validated affine fixed-torus engine
to the effective reciprocal actions K=(-1)^a W^(-T). It deduplicates geometric
actions, not the physical group: full centering, nonsymmorphic and antiunitary
cosets remain in every returned original-operation stabilizer. Original group
closure is validated at Gamma before the geometric quotient is formed.

For each connected closure N*k=c modulo integers, an exact Smith certificate
provides a primitive Z-basis for ker(N). For every higher-symmetry specialization,
the constructor finds an operation and integer shift embedding both the child
seed and its tangent lattice in the parent's lifted closure. The integer lift
is solved with the checked lattice solver. Resource/arithmetic failures return
no partial result. The code bounds integer coefficients before products.

`reciprocal_sampling` wraps this into explicit points and directed edges.
The first points are generic stratum seeds; further symmetry-related child
images are evaluated explicitly. Integer destination shifts are retained,
including self-edges for primitive torus cycles. These outputs feed the existing
atomic reference and segment-compatibility APIs without an assumed irrep-label
transport between star arms.

This follows the reciprocal-space-group/Wyckoff correspondence described by
the [Bilbao definitions](https://cryst.ehu.es/cryst/help/definitions_kvec.html).
The implementation uses native geometry, not Bilbao data at runtime.

## Complete Hall-setting evidence

For each of 530 conventional Hall settings, ordinary and grey groups are tested:

- 1,060 geometric groups;
- 18,700 connected fixed-family closures;
- 57,357 specialization incidences;
- 18,730 primitive integer tangent cycles;
- 333,900 rational-grid point classifications;
- 2,481,572 scalar/spinful atomic-column directed-edge comparisons.

The rational meshes have denominators 2, 3, 4 and 6. A separate direct-transpose
test determines each point's stabilizer; its orbit must match exactly one
generated family of the same order. This is independent finite-grid evidence,
not an exhaustive numerical sampling of the continuum.

Full spinful little groups independently verify all generic stabilizers.
Integer cycle lattices are checked for rank and primitivity, incidence lifts
must preserve the parent's stabilizer pointwise, and simultaneous unimodular
shear/origin/coset/operation-order changes must preserve the incidence graph.
Analytic examples cover P1, time reversal, PT, magnetic half translations,
nonsymmorphic screw lines, centering, and orthorhombic intersections. Invalid
groups return no partial family list.

The atomic reference matrices are generated for both scalar and spinful factors
at the graph's seeds/images. Every column must satisfy every generated
unitary restriction and, where present, every antiunitary corepresentation
restriction. Cycles preserve reciprocal shifts, so nonsymmorphic monodromy is
not removed by folding endpoints onto the same point. This tests constructed
atomic reference bands, not new Gaussian material calculations.

Serial, two-rank MPI and instrumented runs pass the complete sweep. MPI repeats
the algebra on each rank, not a distributed construction or scaling benchmark.
The instrumented run compiles the new module and its unit program with
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined
-fno-sanitize-recover=all`; other modules come from the existing library.
An additional test with the unit program's SPGLIB macro disabled passes the
analytic subset, linked against the existing library. This is not a full
no-SPGLIB rebuild.

The official focused drivers pass 74/74 serial and 74/74 two-rank MPI checks,
consisting of five native programs and 69 unchanged Gaussian assertions.
`make_pretty.sh` passes and is byte-identical on the compiled source. No physical
reference values or output suppression thresholds were changed.

## Records and reproduction

- `reciprocal-strata-serial-all-settings.log.gz`
- `reciprocal-strata-mpi-all-settings.log.gz`
- `reciprocal-strata-debug-all-settings.log.gz`
- `reciprocal-strata-without-database.log.gz`
- `reciprocal-strata-serial-regtests.log.gz`
- `reciprocal-strata-mpi-regtests.log.gz`
- `reciprocal-strata-pretty.log.gz`

Run `topology_reciprocal_unittest --all-settings` for the complete sweep.
The default CI selection has nine representative Hall settings plus the
database-free analytic tests. The focused driver restricts directories to
`UNIT/topology_(integer|atomic|wyckoff|band|reciprocal)_unittest|QS/regtest-little-group`.
Local runs use GNU Fortran 16.1.0 on Apple ARM, one BLAS thread, the retained
thread-safe OpenBLAS preload, and OMP_STACKSIZE=128M for Gaussian runs.
Serial regtests use one rank/two threads, MPI two ranks/two threads. Recorded
elapsed times are test costs with concurrent local jobs, not scaling evidence.

## Remaining mathematical and physical work

This completes the fixed-closure/incidence/cycle geometric backend, not a full
general TQC classifier. Open strata may disconnect when higher-symmetry sets
are removed; their equivariant cell complex and full global band connectivity
still require construction. Primitive-cell equivalence, EBR elementarity and
standard labels remain separate. The actual `ATOMIC_SIGNATURES` Gaussian path
still uses explicit NNKP/Wilson/TRIM samples; automatic sampling and spectral
isolation checks along the resulting geometry are not yet connected to it.
