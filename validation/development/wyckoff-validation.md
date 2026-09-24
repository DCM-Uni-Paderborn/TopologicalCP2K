# Native Wyckoff-family validation

Verified 24 September 2026. This is geometric/reference validation, not a
classification of 230 converged materials or a general EBR classifier.

## Implementation

The native `topology_wyckoff` module solves all periodic fixed equations
`(W-I)x=-tau mod Z^3`, forms all nonempty intersection components, and identifies
space-group-equivalent affine families. Primitive integer normal lattices
distinguish connected components. It uses no coordinate table or fixed sampling
denominator. Generic seeds are checked against complete stabilizers before
they are returned. Arithmetic/work limits return an error, never partial output.

The output consists of affine closures, generic seeds, free dimensions,
multiplicities, stabilizer indices and specialization indices. Maximal site
groups have no higher-symmetry specialization. They are not automatically
elementary band representations. Standard Wyckoff letters are not assigned.

## Complete sweep and independent reference

- All 530 Hall settings: 3,467 site-symmetry families.
- Independent spinful induction confirms each generated orbit and stabilizer.
- Every generated family matches exactly one tabulated affine-coordinate family
  with the same dimension and multiplicity, and conversely.
- Largest affine matching residual: 2.220e-16.
- All 5,648 directed specialization relations match the independent affine
  parameterization comparison, including non-immediate specializations.
- A combined arbitrary origin shift, nonorthogonal unimodular shear, integer
  operation-representative shifts and operation reordering preserves every
  family and the entire specialization graph in all 530 settings.

The independent table is
[SPGLIB Wyckoff.csv](https://github.com/spglib/spglib/blob/a6b561fb60cdd021a1ac90852c3c2ae14405b2c8/database/Wyckoff.csv),
SHA256 `d3d786a1f0187e5c6d69a3ade35648ffab34fd1b977d61ad84d8b0434b8b7ca0`.
The table's affine coordinates are parsed using a restricted expression tree,
not evaluated as Python source. External SPGLIB 2.7.0 supplies the operation
database. The CP2K generator reads operations but not the Wyckoff-coordinate
table. Both comparisons share the supplied group setting; this is not an
independent test of structure-to-space-group detection.

## Retained evidence

- `wyckoff-serial-all-settings.log.gz`: complete serial sweep.
- `wyckoff-mpi-all-settings.log.gz`: two ranks, each independently running the
  kernel. This is not a distributed enumeration or a multi-node scaling result.
- `wyckoff-debug-all-settings.log.gz`: bounds checks and invalid/zero/overflow traps.
- `wyckoff-reference.log`: external bijective-coordinate and specialization comparison.
- `check_wyckoff.py`: independent comparison implementation.
- `wyckoff-serial-regtests.log`, `wyckoff-mpi-regtests.log`: 49/49 assertions each.
- `wyckoff-without-database.log`: analytic tests compiled without `__SPGLIB`.

The regular driver combines the bounded Wyckoff unit test, the prior atomic-band
unit test and 47 Gaussian little-group assertions. The new bounded Wyckoff unit
costs 0.38 s serial and 0.51 s with two MPI ranks locally. These are smoke timings,
not scaling measurements. The full sweep is opt-in through `--all-settings`.

Analytic cases include P1, eight distinct inversion centers, two mirror planes,
a fixed-point-free screw, all 27 Pmmm families and their specializations, grey
groups, an antiunitary half translation, and rejected nonclosed operations.
The no-database executable omits the SPGLIB test interface but links the existing
CP2K library; it is not a complete CP2K rebuild without SPGLIB.

## Reproduction

Build CMake target `topology_wyckoff_unittest`, then run the versioned executable
with `--all-settings`. On the local Apple ARM builds, set one OpenBLAS thread and
use the retained thread-safe OpenBLAS library. The MPI sweep uses two ranks and
two OpenMP threads per rank. The focused driver additionally sets
`OMP_STACKSIZE=128M`. The standalone debug kernel and test use
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.

Download the linked table and check its SHA256, decompress the serial log, then
run `python check_wyckoff.py Wyckoff.csv wyckoff-serial-all-settings.log` in an
environment with NumPy and SPGLIB. The standalone script deliberately uses the
small standard tabulated-coordinate representatives for its independent
specialization-lift search, not arbitrary sheared coordinates.

## Remaining requirements

This closes the geometric representative-generation block, not general TQC.
Representation elementarity, complete atomic signature matrices, appropriate
primitive/conventional-cell matching, complete reciprocal connectivity and
comparison with actual Gaussian band data remain distinct requirements. Neither
matching symmetry signatures nor a nonnegative atomic decomposition alone proves
a trivial Bloch bundle. Magnetic electronic-structure detection is not supplied
by tests of antiunitary spatial groups.
