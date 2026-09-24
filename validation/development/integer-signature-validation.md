# Exact integer signature validation

Verified 24 September 2026. This extends the atomic matrices described in
`atomic-signature-validation.md`. It does not add a material classifier.

## Implementation

`topology_integer_lattice` constructs U*A*V=D using checked 64-bit Euclidean
row and column operations. D has positive nonzero entries in divisibility
order. The full certificate is multiplied out before acceptance. The signed
solver checks transformed-target divisibility and all zero rows, then returns
an exactly checked A*x=b witness. The factorization is opaque and reusable.

The separate nonnegative solver exhaustively searches componentwise bounded
coefficients with suffix-gcd and zero-support pruning. It verifies each
witness and has a mandatory node budget. Only status zero is a decision.
An exhausted budget, coefficient growth, invalid data or failed verification
must not be interpreted as a symmetry obstruction. The algorithm can be
exponential and is not intended as a general large-scale integer optimizer.

## Independent checks

- 384 small rectangular/singular matrices, sizes one through four:
  gcds of all minors agree with products of Smith factors; both transformation
  determinants have magnitude one. Rank counts 0..4: 38, 167, 103, 60, 16.
- All 81 four-TRIM and 6,561 eight-TRIM parity patterns for two Kramers pairs
  agree with the independent analytic inversion Fourier criterion.
  Signed/nonnegative counts are 41/33 and 241/129, respectively.
- Every 2x3 matrix over {0,1,2}, with every target over {0,...,4}^2:
  18,225 cases agree with independently enumerated coefficients 0..4.
  A coefficient above four cannot contribute through a nonzero nonnegative
  column to one of these targets; zero columns can have coefficient zero.
- Invalid dimensions, negative data for the nonnegative solver, zero columns,
  deficient rank, overflow-sized input, intermediate coefficient growth and
  zero/exhausted search budgets exercise separate outcomes.

## Generated atomic matrices

All 530 Hall settings, scalar/spinful and ordinary/grey variants, use the six
explicit reciprocal points retained by the preceding validation stage.
Every sweep factors 2,120 matrices, recovers 2,120 signed test targets and
33,199 nonnegative witnesses (each atomic column and one three-column sum
per matrix). The constructed targets have known witnesses; this is not a
classification of occupied DFT bands or independent material calculations.
The pre-existing 217,553 column/segment comparisons remain enabled.

Both serial and MPI complete sweeps pass. The MPI test repeats native algebra
on two ranks: it is not a distributed solver or a multi-node performance test.
A debug build of the new integer kernel and matrix driver also passes the
complete sweep with `-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
-fsanitize=undefined -fno-sanitize-recover=all`. Other CP2K modules are linked
from the existing library rather than rebuilt with these flags.

## Focused integration and runtime

The official driver passes 51/51 assertions in both configurations:
four native unit programs plus QS/regtest-little-group (47 assertions).
MPI uses two ranks, two OpenMP threads and OMP_STACKSIZE=128M. Both runs use
one OpenBLAS thread and the local thread-safe OpenBLAS installation.
Compiler: local GNU Fortran 16.1.0, Apple ARM.

The new integer unit takes 0.19 s in the serial driver and 0.50 s under MPI.
The matrix unit takes 0.17 s and 0.47 s for its regular bounded CI sample.
These are local observed test costs, not performance guarantees.

The analytic matrix driver also passes without its SPGLIB preprocessor option.
It remains linked to the existing CP2K library, so this is not evidence of a
complete no-SPGLIB build. The integer kernel itself has no SPGLIB calls.
`make_pretty.sh` and `git diff --check` pass. No physical reference values changed.

## Retained logs

- `integer-atomic-serial-all-settings.log.gz`
- `integer-atomic-mpi-all-settings.log.gz`
- `integer-atomic-debug-all-settings.log.gz`
- `integer-signature-serial-regtests.log.gz`
- `integer-signature-mpi-regtests.log.gz`
- `integer-signature-debug.log`
- `integer-atomic-no-database.log`
- `integer-signature-pretty.log`

The first release sweeps precede the added printed integer-check totals; their
source already contains all integer checks. The debug sweep prints the totals,
and the focused driver runs the final formatted source.

## Interpretation

Membership refers only to the supplied integer matrix. Smith factors of an
incomplete reciprocal sample are not a complete symmetry-indicator group.
Nonnegative signatures do not prove triviality of a Bloch bundle; signed-only
signatures do not alone prove fragile topology. Complete reciprocal strata,
compatibility, elementarity/equivalence, and physical Gaussian-band mapping
remain distinct requirements. There is no new automatic general-TQC keyword.
