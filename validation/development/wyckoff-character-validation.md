# Actual Wyckoff onsite characters and independent Bloch induction

Validated on 25 September 2026. CP2K local commit
`e32a2c63926b4bfd76a417d33cd90caef96eb3c3` adds an optional development
test path and its documentation. It does not change production numerical
methods or Gaussian reference values, add a runtime database, or introduce
a production input keyword. The CP2K commit has not been pushed upstream.

## Coverage and results

- All 1,731 independently archived conventional Wyckoff families of the
  230 ordinary space groups, including nonmaximal sites.
- Scalar/spinful and ordinary/grey variants: 6,924 cases and 16,957 onsite
  irrep/corepresentation columns before coordinate repetition.
- Two coordinate conventions: 13,848 complete onsite-table matches and
  189,664 onsite character entries per executable.
- Six primitive reciprocal points per case/convention: 83,088 complete
  induced character tables and 1,384,328 induced character entries.
- Largest onsite external residual: 3.1415927e-5 (combined spin/frame/character).
- Largest onsite cross-build difference: 8.21568855656745e-15.
- Largest onsite cross-coordinate difference: 6.217283432986974e-15.
- Largest induced reference difference: 1.4510351e-4 absolute;
  largest dimension-scaled difference: 1.0471976e-5.
- Largest induced cross-build difference: 1.865174681370263e-14.
- Largest induced cross-coordinate difference: 2.6645352591003757e-14.

The native geometry and algebra tolerance stays 1e-9. External printed
reference characters/spin matrices use 2e-4. For induced characters this
bound is scaled by the positive whole-band dimension of each column,
not by the largest column in the table. Absolute and scaled errors are
both retained. Cross-build/coordinate character checks use an absolute
1e-9 bound. No tolerances or physical reference values were relaxed.

The optimized serial run, separately each of two MPI ranks, and the focused
instrumented executable pass the complete sweep. MPI repeats the audit on
each rank; the counts are not distributed speedup evidence. The default
native program also passes in SSMP and two-rank MPI with two OpenMP threads.
The preceding 179/179 serial and MPI driver results are retained in the
onsite-reference report; they were not rerun or counted as new Gaussian
calculations for this test-only extension.

## Independent reference construction

`prepare_wyckoff_characters.py` reads the exact fixtures retained with the
affine/character-reference reports. The affine data are from the pinned
Crystalline.jl archive; the character data are from irreptables 3.1.0. The
generator imports neither CP2K nor their representation-generation code.
It uses NumPy, SciPy and exact integer Hermite arithmetic from SymPy.
The recorded environment is Python 3.14.6, NumPy 2.5.3, SciPy 1.18.1 and
SymPy 1.14.0 on arm64 macOS 15.7.4; the archive retains the environment record.
The full manifest records SHA-256 hashes of both inputs and the resulting
fixture, plus every site's expected labels, orbit size, stabilizer order,
column count and reciprocal-operation membership.

Each family is sampled at its supplied affine origin plus tangent parameters
sqrt(2)/7, sqrt(3)/7 and sqrt(5)/7. These generic positions need not coincide
with the native generator's preferred seeds. Conventional centering is removed
using an explicit primitive basis computed by integer Hermite normal form.
The original fractional translations are retained modulo this primitive
lattice. Expected primitive multiplicities are the tabulated conventional
multiplicities divided by the centering factor.

Complete scalar and spinful Gamma character tables of 32 symmorphic point-group
representatives supply the onsite references. The correspondence between each
site stabilizer and its reference point group is established by the full
Cartesian operation set, using an explicitly recorded proper orthogonal
conjugacy. Grey restrictions are formed independently through conjugate
character pairing and the spin-aware Frobenius--Schur indicator. These are
declared reference axes and labels, not an assertion that the first allowed
conjugacy is Bilbao's standard site-label convention.

For the Bloch oracle, the Python generator explicitly constructs the orbit
and unitary spatial orbit representatives t_j. This suffices for ordinary
groups and their grey extensions, since pure time reversal fixes positions.
For every g fixing site j modulo an integer L_j, it evaluates h_j=t_j^-1 g t_j
in the onsite point group and sums

```text
chi_ind(g,k) = sum_j exp(-2*pi*i*k.L_j) alpha(g,t_j,h_j) chi_site(h_j).
U(t_j)^dagger U(g) U(t_j) = alpha V U_reference(h_j) V^dagger.
```

Physical spin lifts of the full spatial operations come from SciPy's
Cartesian rotation quaternions; the local reference spin matrices remain
the published finite-precision values. The overlap coefficient is checked
to be unit-modulus within the reference precision before normalizing its
phase. The onsite characters themselves are not rounded or fitted.
The native output is normalized only by explicit spin-lift and integer-coset
conventions, independently of the irrep column. The six k points are
(0,0,0), (1/2,0,0), (0,1/2,0), (0,0,1/2), (1/4,1/4,1/4),
and (0.173,0.217,0.319) in the recorded primitive basis.

The second coordinate convention simultaneously applies a unimodular shear,
a 0.73-radian Cartesian rotation, origin shift (0.137,0.219,0.317), reversed
operation order and integer coset shifts. Every reference-normalized onsite
and induced character, not merely a dimension or irrep count, is compared.

## Checks that reject incorrect data

The independent reader checks literal labels, complete case membership,
all operation/label pairs and both MPI ranks. It detected a Fortran
list-directed I/O problem: labels such as `2*GM1` were interpreted as
repetition syntax. The development reader now uses `(A)` for literal
labels. The corrected run retains the full expected names without changing
any numerical reference character.

`corrupt_wyckoff_phase.py` selects a scalar onsite character that is real
on every onsite operation but has a genuinely complex induced Bloch trace.
For space group 100, site c, label GM1 at (1/4,1/4,1/4), conjugating one
induced entry changes -2i to +2i. The native reference test rejects this
with `Independent Wyckoff induced character mismatch`, absolute error 4,
dimension-scaled error 1 and nonzero exit status. This tests a translation
phase, not a change to a complex onsite irrep or to the dimension.

## Retained evidence and replay

`wyckoff-character-records.tar.gz` contains the full independent fixture and
manifest, all serial/MPI/instrumented character output, normal unit logs,
negative fixture/log, source files and build/formatting logs. Archive member
hashes are listed in `wyckoff-character-files.sha256`. The top-level JSON
summary is generated from the archive, not from expected numbers alone.
The replay script supports `archive.tar.gz::member` paths so approximately
800 MiB of raw logs need not be extracted just to verify the results.

```sh
python3 verify_wyckoff_characters.py \
  wyckoff-character-records.tar.gz::build-serial/wyckoff-bloch-full.json \
  wyckoff-character-records.tar.gz::build-serial/wyckoff-bloch-full.log \
  wyckoff-character-records.tar.gz::build-mpi/wyckoff-bloch-full.log \
  wyckoff-character-records.tar.gz::build-serial/wyckoff-bloch-full-debug.log
```

To regenerate the fixture from the retained input data:

```sh
tar -xzf wyckoff-character-records.tar.gz \
  build-serial/character-convention-reference.dat build-serial/affine-reference.dat
OPENBLAS_NUM_THREADS=1 python3 prepare_wyckoff_characters.py \
  build-serial/character-convention-reference.dat build-serial/affine-reference.dat \
  wyckoff-bloch-full.dat wyckoff-bloch-full.json
```

Run the recorded CP2K source's `topology_band_unittest.ssmp` or
`topology_band_unittest.psmp` with `--wyckoff-reference=wyckoff-bloch-full.dat`,
OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The MPI audit uses
`mpiexec --tag-output -n 2`. The existing `build_onsite_instrumented.sh`
recompiles the symmetry/representation modules and test with bounds checks,
floating-point traps and undefined-behavior instrumentation. It links the
existing release CP2K library and omits its SPGLIB database test interface;
it is not a fully instrumented Gaussian/MPI or independent no-SPGLIB build.

## Interpretation boundaries

This tests generic representatives of every listed actual Wyckoff family
and their independently induced band characters. It does not sample every
position in a free family, prove full reciprocal-stratum connectivity,
enumerate all magnetic-space-group sites, assign canonical Bilbao onsite
axes/EBR labels automatically to Gaussian bands, or establish topological
classification of real materials. No new material calculation or runtime
feature is claimed by this validation extension.
