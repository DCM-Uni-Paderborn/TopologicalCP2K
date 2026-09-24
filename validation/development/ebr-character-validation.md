# Site-resolved EBR character validation

Validated on 25 September 2026 with local CP2K commit
`99cde3781dcce8eddee3159045b014ab8e44f724`. The commit adds an optional
native reference test and documentation. It changes no production
Hamiltonian, symmetry-reduction path, input keyword, reference energy or
tolerance. The CP2K commit has not been pushed upstream.

## Reference construction

`prepare_ebr_characters.py` combines the named ordinary-group scalar and
spinful EBR vectors in irreptables 3.1.0 with the separately validated
little-group character fixture. The affine fixture supplies conventional
Wyckoff families and the explicit primitive basis P. A generic seed is
selected in each referenced family; this is not a sweep over every point
in a free family. Labels and dimensions alone are not used as matches.

Every EBR basis label must have one character row with the stated
dimension. Every nonnegative integer vector is multiplied by the complete
local character table. All supplied k points must conserve its band
dimension. The result is an induced-character reference, independent of
the native atomic-band induction routine.

Twenty space groups require opposite-valley references absent from the
older character fixture. The adapter reads their explicit type-I tables,
not inferred complex conjugates: 40 scalar/spinful tables add 62 points,
including KA/HA, PA/WA. It verifies the full operation map, physical spin
lifts and integer coset shifts, then compares every shared labelled
character before importing the additional points. The largest shared-table
residual is 3.0044025837014233e-16. No phase is fitted separately to an
irrep. The convention remains k_native=-transpose(P)*k_table; printed
characters retain their source rounding. The fixture, original EBR JSON
files, supplementary type-I tables and their hashes are archived.

## Native comparison

`topology_band_unittest --ebr-reference=FILE` constructs all local columns
at each supplied site. The test compares every unitary induced character
against each named reference. A candidate must match the entire k-point
sequence, not a different reference at each point. Original and combined
Cartesian rotation, lattice shear, origin shift, operation reversal and
integer coset-change conventions are tested. Spin-lift and Bloch-coset
phases are removed explicitly before the comparison.

All matching native columns are retained. Identical candidate sets are
compared as multisets, preventing multiple duplicated reference columns
from reusing a single native column. Overlapping but unequal tolerance
neighborhoods, incomplete little groups and malformed records are rejected.
The character bound is 2e-4 after scaling each column by its own positive
band dimension. The native construction tolerance remains 1e-9.

## Results

Each executable/rank checks both coordinate conventions:

| Quantity, per convention | Result |
| --- | ---: |
| Ordinary space groups | 230 |
| Site/spin cases | 1,599 |
| Named reference columns | 5,641 |
| Reciprocal queries | 9,667 |
| Induced character entries | 297,529 |
| Unique candidates within their specified site | 5,263 |
| References with multiple candidates | 378 |
| Native local columns at the selected sites | 5,704 |
| Largest dimension-scaled character residual | 1.0471976e-5 |

Optimized SSMP, both processes of a two-rank PSMP run, and a focused
instrumented executable all pass and agree on the candidate multiplicities
for every label and convention. The MPI processes independently repeat
the algebra; this is not a speedup or distribution benchmark. Two OpenMP
threads, one OpenBLAS thread and OMP_STACKSIZE=64M were used.

The instrumented recipe recompiles the test, symmetry, little-group,
corepresentation, character-matching and induction kernels with
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined
-fno-sanitize-recover=all`. Other modules are linked from the existing
serial release library. This is not a full instrumented Gaussian/MPI or
standalone no-SPGLIB build. The separate default native driver passes
8/8 checks in each optimized build. The earlier 179/179 Gaussian drivers
are not counted as new DFT runs for this test-only extension.

`check_ebr_rejections.py` runs an unchanged P4 spinful baseline and six
negative inputs against both optimized executables. The following fail
with the expected diagnostic: altered nonidentity boundary character;
shifted seed with unchanged Gamma reference; duplicate character column;
zero spin lift; duplicate operation; trailing data. The last two also test
input completeness. The baseline passes. These are 12 successful negative
checks, not failed physical calculations.

## Interpretation

This is a character-sequence comparison to named external EBR references.
It is stronger than dimensions, Gamma-only matching or separately permuting
the columns at different points. It is not a proof of canonical onsite
axes, EBR identity as a Bloch bundle, full connectivity or occupied-band
topology. In particular, the 378 multi-candidate references remain
explicitly ambiguous. The 63 native local columns not represented in the
supplied EBR data are not assigned an elementarity verdict by this test.
Grey-group EBR reference matching and automatic reference labels in the
Gaussian analysis remain separate work. No external table is a CP2K
runtime dependency.

## Evidence and replay

`ebr-character-records.tar.gz` contains 319 files (27,484,519 bytes before
compression): fixtures, original reference data, manifest, logs,
instrumentation/formatting records, negative inputs, changed source files
and the local commit patch. Every archived file was checked byte-for-byte
against its source using SHA-256. See `ebr-character-files.sha256`,
`ebr-character-archive-check.json` and `ebr-character-comparison.json`.
This is targeted validation evidence, not a complete release source tree.
After extraction into a fresh scratch directory, all 319 hashes were
checked again and the archived-log audit was byte-identical to the
original summary. A fresh native serial run with the extracted fixture
also reproduced every candidate multiplicity. The replay record is
`ebr-character-replay.json`.

Extract the archive into a scratch directory. The Python audit uses only
the standard library; invoke the script from this development directory:

```sh
tar -xzf ebr-character-records.tar.gz
shasum -a256 -c /path/to/development/ebr-character-files.sha256
python3 /path/to/development/audit_ebr_characters.py \
  build-serial/ebr-character-full.json \
  build-serial/ebr-character-verified-serial.log \
  build-mpi/ebr-character-full.log \
  build-serial/ebr-character-instrumented.log
```

To rerun the native calculations, use a matching CP2K checkout and built
executables, with the archived fixture supplied by absolute path:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  build-serial/bin/topology_band_unittest.ssmp --ebr-reference=/path/to/ebr-character-full.dat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  mpiexec --tag-output -n 2 build-mpi/bin/topology_band_unittest.psmp \
  --ebr-reference=/path/to/ebr-character-full.dat
```

`build_ebr_instrumented.sh` rebuilds the focused instrumented executable
from that checkout. `check_ebr_rejections.py BINARY FIXTURE OUTPUT_DIRECTORY`
replays the negative checks. Regenerating the reference requires NumPy,
SciPy, SymPy and irreptables 3.1.0, plus the sibling
`prepare_wyckoff_characters.py`; exact versions are in the archived
environment record. The generator accepts the character fixture, affine
fixture, output fixture and manifest, followed optionally by space-group
numbers. The reference package's full character/EBR data, not any CP2K
induction output, provides the numerical reference.
