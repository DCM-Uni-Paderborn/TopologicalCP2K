# Star identifications in the integer compatibility lattice

Verified 24 September 2026 in the local CP2K research branch. The earlier
`compatibility-lattice-validation.md` and its records are a retained baseline,
not the current star-aware output. No CP2K pull request is published by this update.

## Mathematics and tests

For s taking k_a to k_b+G and h_a=s^-1 h s, the representative difference is
L=tau_h+W_h tau_s-tau_s-W_s tau_ha. The transported character is multiplied by
z(h,s)/z(s,h_a)*exp(-2*pi*i*k_b.L) and is conjugated for antiunitary s.
The cached Gamma factor z includes the physical spin and time-reversal factor.
Both full product conjugation and the twisted factor identity are checked.
Unitary irrep maps, and the full corepresentation maps reconstructed from their
unitary restrictions, must be bijective and dimension-preserving.

Star equations M_s*n_a-n_b=0 enter the same integer matrix C as line restrictions.
A spanning set suffices; edge reciprocal winding remains a separate constraint.
Full-group input operation IDs and reciprocal shifts are exported for each link.
No new keyword, external character table or runtime solver dependency is required.

Tests include a complex threefold screw with time reversal, separately shifted
endpoint representatives, reciprocal images, corrupted factors, nonfinite
characters and unsafe coordinates. A star-only graph with no line edges reduces
nine local coordinates to three and rejects a mismatched target. For every link
in all 530 ordinary/grey Hall settings, scalar and spinful, the map agrees with
independently induced atomic columns. Its inverse and every alternative operation
connecting the same pair yield the same full-corepresentation permutation.

## Independent published reference

`compare_star_reference.py` reads the layout-preserving text extraction of
arXiv:1703.00911v2, Tables I-IV and IX-XII. Each rank table must contain all 230
space groups exactly once. Unlisted quotient entries are trivial, as in the
reference. SPGLIB 2.7.0 metadata maps all 530 Hall IDs to space-group numbers.
The native test output is parsed separately; complete coverage and a successful
test terminator are required. The script contains no fitting or expected-value
correction based on the native results.

All 2,120 comparisons match BOTH the compatible rank and the ordered nontrivial
invariant factors; all these automatic graphs have free quotient rank zero.
The four settings are ordinary scalar, ordinary spinful, grey scalar and grey
spinful. These are algebra tests, not 2,120 electronic-structure calculations.
The corresponding JSONL records retain the observed and reference values.
This agreement does not prove EBR elementarity, canonical labels, arbitrary
user-graph completeness or physical band isolation.

Reproduce using the cited authors' PDF (not redistributed here):

```sh
pdftotext -layout symmetry-indicators-reference.pdf reference.txt
python3 compare_star_reference.py reference.txt all-settings.log /path/to/libsymspg.so
```

On macOS the shared library suffix is `.dylib`. The small ctypes structure follows
the SPGLIB 2.7 public header; adapt it if that public ABI changes. This reader is
development validation, not a CP2K runtime dependency.

## Gaussian and parallel checks

The official drivers pass 170/170 assertions in SSMP and two-rank MPI: five native
unit programs, 129 little-group assertions and 36 pre-existing exporter/k-point
assertions. Five added assertions count star links. Existing physical reference
values are unchanged. The latest MPI analytic run additionally exercises the
star-only graph added after the first MPI driver was started.
The complete 530-setting output is byte-identical between optimized SSMP and
the instrumented executable. After removing rank prefixes, each MPI rank gives
the same 7,438 lines. The instrumented build covers the integer, star-transport,
atomic-signature, compatibility-lattice and reciprocal modules plus the unit
program, with `-fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined`.
Other dependencies use the existing CP2K library. Its database-free analytic
variant also passes; this is not a complete no-SPGLIB build.

`verify_star_quotients.py`, using SymPy 1.14.0, independently verifies eleven
Gaussian certificates per build. It checks rank(C), saturated Smith(K), Smith(Y),
C*K=0, K*Y=A, C*A=0, unimodular quotient maps, their integer image and target
classes. It also checks the reciprocal geometry of each exported star link with
exact rational arithmetic. Serial and MPI summaries are identical.

| Fixture | Equations | Coordinates | Kernel rank | Free rank | Finite factors | Star links |
|---|---:|---:|---:|---:|---|---:|
| Automatic He/GPW, He/GAPW, sheared He | 109 | 91 | 6 | 0 | 2 | 0 |
| Automatic Ne/SOC | 81 | 67 | 5 | 0 | 2, 2, 4 | 0 |
| Explicit He screw, GAPW screw | 10 | 9 | 3 | 0 | none | 1 |
| Explicit Ne/SOC screw | 5 | 5 | 2 | 0 | none | 1 |
| Partial He, sheared He | 1 | 8 | 8 | 4 | none | 0 |
| Refined He Wilson | 76 | 75 | 23 | 4 | none | 6 |
| Single He screw branch | 1 | 2 | 1 | 0 | none | 0 |

All targets except the single screw branch are compatible and class zero. That
branch remains incompatible and has no class. The partial graphs are still not
complete indicators: four free factors remain in the Wilson graph, down from 23
in the baseline. Its kernel rank changes from 42 to 23. No free factor is assigned
a physical topological interpretation.

Use the restriction
`UNIT/topology_(integer|atomic|wyckoff|band|reciprocal)_unittest|QS/regtest-little-group|QS/regtest-topology$|QS/regtest-kp-1$`
with the official test driver. SSMP uses one rank/two OpenMP threads; PSMP uses
two ranks/two threads. Gaussian runs use one BLAS thread, OMP_STACKSIZE=128M and
the retained thread-safe OpenBLAS preload. Timings under concurrent local tests
are not scaling benchmarks. All integer algebra is replicated, not distributed.

Compressed outputs and certificates use the `star-transport-` prefix. No SCF
restart, executable, downloaded reference PDF or large build is included in this
paper repository.
