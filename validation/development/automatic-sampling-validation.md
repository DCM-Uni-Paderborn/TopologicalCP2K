# Automatic Gaussian symmetry sampling

Verified 24 September 2026. Research-branch results, not a released general TQC
classifier or a proof of global spectral isolation.

## Implementation

`KPOINTS_SOURCE SYMMETRY` uses the actual input-cell SPGLIB operations, including
nonsymmorphic/centering cosets and optionally the grey extension. It constructs
reciprocal fixed-family seeds, specialization incidences and primitive cycles,
then uniformly subdivides the lifted segments. Wrapped interior points are
deduplicated, but each new edge uses the difference of the endpoint cell offsets.
Thus total reciprocal winding is retained. `LITTLE_GROUP_PATH_POINTS` defaults
to three, including endpoints. `LITTLE_GROUP_MAX_POINTS` defaults to 4096 and
bounds the conservative refinement work before deduplication, not the preceding
geometric enumeration. Failure leaves the geometric sampling object unchanged.

At every generated point CP2K diagonalizes the frozen SCF Hamiltonian, then
checks selected/excluded band separation, metric covariance, subspace closure,
Hamiltonian commutation and the full sewing law. Native little-group and
compatibility analysis are implicit; `ATOMIC_SIGNATURES` is optional. The
existing signed/nonnegative membership decisions apply only to the sampled
signature matrix. `REQUIRE_GLOBAL_GAP` checks a sampled indirect gap for a
lowest-band prefix, not a continuum gap theorem.

The `.little_group` output begins with the family, cycle, point and shifted-edge
records, followed by the existing band representations and restrictions.
Eigenvalues remain in `.eig`. The native path creates no `.mmn`: its graph is
neither a Wannier neighbour mesh nor a Wilson surface. Existing NNKP, Wilson,
TRIM and Wannier export modes retain their original behavior.

## Gaussian and independent replay checks

Four fixtures use the same two-site screw geometry, with an appropriately
transformed second-atom coordinate in the sheared case:

| Case | Selected bands | Families | Points | Directed segments |
| --- | ---: | ---: | ---: | ---: |
| He/GPW, UZH DZVP MOLOPT/GTH | 2 | 15 | 56 | 82 |
| He/GAPW, UZH all-electron TZVPP | 2 | 15 | 56 | 82 |
| Ne/GTH-SOC, full scalar AO space in second variation | 16 | 15 | 56 | 82 |
| Unimodularly sheared He/GPW | 2 | 15 | 56 | 82 |

All 82 segments are compatible, including their antiunitary restrictions.
All four cases have verified signed and nonnegative atomic witnesses. The
file-level replay constructs explicit NNKP inputs from the native `AUTO_POINT`
and `AUTO_EDGE` records. To satisfy NNKP's rectangular neighbour format, it adds
distinct extra connections, then compares restrictions on the original edges.
It does not rely on a particular ordering of graph edges. This is differential
testing against the existing explicit-point path, not an independent DFT code.

The largest eigenvalue difference is 2.0962e-13 eV (rounded upward), and the
largest band-character difference is 2.1317e-14. Point/operation identities,
irrep/corepresentation counts, original segment restrictions, full atomic
matrix and signed/nonnegative witnesses agree. The integer records are compared
exactly. Serial/MPI eigenvalue differences are at most 2.0962e-13 eV.
Raw per-case numbers are in `automatic-sampling-results.json`.

With five points per original segment, He/GPW gives 138 points and 164 compatible
segments. Without time reversal it retains valid unitary restrictions and has
zero antiunitary segments. Four deliberate failures reject insufficient point
budget, one point per segment, a conflicting Wilson-loop request, and a single
screw branch that meets a degeneracy at a generated sample. No reference values
or tolerances are changed to hide these failures.

## Native and regression evidence

- 106/106 focused serial assertions and 106/106 two-rank MPI assertions:
  five native programs, 69 existing Gaussian assertions and 32 new assertions.
- The existing topology and k-point export directories additionally pass
  36/36 assertions in serial and 36/36 with two MPI ranks.
- The complete 530-Hall-setting ordinary/grey sweep (1,060 groups) passes with
  independent subdivision checks for every generated path. Original counts
  remain 18,700 closures, 57,357 incidences, 18,730 primitive cycles and 2,481,572
  induced-column compatibility comparisons.
- The same complete sweep passes with bounds checks, floating-point traps and
  undefined-behavior instrumentation on the reciprocal module and unit program;
  their dependencies use the existing CP2K library.
- Analytic subdivision tests cover interval counts 1 through 7, zero-length
  paths, positive/negative/multiple boundary crossings, exact accumulated winding,
  nonfinite/unwrapped inputs, invalid source indices and unchanged data on a
  budget failure.
- `make_pretty.sh` passes. All-electron, SOC and k-point comparison tests run
  with one BLAS thread and two OpenMP threads on Apple ARM/GNU Fortran 16.1.
  Test timing under concurrent local jobs is not a scaling benchmark.

## Reproduction and retained records

The native all-setting test is `topology_reciprocal_unittest --all-settings`.
Focused regtest directories are
`UNIT/topology_(integer|atomic|wyckoff|band|reciprocal)_unittest|QS/regtest-little-group`;
the exporter regression directories are `QS/regtest-topology|QS/regtest-kp-1$`.

The retained replay driver is `check_automatic_sampling.py`. For the local layout
used here, copy it as `build-serial/automatic-symmetry-smokes/validate.py` in the
CP2K checkout, copy `automatic-*.inp` and `automatic-common.inc` from
`tests/QS/regtest-little-group` into that directory, and run it with Python/NumPy.
It selects `build-serial/bin/cp2k.ssmp`, supplies `CP2K_DATA_DIR`, and preloads the
retained local thread-safe OpenBLAS. Adjust the executable/runtime paths when
porting the test. Do not reuse a directory with stale `.mmn` files for the native
seeds, since their absence is itself checked.

Compressed logs retain the serial/MPI focused and exporter drivers, full native
and instrumented sweeps, replay output and prettify output. The records archive
contains the fixture/replay inputs and symmetry/eigenvalue results, excluding
regenerable wavefunction files and large overlap exports.

## Boundaries

This integration does not create an equivariant open-cell decomposition,
classify complete global connectivity, establish EBR elementarity or standard
labels, or implement primitive-cell equivalence. Finite sampling can miss a
gap closure. The post-SCF SOC space must be converged independently. Retained
little-group/reference data and integer solvers are replicated on MPI ranks;
the tests are not a multi-node scalability claim.
