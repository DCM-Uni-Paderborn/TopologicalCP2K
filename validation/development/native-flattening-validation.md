# Native metric-covariant spectral flattening

This record validates the opt-in `SPECTRAL_FLATTENING` adapter separately
from the exact-sign full-band material scan in `stanene-flattening.tar.gz`.
No earlier archive is overwritten or relabelled as native output.

## Method and controls

- The physical AO/spin Hamiltonian and metric are assembled before flattening.
  GTH SOC is included once. Positions and periodic Gaussian coordinate links
  remain unchanged; neither occupied-space projection nor coordinate-link
  unitarization is used.
- The paired original electronic pencil `diag(H-E*S,-H+E*S)` is tested for a
  positive metric and a resolved gap. Dense generalized diagonalization or
  MUMPS shifted inertia provides the check. Its half-signature must be zero.
- DBCSR evaluates `Delta*S*sign(S^-1*(H-E*S))`. Metric-inverse, sign-square,
  complex-structure and Hermiticity residuals must pass before any index is
  printed. The defaults are Delta=1 Ha, electronic-gap threshold 1e-8 Ha and
  sign-iteration tolerance 1e-9. The option is off by default.
- The operator is constructed once per energy and reused across positions and
  localization scales. Original SCF matrices are not modified. Distributed
  metric inversion and sign iterations can fill in; no linear-scaling claim
  follows. Dense preflight gathers matrices; sparse preflight does not.

## Numerical evidence

Two-rank stanene runs use the same inputs as the independent full-band scan.
DZVP at N=3 gives the expected indices 1, 0, 0 at eta/Delta=0.5, 1, 1.25;
the largest gap discrepancy is 5.7e-13 Ha. At N=6 the eta/Delta=0.75 and 1
gaps are 0.204263047335511 and 0.062924758531492 Ha, both with index one;
they differ from the complete Bloch reference by less than 2e-14 Ha.
The N=6 TZVP/400 Ry native control gives gaps 0.205441490284396 and
0.064209613692720 Ha and indices one at both scales. Its maximum gap
discrepancy is 2.1e-13 Ha, inverse residual 5.93e-11, and other residuals
below 8.1e-14. Across both bases all seven dense native material queries
agree within 5.7e-13 Ha with their independent full-band references.
The N=3 MUMPS/Tacho path retains both the trivial and nontrivial queries.
Both independent gaps are inside the sparse brackets, whose widths are below
9.2e-11 Ha. No dense fallback is used in that path.

Finite H2 controls use UZH bases/potentials, two energies, two positions and
two localization scales. All 32 queries are re-evaluated from exported H, S,
X and Y using independent generalized eigensystems. Maximum discrepancies:

| Calculation | Queries | Largest gap discrepancy (Ha) |
| --- | ---: | ---: |
| GPW, dense | 8 | 2.84e-12 |
| GAPW, dense | 8 | 2.89e-12 |
| All-electron GAPW, dense | 8 | 1.79e-12 |
| GPW, MUMPS lower bound | 8 | 4.70e-11 |

All molecular Chern indices are zero. The sparse brackets contain every
independent finite reference. Exactly two flattening operations occur per
eight-query run. Additional native controls reject an unresolved original
electronic gap, an unsuitable metric and a negative scale before an index.
Six constructed 12x12 complex pencils, including nonorthogonal metrics and
Kramers-mixed states, agree with exact spectral constructions within 6.1e-15
on serial and two/four MPI ranks. Negative tests include an indefinite metric
and a near-zero spectral mode. The four-rank sparse test uses a column grid.

The final directory driver passes 45/45 serial and 74/74 two-rank MPI
assertions, the latter including MUMPS and Tacho. Formatting passes all 54
selected files. Existing unflattened regression references are unchanged.

## Reproduction

Archive: `native-flattening.tar.gz`, 174 members plus manifest, 1,223,207 bytes.
SHA-256: `f3c5f9eea0e17a700afbf78944ef2f54f2898eb5337b46229cada13e2b1e96fc`.
The companion reference archive remains unchanged, with SHA-256
`4d372410f3eaec7d356078c64243665e6dca45f8e23cb0e6dc01fb37b65b9129`.
The replay verifies all 174 native and 70 reference members, compares nine
native material queries and recomputes all 32 finite-AO queries. The native
source snapshot is the local implementation commit `be85e0b8eb`.

The native archive contains inputs, complete text outputs, method source
snapshots, basis/potential files, CMake caches, unit and regression logs, and
a per-member SHA-256 manifest. Runtime banners and per-run commit metadata
may predate incremental rebuilding; the manifest identifies the tested
source snapshot. Runs overlap in wall time, so timings in their logs are
not used as benchmark measurements.

With Python 3.12+, NumPy and SciPy, recheck both archives and independently
recompute all finite-AO references:

```bash
OPENBLAS_NUM_THREADS=1 python validation/development/replay_native_flattening.py \
  validation/development/native-flattening.tar.gz \
  validation/development/stanene-flattening.tar.gz \
  /tmp/native-flattening-summary.json
```

This comparison reuses the retained independent full-band material values.
Those values can separately be regenerated from all exported Bloch states
using `replay_stanene_flattening.py --recompute all` as documented in the
preceding validation record. A CP2K executable is not needed for either replay.

Numerical cross-path agreement is not a bulk-material convergence result.
The small-torus/scale counterexamples, nonmonotonic gaps, longer-ranged sign
operator and nonunitary projected coordinates remain relevant. No general
finite-range theorem certificate, large-flake result or multi-node scaling
measurement is claimed.
