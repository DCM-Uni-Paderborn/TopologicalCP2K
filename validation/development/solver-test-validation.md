# Quantitative localizer and quadratic test evidence

The manuscript/SI update reanalyses retained runs from 24 September 2026.
It does not launch new CP2K calculations or attribute these results to the
current, partly uncommitted representation-development source.

## Reproduction and provenance

Run `python3 validation/development/summarize_solver_tests.py` from the paper
repository with NumPy available. The script reads `solver-test-records.tar.gz`
without extracting it and regenerates `solver-test-summary.json`. SHA-256
checksums identify the archive and every member. The summary verifies:

- Nine final QWZ timing records and their nonorthogonal comparisons, including
  the dense-gap/sparse-bracket agreement at the original driver's tolerance.
- Three accepted AII model records and their sparse solve-probe residuals.
- Resolved, unresolved and singular gap examples, without relabelling rejected
  cases as trivial indices.
- All four Bi2 dense/iterative GPW/GAPW output/state pairs: normalized states,
  metric-subspace singular values, residual decompositions, eigenvalues and
  reported eigenpair residuals.
- Primitive-mesh/supercell Ne torus gaps, indices and per-cell energies.

The input and output files in the archive are unchanged. They contain the
original development-version and compilation headers where CP2K printed them.
The model logs themselves have no executable fingerprint; the accompanying
`spectral_localizer_sparse_unittest.F` is a source snapshot at the time of this
paper update, not proof of an identical historical full build. The model index
checks are assertions in that driver, not separately printed output columns.
The existing `mumps.md` records the benchmark host, libraries and thread layout.
Those metadata are not inferred from the timing values.

The torus inputs retain their original absolute include path. To rerun them
elsewhere, point that include at the archived `torus-common.inc` in a compatible
CP2K checkout. Basis and potential filenames refer to the CP2K data directory;
this archive is a compact evidence record, not a complete frozen public
reproduction environment. No full basis/potential files or executables were
duplicated for this manuscript edit.

## What the numbers establish

### Dense and sparse class A

The table uses `localizer-final-{10,20,30}-{1,2,4}.log`, not the earlier
diagnostic/repetition logs. Each entry is one final run. The dense reference
uses one source rank, not a distributed dense eigensolver. MUMPS timings include
symbolic analysis, metric validation and all shifted gap factorizations, but
exclude fixture and coordinate-list construction. `factor_entries/input_entries`
is fill-in; `solver_memory_sum_mb` excludes CP2K and DBCSR allocations. Dense
oracle arrays are deliberately retained in the test. No process-RSS comparison
or multi-node scaling claim is supported.

All tested model sizes have index -1. At side 30, dense and four-rank MUMPS times
are 23.23512 and 1.36989 seconds. One-rank MUMPS instead takes 2.40084 seconds;
its reported memory sum grows from 19 to 45 MB when moving to four ranks.
At side 10, dense remains faster. These findings support optional backend
selection, not a universal sparse speedup.

### Class AII and threshold handling

`tacho-unit-final-large.log` contains accepted nontrivial, trivial and
nonorthogonal spin-mixed cases of complex localizer order 800. Their gaps are
0.79004778, 1.1109451 and 0.79004778 in model units, respectively. The separate
Pfaffian and MUMPS checks matter: ordinary inertia alone cannot distinguish
these AII indices. The source-rank solve-probe residual is at most 2.28023e-14
for the selected three cases, not a rigorous full-factor backward error.

Two zero-pivot messages in the Pfaffian log belong to intentional singular
inputs. The original driver checks their rejected status and proceeds to its
successful polynomial/permutation/scale/component result. The near-gap data
come from `periodic-sparse-unit.log`; despite that filename, the selected
fixtures are explicit finite matrix pencils, not periodic material runs.

### Periodic localizer

`torus-converged-{kpoints,supercell}.out` compares the same six-cell torus,
constructed from a primitive full mesh or a Gamma supercell. Both Ne/SOC
localizers have order 312 and Z2 index 0. Their generalized gaps differ by
4.27962e-10 Ha and the separately converged total energies, normalized to one
primitive cell, differ by 6.59982e-9 Ha. The latter is not zero and must not be
confused with a raw total-supercell/primitive energy comparison.

### Molecular quadratic states

Both GPW and GAPW are GTH-pseudopotential fixtures, not two independent
all-electron material predictions. Their recorded results coincide at printed
precision. The largest solver-to-solver gap difference is 5.99521e-15 Ha.
The iterative solver takes 145 Ritz steps and reaches a maximum eigenpair
residual of 8.21662e-10 Ha^2; the dense maximum is 2.82115e-15 Ha^2.
Recomputed metric normalization errors are at most 1.08802e-14, and the four
cross-subspace singular values differ from unity by less than 1e-10.
Near-unit singular values do not imply equally precise individual eigenvectors
inside a degenerate Kramers pair. The test compares the whole four-state space.

## Placement and exclusions

The main paper now summarizes the torus and Pfaffian comparisons and the
size-dependent solver crossover. The SI carries the full nine-row timing table,
model parameters, threshold outcomes, three-row AII table, torus settings and
Bi2 state diagnostics. Existing Kubo, Wilson and crystalline algebra results
retain their own provenance. No new 230-space-group material classification,
nontrivial converged DFT localizer benchmark, distributed Pfaffian factorization,
or current onsite-reference implementation is claimed by this update.
