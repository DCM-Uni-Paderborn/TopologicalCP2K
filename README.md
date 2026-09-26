# Topological Systems with CP2K

Working manuscript for Computer Physics Communications, started 24 September 2026.

Repository: https://github.com/DCM-Uni-Paderborn/TopologicalCP2K
The repository is private during manuscript preparation. This is a manuscript
and validation repository, not the CP2K implementation repository.
The Dropbox/Overleaf project is the manuscript source of truth. Verified changes
are mirrored to GitHub without rewriting history; divergent remote edits must be
reconciled before publishing a new snapshot. The Git checkout is kept outside
Dropbox to avoid synchronizing a live Git database through file storage.

Authors, in the requested order:
Thomas D. Kuehne; Alexander Cerjan; Vladislav Efremkin;
Hermann Schulz-Baldes; Emil Prodan.
Corresponding author: Thomas D. Kuehne, tkuehne@cp2k.org.
Affiliations and submission declarations remain subject to author approval.

## Manuscript

- main.tex: title, abstract, introduction, conclusions and section includes.
- sections/operators.tex: nonorthogonal AO operators, Bloch links and SOC.
- sections/noncommutative.tex: Bellissard, covariance, trace per volume and Chern pairings.
- sections/kubo.tex: Liouvillian response, divided differences, currents and property symmetry.
- sections/band_topology.tex: native Wilson/Z2/Chern with retained neon/stanene figure.
- sections/space_group.tex: inversion TQC, native projective irreps, antiunitary
  corepresentations, explicit-segment compatibility, generated Wyckoff families, site-induced atomic bands,
  atomic signature matrices and exact signed/nonnegative integer decomposition
  and phasons.
- sections/localizers.tex: finite/periodic localizers, metric gap, LDL and Pfaffian.
- sections/quadratic.tex: finite/periodic quadratic pseudospectrum and matrix-free iteration.
- sections/applications.tex: TopoHSE-DB and the submitted Brocai et al.
  discovery study as application examples; small CP2K and companion Si/Al
  calculations as tests and demonstrations, with separate provenance.
- sections/implementation_validation.tex: numerical backends, evidence and scope.
- references.tex: cited bibliography.
- supporting_information.tex: retained and additional validation records.

The k-point manuscript retains Bloch/Wannier infrastructure, band structures,
densities of states, and metric-consistent Lowdin populations. Its detailed
topology sections, figure and topology-specific SI have been transferred
here. Interband dipoles are now related directly to the existing Bloch
velocity in sections/kubo.tex. sections/band_topology.tex connects the
implemented Berry-curvature expression to the subspace-based Wilson
description, including band degeneracies, phase conventions, and finite-basis
scope. These additions transfer the exposition, not a new validation dataset.
Calculation records are copied, not deleted from the historical source project.

## Build

Use main.tex as the Overleaf main document.
To build and check both documents locally, run `make check` with a TeX Live
installation providing `latexmk` on PATH. The generated PDFs are retained in the
repository; temporary build files stay in the ignored `.build` directory.
Alternatively, run latexmk -pdf main.tex and, separately,
latexmk -pdf supporting_information.tex.
The latest checked local build directory is outside the Dropbox project at
/Users/tkuehne/paper-revisions/topology-progress-20260925/bismuth-neutral-window/.
The earlier topology-progress-20260924 directory retains the pre-update manuscript and a compressed audit
of the byte-identical unused regression copies removed during cleanup.

## Evidence and boundaries

validation/topology contains the original Wilson, inversion and phason
records. validation/development contains implementation notes and the
Bloch-response test summary. These historical notes can describe earlier
scope limits than the current source; the manuscript's dated method and
validation sections are authoritative for this draft. Raw implementation notes are not new
numerical data and must not be substituted for material benchmarks.

The quantitative localizer/quadratic evidence is now also represented directly
in the manuscript and SI: nine dense/MUMPS model timings, finite-gap rejection
checks, three spin-mixed Pfaffian comparisons, a periodic Ne/SOC torus comparison,
and Bi2 GPW/GAPW dense/iterative state diagnostics. The compact raw-data archive
and the NumPy reanalysis script are described in
`validation/development/solver-test-validation.md`. The script regenerates
checksummed numerical summaries without rerunning CP2K. Single-run solver costs
are distinguished from material benchmarks and whole-process memory measurements.

The native little-group source now has actual GPW/GAPW/SOC integration
tests and all-530-Hall-setting checks of projective characters,
antiunitary corepresentations and explicit-segment compatibility.
Native generation also covers 3,467 Wyckoff families and 5,648 specialization
relations across all 530 Hall settings, with an independent tabulated-coordinate
comparison. The generated sites are connected to the spinful site-induction
tests. Real-space specialization now also produces onsite-induction
matrices, singleton equivalences and certified composite expansions.
All 22,592 relations pass 180,736 reciprocal character checks; independent
reference comparisons match exceptional-space-group coverage and archived
Bilbao EBR generator counts in all 2,120 variants. These are not yet
canonical-irrep-by-irrep identifications. Reproduction details are in
`validation/development/site-induction-validation.md`.
A subsequent affine-reference comparison matches all 1,731 named Wyckoff
families from independent generators/coordinates, including changed cell,
origin and operation conventions. Across four symmetry variants, all 6,924
site-resolved EBR dimension multisets agree (10,398 unreduced maximal-site
columns). This is not yet a canonical character match between individual
equal-dimensional onsite irreps. See
`validation/development/affine-reference-validation.md`.
The combined atomic-signature layer has 2,120 checked matrices over ordinary/grey
and scalar/spinful variants on six explicit reciprocal points. Its 31,079
atomic columns satisfy the tested rank and compatibility conditions. A separate
checked-integer kernel now factors these matrices and tests signed/nonnegative
membership, including independently checked inversion and small-matrix cases.
It supplies exact decisions or explicit unresolved resource-limit outcomes,
not complete reciprocal-stratum sampling or a general material classifier.
The `ATOMIC_SIGNATURES` keyword connects the reference and integer layers to
physically validated Gaussian-band characters. Unit tests exercise coset-phase,
operation and irrep-order changes; small GPW/GAPW/SOC calculations exercise
the physical integration, including nonmembership and an unresolved search.
The reciprocal sampling backend now generates fixed-torus closures, incidence
embeddings and primitive integer cycles. All Hall settings and grey extensions
pass geometry, coordinate-change and scalar/spinful induced-band checks.
`KPOINTS_SOURCE SYMMETRY` now connects this graph to Gaussian sampling with
subdivision, spectral-isolation and sewing checks. With atomic signatures,
the directed restrictions also define a certified integer compatibility
kernel and quotient by the atomic columns. Exported certificates agree with
independent arbitrary-precision Smith calculations. Free factors and incomplete
sampled graphs are not advertised as physical symmetry indicators. Star-arm
identifications now include spin factors, antiunitary conjugation and Bloch
phases. All 2,120 Hall-setting/symmetry combinations reproduce the published
compatible ranks and finite quotient factors for the 230 space groups. These
paths do not replace a full equivariant open-cell/connectivity construction.
The separate supplied-magnetic-group sweep tests algebra, not magnetic SCF.
The character-reference layer now matches 8,907 named irreps in 2,700
scalar/spinful tables covering all 230 groups. Explicit Bloch and actual
spin-lift conventions, central-projector checks, and full column bijections
are tested in serial, on both MPI processes and in an instrumented build;
the official local drivers each pass 173/173 assertions. This is reference
correspondence, not automatic standard labels in the Gaussian input path.
The Gaussian SOC sewing path now reuses these factor-system spin lifts.
A primitive BCC Ne doublet reproduces the previous product residual of 2
and passes with a residual below 9e-14 after the fix, without changing
physical inputs or symmetry tolerances. Including this test, the serial
and two-rank MPI drivers each pass 179/179 assertions; focused instrumented
spin/frame tests pass as well. See the shared-lift validation section in
the SI and `validation/development/gaussian-spin-lift-validation.md`.
The onsite matcher now retains explicit Cartesian axes and spin lifts.
Across 460 Gamma reference tables and three coordinate variants, all
1,380 complete-table matches pass; 88,602 reference-normalized characters
agree across serial, both MPI processes and an instrumented build within
8.44e-15, and across coordinate variants within 9.77e-15. Native column
numbering can differ between builds and is not a physical identifier.
This sweep removes space-group translations to realize point groups at
a site. Analytic C3, P6_3 and magnetic examples additionally test actual
onsite restrictions and nonsymmorphic induction. It is not a complete
real-Wyckoff-site EBR match or automatic Gaussian standard labeling.
Final serial/MPI drivers both remain 179/179. See the new SI subsection
and `validation/development/onsite-reference-validation.md` for scope,
raw records, reproduction scripts and precision conventions.
An additional actual-site audit retains all fractional translations and
generic representatives of all 1,731 conventional Wyckoff families. Four
symmetry variants and two coordinate conventions give 13,848 onsite-table
matches and 83,088 independently checked induced tables at six k points.
All 1,384,328 induced character entries agree across builds/coordinates
within 2.67e-14. The reference point-group axes are explicitly recorded,
not assumed to be canonical Bilbao onsite axes. This is real-site induction
validation, not complete canonical EBR naming or global band connectivity.
See `validation/development/wyckoff-character-validation.md` and the SI.
The production atomic-signature output also retains its onsite character,
atom-cell and physical spin conventions. Independent replay of 16 Gaussian
outputs per build reconstructs 31,537 induced character entries within
3.58e-14, and matches all 697 local columns to explicitly oriented external
point-group references. Serial and MPI drivers each pass 179/179 checks.
This preserves reference-label provenance, not canonical EBR naming; see
`validation/development/onsite-export-validation.md` and the SI.
A subsequent site-resolved audit compares complete tabulated character
sequences for 5,641 named ordinary-group scalar/spinful EBR columns in all
230 space groups. All columns match in two coordinate conventions and
across serial, two MPI processes and focused instrumentation; 378 retain
multiple candidates at the sampled points. Explicit type-I tables supply
the previously absent opposite-valley references. This is not unique onsite
labelling or a connectivity proof; see
`validation/development/ebr-character-validation.md` and the SI.
A full invariant-torus Fourier comparison now shows that all 378 of those
multi-candidate sets have equal unitary character functions, not merely
equal sampled characters. Denser k sampling alone cannot resolve them.
Conversely, analytic inversion centers and the Gaussian He tests show
that the new diagnostic does resolve genuine sampling aliases. It does
not merge atomic columns or infer Bloch-bundle equivalence. The all-530
Hall-setting audit and both 179/179 regression drivers pass. The SI adds
the quotient-lattice derivation and representative GPW/GAPW/SOC class
counts; full records and a verified extraction/replay are described in
`validation/development/fourier-character-validation.md`.
Periodic quadratic analysis now shares the full AO torus and ordered Berry
integrals with the localizer. Dense and two-rank iterative comparisons cover
He GPW/GAPW, Ne GTH-SOC, reduced SCF, query wrapping and a Gamma supercell.
Additional checks cover rigid skew-cell rotation, all periodic-axis choices,
equivalent MP/MACDONALD/GENERAL SCF meshes, and invalid input rejection.
The shared localizer/quadratic drivers pass 52/52 serial and 87/87 MPI assertions.
Raw data and a replayable comparison script are recorded in
`validation/development/periodic-quadratic-validation.md`. This is a
trigonometric position/energy diagnostic, not band unfolding or an
irreducible property solver. The manuscript and SI separate finite-basis
residuals, eigenpair accuracy and actual material-convergence requirements.
The complex-query extension now also checks nonnormal translations, explicit
adjoint products and complex expectations. An 18-site trimerized model reproduces
1,170 Bloch-block eigenvalues to 1.60e-14 in squared model energy units. An open
chain demonstrates why independent Hermitian real/imaginary squares are not
equivalent. The new distributed unit test passes with dense and MUMPS metric
solves; the affected regression drivers now pass 53/53 serial and 88/88 MPI.
Records, source fixtures and reproduction instructions are in
`validation/development/translation-quadratic-validation.md`. These are numerical
kernel tests, not Gaussian band unfolding.
The subsequent finite Gaussian adapter now evaluates shifted AO overlaps and
retains the continuum projection complement. Independent ghost-basis references
check He/H2/Ne GPW/GAPW spectra, complex means and leakage, with a maximum squared-gap
discrepancy of 6.1e-12 hartree squared. Analytic compressed-unitary tests, identity
translations, opposite momenta, four-rank all-electron/SOC comparisons and invalid
inputs supplement these checks. The rebuilt regression subset passes 57/57 serial
and 96/96 MPI assertions without changing earlier references. Main text and SI
distinguish the full translation residual from the projected energy residual.
Reproducible records are in `validation/development/gaussian-translation-validation.md`.
The periodic extension is documented separately in
`validation/development/torus-translation-validation.md`. It tests periodized
Gaussian translations against independent ghost-overlap matrices and explicit
Gamma supercells, including a corrected GAPW frozen-density handover. Primitive
3x1x1 meshes and supercells agree within 1e-10 hartree; the complete 18-state
lattice-translation spectrum also checks complex phases. These are implementation
checks, not converged material band-unfolding spectra.

The finite Bi neutral-window record now adds a 300/30 K scalar-SCF comparison,
complete-subspace frontier populations and a numerically bounded common
nontrivial localizer-scale interval. A fresh native 30 K export passes all
45 like-for-like index comparisons, and archive-only replay reproduces the
independent scans. Whole-process sampled RSS is distinguished from sparse
factorization buffers. These results validate the finite operator/numerical
paths; the small electronic gaps and strong boundary weights do not establish
a converged insulating bulk. See
`validation/development/bismuth-neutral-validation.md` and the corresponding SI
subsections for raw archives, reproduction instructions and remaining controls.
The separate `bismuth-termination-validation.md` record now tests removal of
only the two singly coordinated corners. The neutral gap grows from 7.735 to
219.6 meV, but frontier boundary weight remains high. Native and independent
localizers match at all 45 queries. A rejected wide kappa interval is retained
alongside the narrower, numerically bounded common nontrivial window. This
is a controlled termination test, not size convergence or edge relaxation.
The `bismuth-size-basis-validation.md` record adds a 16-atom TZVP calculation
and a separate 30-atom DZVP flake. Physical cross-basis overlaps check
occupied and frontier subspaces without identifying unrelated AO gauges.
The gap changes from 219.6 to 231.7 meV with basis, but drops to 72.14 meV
in the larger DZVP patch. Numerically bounded scale intervals can have
opposite indices at different sizes; a finer scan also finds a nontrivial
region missed by the coarse scale grid. All these differences are retained
as convergence evidence, not replaced by fitted reference values.
The `bismuth-temperature48-validation/` evidence adds a controlled 300 to
100 K comparison of the 48-atom patch. Its gap decreases by 2.91%, while
the physical occupied-subspace projector distance is 0.05181. Both sampled
index surveys agree; three joint parameter boxes are independently rechecked
at 100 K in `bismuth-temperature-region-validation/`. Three native sparse
queries match complete-AO references without changing comparison tolerances.
These endpoint checks do not establish the zero-temperature limit.
The optional `HALL_RESPONSE` adds the antisymmetric charge DC tensor to
the existing projected-AO/Bloch currents, using scalar relaxation.
Its nonzero model checks and native time-reversal controls are separate:
it does not introduce magnetic order or constitute a magnetic-material benchmark.
The source tree still has no full general EBR catalogue. Spin Hall response,
three-dimensional AII localizers, Gaussian band unfolding and multi-node scaling
remain open. The working draft is not submission-ready.

The stanene localizer data were subsequently corrected after a complete
Bloch-matrix comparison exposed a redundant SOC atom-block sign. The
current results and replay are in `validation/development/stanene-soc-validation.md`;
older localizer scans are explicitly superseded, while their Wilson
references remain valid. The corrected matrices agree within 1.85e-12 Ha,
but size/basis/scale dependence still prevents a material-localizer claim.
A separate complete-band spectral-flattening scan now compares DZVP tori
of sides 3, 6, 8 and 9 and a TZVP side-6 control, without changing the
projected coordinate links. At eta/Delta=0.75 and 1, every tested larger
torus has index one; small-size and eta/Delta=1.25 counterexamples remain
explicit in SI Table S15. This is not an asymptotic gap extrapolation.
The exact-sign material analysis and the separately tested DBCSR kernel
are distinguished in `validation/development/stanene-flattening-validation.md`.
The subsequent opt-in native adapter is validated separately against those
full-band references and finite GPW/GAPW/all-electron AO exports. It checks
the original positive metric and electronic gap before matrix-sign iteration,
and retains the distinction between electronic and rescaled localizer gaps.
The final regression subset passes 45/45 serial and 74/74 MPI assertions.
See SI S7.13 and `validation/development/native-flattening-validation.md` for
the source snapshot, native logs and hash-checked independent comparisons.
This numerical agreement does not remove material size/scale dependence.

The periodic flattening now uses complete primitive Bloch eigensystems after
verifying translation covariance, rather than iterating a whole-torus matrix
sign. A separate sparse complete-band reference extends the analysis to
12x12, 15x15 and 18x18 tori, plus a 12x12 TZVP and a finer-SCF-mesh control.
At the sampled eta/Delta values 1, 1.25 and 1.5, all three larger DZVP tori
give index one; their localizer gaps are still size dependent. Fifteen smaller
sparse/dense reference comparisons agree within 4.3e-15 Ha. Five other sparse
queries retain unresolved Pfaffian factorizations despite nonzero spectral-gap
estimates. Their failures and all ordering retries are preserved, not replaced
by expected integers. SI S7.14 and `validation/development/stanene-convergence/`
contain the larger-volume controls, source snapshots and hash-checked replay.
The sparse reference is not claimed as a new native production representation.

A separate delayed-skew prototype subsequently resolves all 33 queries without
symbolic restarts or relaxed tolerances. It retains every previously accepted
index, including the trivial small-torus counterexamples, and agrees with all
15 dense references. Two independent synthetic suites resolve 804 nonsingular
queries and reject 432 singular queries; 834 component-factor reconstructions
have relative errors below 1.04e-13. A standalone public-Tacho-interface build,
sanitizer checks, resource-boundary test and immutable-data replay are recorded
in `validation/development/pfaffian-delayed/`. Subsequent native integration is
recorded separately in `validation/development/pfaffian-native/`. It adds a
4 GiB default numerical-buffer limit and collective error handling. A Bismuth
SOC failure exposed one-sided roundoff in symbolic patterns; a consistent
real-skew projection, bounded relative to the independent physical gap,
repairs it without changing reference values. The native C++ implementation
passes 1,236 synthetic queries with runtime instrumentation, and the complete
localizer directories pass 74 checks with four MPI ranks plus 45 in the build
without optional sparse libraries. The native library also reproduces all 33
archived complete-band material queries at unchanged tolerance, with maximum
solve residual 4.07e-11; this does not rerun the large-system SCF calculations.
Numerical Pfaffian factors remain serial;
these are not multi-node scaling results.

Separately varied fresh SCF controls distinguish basis and cutoff effects
on a fixed 12x12 analysis torus, then refine the scalar SCF mesh, relative
grid cutoff and open-axis cell height independently. Their complete exports,
native-factor diagnostics and replay are retained in
`validation/development/stanene-separated-controls/`. The Supporting
Information table is generated directly from the retained results; its
electronic gaps are distinguished from the rescaled localizer gaps.
This is an independent complete-band reference analysis, not a new native
large-torus production representation or a complete-basis-limit claim.

The pre-split manuscript sources are backed up at
/Users/tkuehne/paper-revisions/topology-split-20260924/.
No journal submission or CP2K implementation push is implied by a manuscript update.

## Application References

TopoHSE-DB is cited with its published title and bibliographic record:
H. Mirhosseini, L. Elcoro, A. Knuepfer, T. D. Kuehne,
Machine Learning: Science and Technology 7, 040601 (2026),
DOI 10.1088/2632-2153/ae97d8. The citation is recorded in the
generative--predictive manuscript's local bibliography and was checked
against the publisher-deposited Crossref record on 24 September 2026.
The submitted study is by Luis Brocai, Thomas D. Kuehne,
Andreas Knuepfer and Hossein Mirhosseini; its title and reported counts
were checked against the local Overleaf sources. No journal acceptance,
DOI or public preprint is asserted for that submission.
Both original workflows use VASP and external TQC classification;
they are not represented as CP2K calculations.
