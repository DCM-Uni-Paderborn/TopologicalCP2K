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
- sections/quadratic.tex: quadratic pseudospectrum and matrix-free iteration.
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
The latest checked local build directories are outside the Dropbox project under
/Users/tkuehne/paper-revisions/topology-progress-20260924/.
That directory also retains the pre-update manuscript and a compressed audit
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
The source tree still has no full general EBR catalogue. Hall response,
three-dimensional AII localizers, periodic quadratic analysis and
multi-node scaling remain open. The working draft is not submission-ready.

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
