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
  atomic signature matrices on explicit reciprocal sets
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

The native little-group source now has actual GPW/GAPW/SOC integration
tests and all-530-Hall-setting checks of projective characters,
antiunitary corepresentations and explicit-segment compatibility.
Native generation also covers 3,467 Wyckoff families and 5,648 specialization
relations across all 530 Hall settings, with an independent tabulated-coordinate
comparison. The generated sites are connected to the spinful site-induction
tests; this does not yet establish representation elementarity.
The combined atomic-signature layer has 2,120 checked matrices over ordinary/grey
and scalar/spinful variants on six explicit reciprocal points. Its 31,079
atomic columns satisfy the tested rank and compatibility conditions. Complete
reciprocal-stratum sampling and general integer/nonnegative classification are
not yet supplied by this layer.
The separate supplied-magnetic-group sweep tests algebra, not magnetic SCF.
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
