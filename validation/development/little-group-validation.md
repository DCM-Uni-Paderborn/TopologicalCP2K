# Native Little-Group Validation

Recorded 24 September 2026 in `/Users/tkuehne/cp2k-spectral-localizer`.
These are algebra and operator-integration tests, not a materials survey.

## Executable Evidence

- `little-group-serial-all-settings.log`: release serial, all 530 Hall settings.
- `little-group-mpi-all-settings.log`: release, two MPI ranks, the same suite.
- `little-group-debug-all-settings.log`: bounds checks and invalid/zero/overflow traps.
- `little-group-serial-regtests.log`: 47/47 Gaussian assertions, two OpenMP threads.
- `little-group-mpi-regtests.log`: 47/47 Gaussian assertions, two ranks x two threads.

Each all-setting run reports 12,720 factor cases, 6,360 complete unitary character
tables, 5,324 antiunitary endpoint tables, and 12,720 ordinary/grey segment
comparisons. Of the segments, 4,462 retain an antiunitary coset. Native analytic
tests include nonsymmorphic phases, physical spin factors, Kramers planes,
hourglass partner switching, monodromy, metric covariance and invalid counts.

## Earlier Wider Integration Pass

The endpoint-corepresentation stage also passed 110/110 serial and 111/111 MPI
checks across the topology, Kubo, localizer and quadratic directories. Those
totals and the 47/47 segment-focused assertions describe different test scopes.
The retained source logs are `build-serial/corepresentation-final-validation.log`
and `build-mpi/corepresentation-final-validation.log` in the development tree.

## Independent and Magnetic Algebra Checks

SPGLIB 2.7 supplied all 1,651 magnetic UNI settings. The six sampled k-points
and scalar/spinful factors give 19,812 cases, including 14,088 antiunitary
little-group tables and 5,724 cases without an antiunitary stabilizer.
Of those tables, the 10,154 primitive-cell cases accepted by the Spgrep 0.7.0
comparison agree in unitary restrictions, Wigner types and dimensions, with
maximum character residual 7.105e-15.

This reference comparison does not validate every explicit antiunitary matrix
returned by Spgrep: a separate product test found nonzero residuals in some
reference matrix corepresentations. Those matrices are never imported into CP2K.
The discrepancy has not been established as a library defect. Actual CP2K
Gaussian sewing matrices pass the full semilinear multiplication check before
classification. The reference sweep log remains in
`build-serial/corepresentation-magnetic-reference.log`.

## Boundaries

No complete EBR catalogue, conventional-label assignment, automatic complete
Brillouin-zone graph, magnetic symmetry search or magnetic noncollinear SCF is
claimed. Compatibility of sampled endpoints does not establish an insulating
gap between them. The papers cited as application examples used their own VASP
workflows and are not recast as CP2K calculations.
