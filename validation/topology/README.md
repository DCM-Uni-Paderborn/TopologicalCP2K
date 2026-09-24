# Native Wilson-Loop Figure and Evidence

Added to the k-point paper on 17 September 2026. These are retained results
from the topology development task, not newly executed DFT calculations in
the manuscript-editing task.

## Files and Provenance

- `neon/`: unmodified input, completed CP2K log and final native Wilson spectrum.
- `stanene/`: corresponding input, log, final spectrum and native/Z2Pack comparison report.
- `checks/`: copied Neon, adaptive Stanene, geometry and multicentre reports.
- `plot_wilson.py`: validates completed runs, dimensions, center ranges and final invariants before generating the figure.
- `figure-summary.json`: source-file SHA-256 checksums and figure sampling metadata.
- `*/figure-data.csv`: full-precision numerical data used for the two panels.
- `source-inspection/`: copies of the relevant source modules inspected for the manuscript, with SHA-256 checksums. This is an inspection record, not a claim that the changing worktree exactly reproduces the earlier executable.

The underlying calculations are in
`/Users/tkuehne/work/cp2k-z2pack-validation/final-serial/soc-native` and
`/Users/tkuehne/work/cp2k-z2pack-validation/stanene-converged`.
The independently adaptive comparison is in
`/Users/tkuehne/work/cp2k-z2pack-validation/stanene-z2pack-verified`.
Geometry checks come from `geometry-final` and `multicentre-final` beneath
the same validation root. These local locations are provenance pointers,
not public data-repository identifiers.

The checkout `/Users/tkuehne/work/cp2k-z2pack`, branch
`z2pack-overlap-loops`, is based on
`3919fb7b1c1b8804bc2b603e2cadadc551023053`. The extension is uncommitted
at inspection. The outputs identify this base hash, which by itself does
not identify the added code. The mainline Si/Al source and optional
Wannier90-library PR are separate implementations and calculation sets.
The final integrated source, compiler/build configuration and matching
reruns must be frozen together before submission; see the author checklist.

## Regeneration

Run from any directory with Python, NumPy and Matplotlib:

```sh
python3 /path/to/paper/validation/topology/plot_wilson.py
```

The script writes the vector PDF and a preview PNG under `figures/`, as
well as CSV data and the checksum manifest. It does not rerun CP2K,
refit curves, jitter degenerate centers or change any retained spectrum.
The transverse coordinates follow the input defaults verified in
`input_cp2k_print_dft.F`: origin `(0,0,0)`, loop winding `(1,0,0)` and
transverse vector `(0,0.5,0)`. Both endpoint loops are retained.

Neon has 5 loops with 8 points each, and Stanene 193 loops with 192 points
each. Both contain 8 occupied spinor states. Results refer to the specified
finite basis and cutoffs, not a new fully converged material prediction.
Native/Z2Pack agreement on identical overlaps tests the loop algebra;
operator identities provide complementary checks of the AO representation.

The complete `.mmn` matrices and adaptive Z2Pack checkpoint remain in the
original validation directories. They are not needed to reproduce this
figure, but must be included in the final deposited numerical archive to
recompute the full cross-implementation comparison. No scaling data,
automatic four-index 3D classification or noncollinear self-consistent
calculation is inferred from this evidence.
