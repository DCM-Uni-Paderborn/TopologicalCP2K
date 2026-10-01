# Source-Frozen Reproduction Checks

Completed 1 October 2026. This directory supports the corresponding SI section.
The original data and failed controls are retained, not replaced by new reference
values. None of these archives is a released CP2K distribution.

## Evidence

- `wilson-overlap-evidence.tar.gz` holds the historical Neon/Stanene inputs and
  overlaps. `wilson-frozen-build-evidence.tar.gz` contains the pinned-source
  reruns and scalar-space/thread controls. Both are needed for the comparison.
  The historical five-state Neon window cuts a degenerate triplet. The separate
  seven- and thirteen-state controls restore thread-independent fine spectra
  for these tested inputs. They do not establish material convergence.
- `localizer-cold-build-evidence.tar.gz` contains the original failing build
  records, the corrected source boundary, 12 successful unit executions,
  and 124/124 regression checks in each of two MPI/OpenMP configurations.
  Only the Ritz restart, GAPW scalar aliases, and quadratic unit-test file
  differ from the reconstructed snapshot. Runtime checks stay enabled.
- `pfaffian-frozen-build-evidence.tar.gz` records the instrumented native
  adapter rebuild, synthetic tests and 33 material-matrix replays. These use
  archived exports, not new material SCF calculations. Existing dependencies
  were fingerprinted, not rebuilt or instrumented. Pfaffian factors remain
  local to one rank.
- `topology-source-notices.tar.gz` preserves dependency notices and source
  correspondence. It does not assign a new license to the manuscript or
  analysis scripts. Unassigned notices and author licensing decisions remain
  explicit in `license-review.txt`.

The `*-results.txt` files are the dated records written for each independent
check. Statements about other blocks still being in progress refer to those
records' creation times. They do not supersede the consolidated SI.

## Verify Without CP2K

The first two commands use only the Python standard library.

```bash
python3 verify_localizer_build.py localizer-cold-build-evidence.tar.gz
python3 verify_pfaffian_bundle.py pfaffian-frozen-build-evidence.tar.gz
```

Wilson replay additionally needs the packages pinned in
`wilson-requirements.txt`. Run in a separate environment rather than altering
an existing CP2K or system installation.

```bash
python3 -m venv /tmp/topology-wilson-replay
/tmp/topology-wilson-replay/bin/python -m pip install -r wilson-requirements.txt
/tmp/topology-wilson-replay/bin/python -m unittest test_frozen.py test_replay.py
/tmp/topology-wilson-replay/bin/python verify_frozen.py \
  wilson-overlap-evidence.tar.gz wilson-frozen-build-evidence.tar.gz \
  --report /tmp/topology-wilson-comparison.json
```

These commands check recorded execution evidence and repeat the Wilson
matrix analysis. They do not launch a new CP2K build or SCF calculation.
The archived source-reconstruction patches, build options, dependency hashes,
and command records document those separate computations. Machine-local paths
in the historical commands must be adapted for a new installation.

`verification.json` records successful checks of these staged copies.
The separate negative-control records identify the deliberately altered
evidence rejected by each reader. Checksums establish correspondence to the
retained records, not cryptographic authenticity of scientific results.

No binary executables, object files, compiler caches, credentials, or Git
metadata are included. Source snapshots are method-specific, not a final
source freeze of every method discussed in the paper. Repository visibility
and author declarations are unchanged.
