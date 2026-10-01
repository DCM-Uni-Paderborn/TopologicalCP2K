CONCENTRIC 96-ATOM BISMUTH CONTROL

This record accompanies SI S12. It concerns a bare, unrelaxed, isolated
Bi(111) patch with scalar restricted PBE followed by GTH SOC. It is not
self-consistent spinor DFT, an edge-relaxation study, or bulk convergence.

The 48- and 96-atom DZVP patches share the same central 48 atomic positions
and Gaussian definitions. The geometry, neutral spectral gaps, complete
frontier doublets, physical shared-AO comparisons and sampled localizers
are retained without adjusting their values to a desired classification.
The first transferred joint energy/position/scale box is numerically
resolved as nontrivial. The middle-scale interval has distinct gapped
endpoint indices and must not be described as a uniform nontrivial window.

CONTENTS AND PROVENANCE

The numbered .tar.gz.partNNN files are sequential byte segments of one
gzip-compressed tar archive. index.json records every part's size/hash,
the combined archive hash and the hashes of the selected archived files.
All native run.json records retain their original, more extensive file
inventory. Only the subset needed for this finite-system analysis is
archived here. Missing .mmn/.win/.eig products in that inventory are
deliberate exclusions, not missing inputs to the reported analysis.

The bundle contains native inputs, outputs, full spinor exports, relevant
96-atom restarts, Gaussian basis/potential data, linked-library provenance,
build configuration, method scripts and analytic tests. The qualified
source archive and source-file hashes match the separately documented
Linux platform controls. CP2K build products are not included. Moment
caches are omitted because the replay reintegrates the printed Gaussians.
Pending TZVP results are not included.
The source BASIS_MOLOPT_UZH symlink is archived as its regular-file
contents, with the original versioned target recorded in the manifest.
The small region-transfer supplement was archived after all three box
tests completed. Its index binds it to the unchanged main archive hash.
It retains all 51 complete-spectrum center/corner checks. The low-scale
box has index one with lower gap bound 2.0727987627e-4 hartree, the
middle box has both gapped indices, and the high-scale box has index zero
with lower gap bound 5.0390483171e-4 hartree. No uniform bound is assigned
to the middle box.

VERIFY WITHOUT EXTRACTING

Run with Python 3.11 or later from this directory:

  python verify_size_checks.py .
  python verify_size_checks.py region-transfer
  python -m unittest -v test_size_archive

The verifier needs only the standard library. It checks part order and
hashes, archive-member paths, native file fingerprints and analysis
provenance. This integrity check is not a new scientific calculation.

INDEPENDENT NUMERICAL REPLAY

With Python 3.12 or later, NumPy and SciPy installed, choose a new output directory on a machine
with several GiB of available memory and disk space:

  OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=1 python replay_size_checks.py \
    . /absolute/new/path/bismuth-size96-replay

The reader extracts the verified archive into that new directory,
reintegrates the complete AO moments and repeats the nine spectrum queries,
the three native-run queries and the separate transition witness. It also
repeats the first box's norm bound, adaptive cover and anchor index, then
the shared-subspace comparison and analytic tests. The 17 original
complete-spectrum center/corner audits are retained in box-0.json. They are
not all repeated by this shorter replay. No SCF is rerun.

The native-comparison threshold is unchanged at 2e-8 hartree. Replay of
identical mathematics requires agreement within 1e-10 hartree, not a change
of the physical comparison tolerance. replay-report.json records actual
errors and distinguishes fresh replay from the retained native results.

The retained replay completed successfully with Python 3.13.5, NumPy
2.4.2 and SciPy 1.17.0. All 13 query gaps and indices, the shared-subspace
comparison and the first box's lower bound reproduce without numerical
change. All 17 analytic/archive tests pass. The other two box analyses
are completed original calculations, not part of this shorter replay.
To rerun the entire box analysis, use the main archive's extracted
analyze_material_regions.py with its flake96-spectrum directory and a new
output directory. This repeats the full 51-spectrum audit and all bounds.

INTERPRETATION

The shared-space projector is the compression of a zero-temperature
neutral SOC spectral projector, not the thermal SCF density. Its trace
is occupation projected onto a Gaussian subspace, not an atomic charge or
the number of electrons inside a sharp spatial region. It need not be
idempotent. The norm comparisons use exactly shared physical functions,
not coefficient padding or unrelated Lowdin gauges.

All timings are individual runs. Process-tree RSS sums can count shared
pages more than once. The three-query native run was briefly debugger-
sampled and is not a clean timing benchmark. Concurrent analysis and basis
control also preclude interpreting the transition run as a scaling test.
These results do not establish basis, cutoff, size or temperature convergence.
