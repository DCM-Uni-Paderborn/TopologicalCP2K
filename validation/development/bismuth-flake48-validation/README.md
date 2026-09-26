# Finite 48-atom Bi control

These are fixed, bare, unrelaxed Bi(111) flakes, not a converged bulk
classification. The five-by-five primitive patch loses its two singly
coordinated corners, leaving 48 atoms, 240 electrons, 624 UZH DZVP scalar
AOs and 1,248 post-SCF GTH-SOC spinors. Scalar restricted PBE uses 400/40 Ry
grids, 20 angstrom vacuum padding, 300 K smearing and EPS_SCF 1e-9.

## Results and scope

- The neutral electronic gap is 0.0013933545044090834 hartree (37.915 meV),
  compared with 219.6 and 72.14 meV in the earlier 16- and 30-atom DZVP
  controls. It is only 1.4666 times the scalar SCF thermal energy.
- Boundary atoms make up 37.5% of the patch; complete occupied and empty
  frontier doublets carry 65.14% and 64.49% of their Lowdin weight there.
- The independent 45-query survey contains seven nontrivial results.
  The separate 20-scale survey shows that the 30-atom nontrivial windows
  cannot simply be transferred to this larger flake. Individual sampled
  indices are not proofs of a continuous window.
- A separate native restart checks exactly three centroid queries at
  kappa = 0.0003, 0.0015 and 0.003 hartree/bohr. All indices (1, 1, 0)
  agree with independent full-AO reconstructions from its own fresh SOC
  export. Every independent gap lies inside the native MUMPS bracket;
  the comparison tolerance remains 2e-8 hartree. This is not a 45-query
  native test.
- Full energy/position/scale boxes, including failed transfers from the
  smaller flake, are retained separately in the sibling
  `bismuth-region48-validation` bundle. They reuse these complete snapshots.

## Runtime qualifications

The four-rank attempt was deliberately stopped under memory pressure after
step 25, not completed successfully. Its restart was resumed at one MPI
rank and two OpenMP threads, without changing physical or convergence
settings. The final 91-step spectrum stage took 2310.61 s; that excludes
earlier interrupted work. The separate six-step native restart took
477.44 s, including three localizers and a complete SOC export.

Process-tree sampled RSS maxima were 20.99 and 21.76 GiB. The separate
macOS spectrum-run snapshot reported a 42.5 GiB footprint including swapped
memory, with 45.2 GiB peak. RSS is therefore not a complete memory bound.
Neither these runs nor the comparison with a 12.63 s dense postprocessing
analysis constitute a solver or MPI-scaling benchmark.

## Retained evidence

`index.json` identifies all three checksummed archives; `summary.json`
retains the original run and independent-analysis metadata. Each successful
case includes its input, original runner, initial restart, full eigenstate
export, analytic AO moments, independent method and native output.
`methods-results.tar.gz` contains the shared methods, basis/potential
data, the 32-test log, both independent scans, the deliberate stop record,
the interrupted run/output and the two memory snapshots. Failed resource
attempts are never counted as scientific results.

`replay.json` records an archive-only replay. Both AO-moment caches are
reintegrated from the retained Gaussian contractions rather than reused.
All moments, independent localizers, native comparisons, scans, frontier
weights and report metadata reproduce with zero numeric difference.
No SCF is rerun. From the manuscript repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  validation/development/replay_bismuth_material.py \
  validation/development/bismuth-flake48-validation \
  .build/flake48-fresh-replay.json --scratch .build
```

Use a new output path. The Python environment needs NumPy and SciPy.
The replay method itself is versioned alongside this bundle, and its
SHA-256 is recorded in `replay.json`.
