# Controlled electronic-temperature comparison of the 48-atom Bi flake

The 100 K calculation changes only `ELECTRONIC_TEMPERATURE` relative to
the retained 300 K spectrum run. Geometry, cell, UZH DZVP basis, 400/40 Ry
grids, post-SCF GTH SOC, `EPS_SCF 1e-9`, executable, shared library and
runtime are unchanged. The initial guess is the converged 300 K scalar
wavefunction. The comparison checks actual input text as well as recorded
options and hashes; a simultaneous basis or cutoff change is rejected.

All 624 scalar AOs and 1248 spinors are retained. This remains a bare,
unrelaxed finite flake, not a size-converged bulk classification. Smearing
controls the scalar-reference SCF, not ionic temperature or a finite-T Z2
definition. No native CP2K implementation was changed in this comparison.

## Results

| Quantity | 300 K | 100 K |
| --- | ---: | ---: |
| Neutral SOC gap (meV) | 37.9151043 | 36.8113022 |
| Gap / kBT | 1.4666201 | 4.2717695 |
| Occupied frontier boundary weight | 0.6514272 | 0.6586375 |
| Empty frontier boundary weight | 0.6448893 | 0.6417633 |
| Nontrivial queries in the 45-query survey | 7 | 7 |

Every paired survey index and every paired point in the separate 20-scale
centroid scan agrees. Energies are fractions of each case's own neutral
gap, not identical absolute energies. The scans are not continuous bounds.

The occupied-subspace projector distance is 0.0518125620, with mean
retained occupied weight 0.9999553044. Occupied/empty frontier projector
distances are 0.0508407063 and 0.0471874816. These comparisons retain whole
degenerate subspaces and use physical cross-AO overlaps, not eigenvector
differences. The absolute Hamiltonian change is 8.76906114e-4 Ha. Subtracting
the neutral-midpoint shift times the AO metric gives a physical operator
norm of 8.65657911e-4 Ha. An unchanged-export control gives zero Hamiltonian
change, with the expected square-root roundoff floor in principal angles.

Three native Tacho/MUMPS queries at kappa = 0.0003, 0.0015, 0.003 Ha/bohr
have indices 1, 1, 0. Independent gaps from this restart's fresh spectrum
all lie inside the native brackets, whose widths are below 9.01e-11 Ha.
The comparison tolerance remains 2e-8 Ha. Maximum Pfaffian solve residual
is 1.50811e-13, with no workspace retries. This is three native queries,
not a native evaluation of the entire 45-query survey.

The spectrum run converges in 172 steps and takes 4609.68 s. The six-step
native restart takes 553.13 s, including SCF and export. Sampled process-tree
RSS peaks are 20.84 and 18.95 GiB, respectively; these are not full physical
footprints or factorization-only benchmarks. Both runs finish without warnings.

## Retained evidence and reproduction

`index.json` hashes the two complete native run archives and frozen methods.
`summary.json` contains the original native metadata and independent reports.
The methods archive retains the controlled comparison, both surveys,
unchanged-export control and 37 method-test results. The existing 300 K
archive is referenced, not recomputed or duplicated in this bundle.

Archive-only replay verifies 227 members, reintegrates both Gaussian moment
caches and reproduces every analysis, survey and temperature comparison
without numerical change on the recorded local stack. It does not rerun SCF.
From `validation/development`, with NumPy/SciPy available:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python replay_bismuth_material.py \
  bismuth-temperature48-validation ../../.build/temperature-replayed.json \
  --scratch ../../.build --reference-bundles bismuth-flake48-validation
```

The output path must not already exist. Continuous energy/position/scale
boxes are retained separately in `../bismuth-temperature-region-validation/`.
These two temperatures do not establish the zero-temperature limit or bound
a continuous temperature interpolation. The finite-size and basis controls
remain separate and are not superseded.
