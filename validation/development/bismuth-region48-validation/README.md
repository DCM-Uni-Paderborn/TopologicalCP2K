# Joint finite-Bi regions at 48 atoms

This bundle retains the same joint energy/position/scale analysis as the
earlier 16- and 30-atom controls, using the complete 48-atom DZVP spectrum
in `../bismuth-flake48-validation`. The Hamiltonian stays frozen in every
box; these are analysis-parameter changes, not new SCF calculations.

All boxes use energy midpoint +/- 2% of the full neutral electronic gap
and centroid +/- 0.05 angstrom in each in-plane direction. Kappa is in
hartree/bohr. The position norm is 29.354648338606204 bohr. The gap margin
remains 1e-10 hartree; no numerical tolerance is relaxed.

| Kappa interval | Status | Index | Evaluations | Covered/pending leaves | Minimum covered-leaf bound (Ha) |
| --- | --- | --- | --- | --- | --- |
| 0.003--0.0035 | resolved | 0 | 5 | 3/0 | 7.6328155856e-4 |
| 0.00025--0.0004 | nonuniform | none | 15 | 6/4 | not a whole-box bound |
| 0.0015--0.0025 | nonuniform | none | 15 | 7/2 | not a whole-box bound |
| 0.00025--0.0003 | resolved | 1 | 9 | 5/0 | 6.1656232760e-6 |
| 0.0015--0.002 | resolved | 1 | 13 | 7/0 | 6.1624392231e-5 |

The two broader lower-scale boxes worked at 30 atoms but have both gapped
indices among their corner checks at 48 atoms. This proves nonuniformity;
it is not merely a consequence of a small subdivision budget. No common
index is assigned to either failed transfer. The narrower positive boxes
were chosen using the retained 20-scale scan, and verified separately.

For all five boxes, center and corner audits independently check full
divide-and-conquer spectra, spectral pairing, the central-four-eigenvalue
gap, actual perturbation norms and Pfaffians. The adaptive cover, not
corner sampling, establishes the whole-box numerical lower bound.
These are fixed-margin numerical Weyl bounds, not interval arithmetic
or size/temperature/basis convergence of the material.

`index.json` references the unchanged spectrum archive by relative path
and SHA-256. `methods-results.tar.gz` contains all five reports, logs,
the source scan, mathematical methods and the 32-test log. Archive-only
reproduction is recorded in `replay.json`. From the manuscript root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  validation/development/archive_bismuth_regions.py . . \
  validation/development/bismuth-region48-validation --replay \
  --output .build/region48-fresh-replay.json
```

Use a fresh output filename. Expected unresolved reports are reproduced
as such, not skipped. The snapshots and original raw native calculations
remain in the sibling material bundle; no DFT calculation is repeated.
