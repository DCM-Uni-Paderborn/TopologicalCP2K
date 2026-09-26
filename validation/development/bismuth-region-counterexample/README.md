# A parameter box that corner sampling cannot certify

This is a deliberate negative test of joint parameter stability, not a
failed CP2K SCF or a failure of the Pfaffian factorization.

For the same retained 30-atom DZVP Hamiltonian, combine kappa = 0.0003 to
0.0025 hartree/bohr with the previous energy half-width (0.02 times the
neutral gap) and both centroid-relative position half-widths (0.05 Angstrom).
All 16 corners and the center have index one. Nevertheless, the interior
neutral-centroid point at kappa = 0.0008 has index zero and a localizer gap
of 0.0008704938737075379 hartree. The two reports share identical native
input and mathematical-method hashes, and the witness is inside the box.

The wider box cannot have a single index. The adaptive method correctly
returns `resolved: false` and `region_z2: null`; eight leaves are unresolved
after its 63-evaluation budget. Its successful sub-boxes must not be read
as coverage of the full requested box. The separate interior-point check
is resolved and returns zero. Complete-spectrum and Pfaffian audits are
retained for both reports.

`methods-results.tar.gz` and `index.json` follow the same format as the
adjacent `bismuth-region-validation` bundle. Replay accepts the deliberate
nonzero exit status only for the explicitly unresolved report, compares
all discrete outcomes exactly and checks all numerical fields.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  validation/development/archive_bismuth_regions.py . . \
  validation/development/bismuth-region-counterexample --replay \
  --output .build/bismuth-counterexample-replayed.json
```

This uses the previously archived complete native snapshot; no new SCF,
native assembly, or material-size convergence is claimed.
