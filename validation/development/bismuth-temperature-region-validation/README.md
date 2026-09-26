# Rechecked 48-atom localizer boxes at 100 K

These are fresh numerical gap bounds for the converged 100 K full-AO Bi
Hamiltonian, not a transfer of the earlier 300 K bounds. The latter are
smaller than the measured physical Hamiltonian perturbation norm.

The energy range is the new neutral midpoint plus/minus 2% of the new
neutral electronic gap. Both in-plane positions vary by plus/minus 0.05
Angstrom around the unchanged centroid. The numerical margin remains
1e-10 Ha; these are floating-point Weyl bounds, not interval arithmetic.

| Kappa range (Ha/bohr) | Z2 | Gap lower bound (Ha) | Evaluations | Covered leaves |
| --- | ---: | ---: | ---: | ---: |
| 0.00025 to 0.0003 | 1 | 1.43644681e-6 | 13 | 7 |
| 0.0015 to 0.002 | 1 | 5.02531520e-5 | 13 | 7 |
| 0.003 to 0.0035 | 0 | 7.67926875e-4 | 5 | 3 |

No unresolved leaves remain. All 51 independent full-spectrum anchor/corner
checks agree with their box index and the perturbation bound. Positive
coverage at 100 and 300 K does not establish an unbroken interpolation
between temperatures or temperature/material convergence.

The index references the unmodified complete 100 K spectrum archive in
`../bismuth-temperature48-validation/`. The methods archive contains the
frozen mathematical driver, source scan, all three reports and logs, and
the retained method-test log. No CP2K calculation is rerun by the replay.
From `validation/development`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python archive_bismuth_regions.py \
  . . bismuth-temperature-region-validation --replay \
  --output ../../.build/temperature-regions-replayed.json
```

The first two positional arguments are unused in replay mode. The output
path must be fresh. `replay.json` records the full-partition reproduction,
including numerical fields, provenance, anchors and every covered leaf.
