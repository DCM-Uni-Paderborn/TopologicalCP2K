# Finite Bi neutral-window validation

This record separates three questions: native/reference agreement, stability
with respect to the localizer scale, and physical convergence of the finite
Hamiltonian. Only the first two are established here.

## Retained calculations

The original 300 K, 18-atom calculation remains in
`bismuth-finite-validation/flake3-soc-storage-fixed.tar.gz`.
The new immutable `bismuth-neutral-validation` bundle contains:

- `flake3-dzvp30-spectrum`: the same geometry and DZVP basis at 30 K scalar
  SCF smearing, restarted from 300 K; 168 SCF steps, 2167.64 s.
- `flake3-dzvp30-native`: a seven-step restart at 30 K, with 45 native sparse
  localizer queries and a fresh complete SOC export; 371.28 s.
- `dimer-memory-control`: a small check of the process-tree RSS recorder.
- `methods-results.tar.gz`: exact analysis methods, Gaussian integrator,
  relevant native sources, 16 method tests, neutral/scale scans, interval
  bounds, 300 K archive replay and the native factorization-only comparison.
- `cold-replay.json`: independent replay from the new case/method archives.

Case archives contain the actual per-run runner and analysis methods, inputs,
native outputs, complete spinor states and integrated AO moment cache. The
index records archive SHA256 values; each archive contains a member manifest.
Neither failed earlier diagnostics nor existing evidence archives are replaced.

## Quantitative results

Both 18-atom calculations retain all 234 scalar AOs and 468 spinors, with 90
electrons, UZH DZVP, PBE and post-SCF GTH SOC. The electronic temperature controls
the restricted scalar SCF, not ionic motion or self-consistent SOC occupations.

| Quantity | 300 K | 30 K |
| --- | ---: | ---: |
| Neutral SOC gap (hartree) | 0.000284260227248756 | 0.00008003364849548666 |
| Gap / kBT | 0.29920724 | 0.84241991 |
| Occupied frontier boundary weight | 0.72293164 | 0.72547654 |
| Empty frontier boundary weight | 0.79652992 | 0.80121284 |
| Nontrivial queries / all queries | 8 / 45 | 3 / 45 |
| Minimum covered gap bound (hartree) | 0.00009098526043833004 | 0.00007947715940038735 |

The energy queries are the neutral midpoint plus -1/4, 0, +1/4 of each gap.
Position offsets are 0, 0.75, 1.5 times the atomic x half-width; the five
scales are 0.0001, 0.0003, 0.001, 0.003, 0.01 hartree/bohr. Boundary weights
are Lowdin populations averaged over each complete frontier Kramers doublet.
Ten of 18 atoms have fewer than three neighbors within 3.3 Angstrom. The
weights are invariant under rotations inside each degenerate subspace but
remain dependent on the AO basis and on this boundary definition.

All 45 native 30 K indices match the dense analysis of that same run's export.
The largest distance of a dense gap outside a native inertia bracket is
6.55584e-11 hartree (unchanged tolerance 2e-8 hartree); the maximum native
Pfaffian solve residual is 6.14483e-12. The preceding spectrum and restarted
native calculation have slightly different converged scalar potentials; their
exports are not interchanged in this test.

At the centroid and the neutral midpoint, both spectra have index one on
the whole interval [0.0025, 0.0035] hartree/bohr. An adaptive Lipschitz gap
bound covers the interval with 11 subintervals and 21 midpoint evaluations,
using ||D|| = 19.89763811125872 bohr and a 1e-10 hartree roundoff margin.
This is a floating-point numerical bound, not interval arithmetic. Its
scope is this finite Hamiltonian at this energy and position. It does not
assert scale independence at arbitrarily small or large kappa.

The native 30 K run sampled process-tree RSS once per second. All 365
observations succeeded, with a maximum sum of 20,534,128 KiB (19.58 GiB).
Shared pages can be counted more than once. This is neither an exact memory
peak nor a measurement of the Pfaffian buffers alone. The earlier 168-step
SCF did not use this telemetry. These timings are whole calculations, not
a dense/sparse or one-/two-rank speed comparison.

## Reproduction

From the repository root, with NumPy and SciPy available:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p /tmp/bi-neutral-methods
tar -xzf validation/development/bismuth-neutral-validation/methods-results.tar.gz \
  -C /tmp/bi-neutral-methods
python validation/development/replay_bismuth_neutral_window.py \
  validation/development/bismuth-neutral-validation/flake3-dzvp30-spectrum.tar.gz \
  validation/development/bismuth-neutral-validation/methods-results.tar.gz \
  /tmp/bi-neutral-methods/additional/flake3-dzvp30-frontier.json \
  /tmp/bi-neutral-methods/scripts/scan_bismuth_neutral_window.py \
  /tmp/bi-neutral-replay.json
```

The output path must not already exist. The retained replay verifies 191
archive members, then reproduces all 45 gaps, indices and frontier weights
without numerical change. It does not rerun native CP2K. The earlier 300 K
replay verifies 184 members and likewise has zero discrepancy. The separate
300 K call through the native Tacho C adapter matches 45 indices, with
maximum solve residual 5.40e-12; it tests factorization, not native AO assembly.

For a fresh native calculation the research CP2K build is required. The
retained runner, its CLI arguments, source revision, basis/potential hashes,
executable/library hashes and environment are recorded in each case. The
working paper is not yet a frozen public source distribution.

## Remaining physical controls

The bare flake is unrelaxed and not size- or basis-converged. Both gaps remain
smaller than kBT of the scalar reference, and frontier weights remain boundary
enhanced. Removing degree-one corners, changing the flake size, enlarging the
basis, testing vacuum/grid convergence and further checking occupations are
separate controls, not consequences of the successful native comparison.
