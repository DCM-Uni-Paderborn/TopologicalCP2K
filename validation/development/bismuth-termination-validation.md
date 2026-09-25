# Finite Bi termination control

This control removes the two singly coordinated corner atoms from the previous
18-atom Bi(111) patch. The simulation cell and all 16 surviving positions stay
unchanged. The resulting patch is still bare and unrelaxed, not passivated.
The neutral electron count changes from 90 to 80. Absolute energies of these
different compositions are not used as a convergence criterion.

## Completed calculations

The immutable `bismuth-termination-validation` bundle retains both calculations,
their per-run methods, inputs, native outputs, complete SOC exports, moment
integrals and source/data hashes. It also contains the successful and rejected
scale-window reports and a separate host/environment record.

| Quantity | Value |
| --- | ---: |
| Atoms / scalar AOs / spinors | 16 / 208 / 416 |
| Scalar reference | restricted PBE, 300 K smearing |
| Basis / grid | UZH DZVP, 400/40 Ry |
| Neutral post-SCF SOC gap | 0.008071023745253536 hartree |
| Gap / kBT | 8.495415604579117 |
| Occupied / empty frontier edge weights | 0.8261641569680339 / 0.7969570762341022 |
| Edge atom fraction | 10/16 = 0.625 |
| Atomic-guess SCF steps / SCF+export time | 131 / 1794.7667 s |
| Native restart steps / full run time | 5 / 301.5002 s |
| Native/reference matching indices | 45/45 |
| Maximum gap distance outside native bracket | 3.473514e-10 hartree |
| Unchanged comparison tolerance | 2e-8 hartree |
| Maximum native Pfaffian solve residual | 2.35828e-13 |

All calculations use one MPI rank, two OpenMP threads and one BLAS thread.
Sampled process-tree RSS maxima are 17,305,040 KiB for the atomic-guess run
and 19,415,776 KiB for the restarted native query/export. The 1751 and 296
one-second observations have no failures. These sums may double-count shared
pages and are not exact memory peaks. The host is an Apple M4 Max with 36 GiB
RAM. The standalone `runtime-host.json` retained in the methods archive was
captured after the calculations, not as a hardware-counter trace.

## Physical interpretation and rejected interval

Removing only two corner atoms changes the finite SOC gap from 7.735 meV
to 219.6 meV at the same scalar smearing temperature. The occupied frontier
pair remains strongly boundary weighted, so the larger gap does not make
this a bulk band-gap calculation. Changing termination is not a size or
basis convergence study and does not establish a magnetic ground state.

The 45-query scan uses each patch's own neutral-gap energy fractions and
x half-width. Three queries are nontrivial, all at the centroid and kappa
0.003 hartree/bohr. Counts alone do not compare identical off-center points
between the two geometries.

The wider interval [0.0025, 0.0035] hartree/bohr is **not** certified for
this termination. The independent midpoint query has index zero at 0.0025
and one at 0.00275. The adaptive cover retains an unresolved neighborhood
of 0.0026215985 after 249 evaluations. This expected rejection is preserved
in `additional/flake3-trim-dzvp300-interval.json` and its log; the certifier
exits nonzero for this record. A positive minimum among the individually
covered intervals must not be mistaken for a resolved whole interval.

The smaller interval [0.00275, 0.0035] is resolved at the centroid and the
neutral midpoint with index one. Its seven covered intervals and 13 midpoint
evaluations use norm bound 15.227198360608146 bohr and a 1e-10 hartree roundoff
margin. The minimum lower bound is 4.645385910375823e-5 hartree. This interval
is contained in both previously resolved 18-atom temperature windows, so it
is a common finite-model interval, not a universal material parameter.

## Reproduction

The generic replay tool checks manifests and operators before evaluating
the neutral scan and the continuous bound. With NumPy and SciPy installed:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p /tmp/bi-termination-methods
tar -xzf validation/development/bismuth-termination-validation/methods-results.tar.gz \
  -C /tmp/bi-termination-methods
python validation/development/replay_bismuth_neutral_window.py \
  validation/development/bismuth-termination-validation/flake3-trim-dzvp300-spectrum.tar.gz \
  validation/development/bismuth-termination-validation/methods-results.tar.gz \
  /tmp/bi-termination-methods/additional/flake3-trim-dzvp300-frontier.json \
  /tmp/bi-termination-methods/scripts/scan_bismuth_neutral_window.py \
  /tmp/bi-termination-replay.json \
  --interval-report /tmp/bi-termination-methods/additional/flake3-trim-dzvp300-stable-interval.json \
  --interval-method /tmp/bi-termination-methods/scripts/certify_bismuth_scale_window.py
```

The output must not already exist. This replays the independently reconstructed
operators, not native CP2K; the fresh native 45-query comparison is separately
retained in the native case archive.
The retained `archive-replay.json` verifies 183 archive members and reproduces
all 45 queries, frontier weights, 13 interval evaluations and seven covered
subintervals without numerical change.
