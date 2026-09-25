# Finite Bi basis and size controls

These are isolated, bare, unrelaxed neutral Bi(111) patches with restricted
PBE and post-SCF GTH SOC, not self-consistent spinor DFT or a converged bulk
classification. One control changes the basis at fixed 16-atom geometry;
the other changes the size at fixed DZVP basis. All retain UZH q5 bases,
400/40 Ry grids, 20 Angstrom vacuum and 300 K scalar smearing.

## Spectra and physical subspaces

| Patch | Scalar AOs / spinors | Neutral SOC gap (hartree) | SCF steps | SCF/export time (s) |
| --- | ---: | ---: | ---: | ---: |
| 16 atoms, DZVP | 208 / 416 | 0.008071023745253536 | 131 | 1794.7667 |
| 16 atoms, TZVP | 272 / 544 | 0.008515556694485366 | 119 | 1628.8112 |
| 30 atoms, DZVP | 390 / 780 | 0.0026508955334406947 | 132 | 2601.8473 |

DZVP and the 30-atom case use atomic guesses. TZVP uses the converged
16-atom DZVP restart only as an initial guess and performs a new SCF in the
larger space. All use one MPI rank, two OpenMP threads and one BLAS thread.
Elapsed times are single-run whole-task times, not isolated solver benchmarks.
Sampled process-tree RSS sums for TZVP and 30 atoms peak at 17,195,504 and
20,289,184 KiB, respectively (1588 and 2530 successful one-second samples).
These sums can double-count shared pages and are not exact physical-memory peaks.

`compare_bismuth_basis.py` integrates the rectangular physical AO overlap
between both bases. It compares the occupied spaces and complete frontier
Kramers pairs through overlap singular values. Independent orbital gauge
rotations and nonsingular AO-coordinate changes leave these quantities
invariant. It does not equate independently orthogonalized AO frames.

The DZVP space embeds into TZVP with metric residual 2.85511e-14 and maximum
position residual 1.16355e-12 bohr. The self-consistent Hamiltonian changes
by 0.0053732733 hartree in this common subspace. The occupied projector
distance is 0.10480657, with mean retained weight 0.99686645; occupied/empty
frontier-pair distances are 0.09276768 and 0.06853998. These are real basis
effects, not a basis-convergence certificate. A DZVP self-comparison gives
a Hamiltonian difference of 4.5058e-14 hartree. The new four analytic tests
and the previous 17 method checks pass (21 total).

## Scale dependence and counterexamples

The 45-query scans use each case's own neutral-gap fractions (-1/4, 0, 1/4),
centroid-relative x-half-width fractions (0, 3/4, 3/2), and scales
0.0001, 0.0003, 0.001, 0.003 and 0.01 hartree/bohr. Thus changed gaps and
sizes do not represent the same absolute energy/off-center coordinates.

TZVP has two nontrivial queries instead of three in DZVP. The lower-quarter
energy at the centroid and scale 0.003 becomes trivial. The originally
attempted [0.00275, 0.0035] scale bound fails near 0.002851703 after 271
evaluations. That failed interval is retained, not overwritten. The smaller
[0.003, 0.0035] interval is resolved with index one, 13 evaluations,
seven covered intervals and minimum lower bound 1.1299836574e-4 hartree.

The 30-atom coarse scan has three nontrivial queries, all at scale 0.0003.
A finer, 20-scale centroid/midpoint scan finds a second nontrivial region
missed by the coarse grid. Numerically resolved intervals are:

| Scale interval (hartree/bohr) | Index | Evaluations / covered intervals | Minimum lower bound (hartree) |
| --- | ---: | ---: | ---: |
| [0.00025, 0.0004] | 1 | 17 / 9 | 1.8233565706e-5 |
| [0.0015, 0.0025] | 1 | 15 / 8 | 7.6723922380e-5 |
| [0.003, 0.0035] | 0 | 9 / 5 | 8.7367624335e-5 |

The last interval has the opposite index from both 16-atom bases. It is
numerically resolved, not a factorization failure. Every bound uses the
physical orthonormalized position norm and a 1e-10 hartree margin; it is not
an interval-arithmetic proof or a bulk material conclusion.

The native Tacho C interface independently matches all 20 finer-scan
Pfaffian indices, with maximum solve residual 2.10807e-12. This is a
factorization check of reconstructed matrices, **not native AO assembly**.
The initial invocation's missing-library-path error is retained separately;
the corrected invocation uses the same libcp2k hash as the native runs.

## Native TZVP comparison

A four-step restart evaluates all 45 queries and exports a fresh spectrum.
Native assembly, sparse Pfaffians and MUMPS gap brackets agree with the
independently reintegrated dense reference for every query. The largest
gap distance outside the native bracket is 5.20085e-10 hartree at the unchanged
2e-8 tolerance; the largest Pfaffian solve residual is 1.07709e-13.
The native run takes 486.4152 s; the independent reference takes 9.4900 s.
Its 478 successful RSS observations peak at 21,348,336 KiB, with the same
memory-accounting caveats as above. These timings cover different work.

## Native 30-atom comparison

A five-step restart evaluates the original 45 queries and exports its own
complete SOC frame. All 45 indices agree with the independent reconstruction.
The largest distance outside a native gap bracket is 1.30319e-10 hartree at
the unchanged 2e-8 tolerance, and the largest Pfaffian solve residual is
2.01609e-13. The full native run takes 1169.9872 s; independent reconstruction
and dense analysis take 24.3136 s. There are 1147 successful one-second RSS
observations with maximum process-tree sum 19,076,160 KiB. The restarted
spectral gap differs from the original by 3.61951e-11 hartree.

## Evidence

`bismuth-size-basis-validation/` retains the case archives, original per-case
runners, original analysis methods, complete spinor exports, restarts,
analytic integrals and data/source hashes. The separate methods archive
contains all additional scans, accepted/rejected intervals, cross-basis
diagnostics, method tests and factorization checks. The previous 16-atom
DZVP native comparison is retained in `bismuth-termination-validation/`.

Reconstruction, norm bounds and subspace comparisons are numerical validation;
the basis and size effects above remain physical convergence requirements.

`complete-replay.json` verifies 270 members in six archives. Both 45-query
neutral scans and frontier weights reproduce exactly. The physical cross-basis
report also reproduces exactly from extracted snapshots, Gaussian integrals
and archived methods. All three accepted 30-atom interval bounds and the TZVP
accepted interval reproduce; the rejected TZVP interval again fails after
271 evaluations with seven unresolved subintervals. No failed interval is
interpreted as a positive certificate. These replays do not rerun native SCF.

For example, after extracting `methods-results.tar.gz` to a scratch `methods`
directory, the TZVP scan and accepted interval can be reproduced with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  methods/scripts/replay_bismuth_neutral_window.py \
  validation/development/bismuth-size-basis-validation/flake3-trim-tzvp300-spectrum.tar.gz \
  validation/development/bismuth-size-basis-validation/methods-results.tar.gz \
  methods/additional/flake3-trim-tzvp300-frontier.json \
  methods/scripts/scan_bismuth_neutral_window.py \
  /tmp/bismuth-tzvp-replayed.json \
  --interval-report methods/additional/flake3-trim-tzvp300-stable-interval.json \
  --interval-method methods/scripts/certify_bismuth_scale_window.py
```

The output path must not already exist. For the physical basis comparison,
extract both 16-atom spectrum archives to directories with their case names.
Then run the retained `methods/additional/compare_bismuth_basis.py`, with
`PYTHONPATH=methods/scripts`, arguments `methods`, the two case directories,
and a fresh JSON output path. It independently verifies all retained native
files and reintegrates both bases before comparing subspaces.
