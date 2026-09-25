# Stanene: complete-band spectral flattening

Retained on 25 September 2026. This is a separate full-band reference
experiment, not a replacement for the corrected unflattened native AO runs.
The physical motivation is the dispersion-reduction discussion in
Doll, Loring and Schulz-Baldes, Math. Phys. Anal. Geom. 28, 13 (2025),
doi:10.1007/s11040-025-09508-0. The finite-range theorem's sufficient
hypotheses are not claimed for these Gaussian-projected operators.

## Results and limits

- Full spinor exports contain all 52 DZVP or 68 TZVP states at each k-point,
  not just the eight occupied bands. Query energies are -0.1592 Ha and
  -0.16050731 Ha, respectively. The restricted scalar SCF uses the same
  full 8x8x1 MP mesh as the corrected reference; the analysis meshes are
  Gamma-centered. Both cell and Poisson periodicity are XY.
- Exact sign flattening maps shifted energies to +/-1 Ha and preserves
  all occupied projectors. Gaussian neighbor links are not unitarized.
- At eta/Delta=0.75 and 1, DZVP mesh sides 6, 8 and 9 and TZVP side 6 all
  give nu=1. The normalized gaps at eta/Delta=1 are 0.06292476,
  0.02421594, 0.18694608 and 0.06420961, respectively.
- Counterexamples are retained: side 3 at eta/Delta=1 and 1.25 is trivial;
  sides 6 and 8 at eta/Delta=1.25 are trivial while side 9 is nontrivial.
  Nonmonotonic gaps and scale dependence remain. These are not electronic
  band gaps, converged material predictions or a theorem certificate.
- Side 3 also retains t=0,0.5,1 interpolation at three scales. Its
  localizer parity can change even though the sampled electronic gap
  stays open. Neither dropping bands nor rescaling a reported gap can
  substitute for rebuilding and testing the whole localizer.

The manuscript adds the covariant sign formula and a short interpretation;
SI S7.12 / Table S15 contains all 15 exact-sign queries and further
coordinate/metric diagnostics. All 21 queries, including the interpolation
controls, are retained in the machine-readable summary.

## Separate native-kernel checkpoint

CP2K source snapshot: `43c53c9e09` (local research commit, no CP2K push).
`localizer_flatten_dbcsr` computes scale*S*sign(S^-1*(H-E*S)) using existing
Hotelling and Newton-Schulz routines and distributed realification.
It requires a separately validated positive metric and resolved gap;
the sign-square residual alone does not certify the input gap.
This kernel was not used to produce the material table and has no
production input keyword yet. Explicit metric inversion can fill in.

Six constructed complex 12x12 pencils test nonorthogonal metrics,
Kramers-paired spin mixing and spectra of either sign. Maximum entrywise
errors versus the known exact answer are:

| Run | Maximum error |
| --- | ---: |
| Serial, two OpenMP threads | 5.37054e-15 |
| MPI2, row distribution | 4.14600e-15 |
| MPI4, row distribution | 4.51930e-15 |
| MPI4, column distribution | 6.05759e-15 |

Invalid/singular inputs are rejected and H,S are unchanged. The ordinary
test driver passes flattening, localizer and quadratic-DBCSR unit programs
in both serial and MPI2 (3/3 each). `make_pretty.sh --no-cache` passes all
five changed source/build-registration files. Logs and source are archived.

## Archive and replay

- Archive: `stanene-flattening.tar.gz` (70 members plus manifest).
- Size: 40,186,558 bytes.
- SHA256: `4d372410f3eaec7d356078c64243665e6dca45f8e23cb0e6dc01fb37b65b9129`.
- Each member has its own size and SHA256 in `manifest.json`.
- Includes complete CP2K inputs/outputs, full .topology/.mmn/.eig exports,
  UZH basis and SOC potential files, CMake caches, analysis scripts,
  source snapshots and native-kernel logs.

From the paper repository, with NumPy and SciPy installed:

```sh
OPENBLAS_NUM_THREADS=1 python validation/development/replay_stanene_flattening.py \
  validation/development/stanene-flattening.tar.gz /tmp/stanene-flattening.json \
  --recompute mesh3
```

All member hashes are verified. `--recompute all` recalculates every scan;
individual labels are `mesh3`, `mesh6`, `mesh8`, `mesh9`, `mesh6-tzvp`.
Large meshes require dense matrices and appreciable RAM/time. The retained
validation reruns mesh3 from the archive and checks its Pfaffian signs and
gaps; larger cases retain their original complete analyses and hashes.
The independent Pfaffian path uses orthogonal Hessenberg reduction and
checks 32 small random matrices against a recursive reference.

No previous archive is overwritten. Material inputs were run with the
same full-band export implementation as the source snapshot; runtime
commit banners can precede incremental rebuilds.
