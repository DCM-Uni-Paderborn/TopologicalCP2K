# Bloch Kubo validation, 2026-09-24

Local research implementation; no push or PR publication.

## Numerical evidence

- Full validation before dimensional additions: serial 78/78, MPI 108/108.
- Kubo dimensional validation: serial and MPI 24/24 (1D/3D initially execution-only;
  numeric references subsequently added and rerun).
- MPI uses four ranks and two OpenMP threads; serial uses two OpenMP threads.
- Ne/SOC/PBE at 2000 Ry, primitive 3x1x2 versus explicit-Gamma 3x1x2 supercell:
  both sigma_iso = 6.21105320e-9 S to all printed digits.
- Total energies: primitive -34.923548219306689 Ha;
  six-cell supercell -209.541289276241002 Ha.
- At regression cutoff 200 Ry, PADE, primitive and supercell potentials differ:
  sigma_iso = 6.42389197e-9 and 6.35401974e-9 S, respectively.
  These references are intentionally distinct, not claimed cutoff-converged.
- Bi2/SOC: serial and MPI sigma_iso = 3.83512768e-6 S.
- 1D Ne/SOC: both sigma_iso = 3.85441018e-18 S*m.
- 3D Ne/SOC: both sigma_iso = 1.07057468e-1 S/cm.

## Local runtime and performance issues

All final runs set OPENBLAS_NUM_THREADS=1 and explicitly preload
build-serial/openblas-thread-safe/lib/libopenblas.dylib using DYLD_INSERT_LIBRARIES.
Without the preload, mixed local BLAS dependencies caused an early CrSBr SCF
ZHEGVD failure with four MPI ranks and two OpenMP threads, before the modified
property routines. The identical input passes with the preload.

The 2000-Ry implicit-Gamma supercell comparison was deliberately terminated after
profiling identified the quadratic Miller-index search in pw_copy_match. It used
about 15 GB RAM. The completed comparison instead specifies an explicit Gamma
KPOINTS section, avoiding that conversion. The underlying lookup bottleneck has
not yet been fixed by the Bloch Kubo change.

## Scope still outstanding

The independent Bloch Kubo mesh now supports validated K290/SPGLIB eigenframe
reconstruction, including scalar/SOC AO actions and cell-image/antiunitary phases.
Wilson/TRIM links and localizer property orbits still need their own production
symmetry paths; they are not supplied by the reduced SCF input.
General 230-space-group TQC requires irreps, compatibility and EBR data, not only
the existing inversion classification. Sparse scaling and broader converged
material tests remain open. No claim of full-goal completion is made.

## Property symmetry validation

- Final directory drivers: serial 132/132; four MPI ranks x two OpenMP threads
  162/162. Logs: build-{serial,mpi}/property-final-validation2.log.
- Coverage includes QS/regtest-kp-1, moments, topology, Kubo, localizers and
  quadratic pseudospectra; MPI additionally covers MUMPS and Tacho directories.
- Mathematical units pass in serial and four-rank launches: Bloch/projected
  currents, topology symmetry/Wilson/curvature/snapshot, dense localizer,
  quadratic pseudospectrum, and MPI sparse localizer.
- Ne/SOC 2x2x2: 4/8 representatives, sigma_iso = 0.107057468 S/cm.
- Bi2/SOC 3x3x3: 6/27 representatives, sigma_iso = 303.572782 S/cm.
  The full nine-component tensor differs from full-grid by at most 2.71e-15
  relative to its largest entry (printed data); K290 uses 12 alternative
  operations, all checked against target H, S and current.
- Primitive diamond/SOC 2x2x2: 3/8 representatives for both backends,
  sigma_iso = 996.335574 S/cm, equal to full-grid at printed precision.
  The full tensor deviation is 2.21e-11 relative to its largest entry.
  These finite-cutoff tests are covariance checks, not converged material data.
- The half-translation unit verifies that two actions produce the full-cell
  Bloch phase, that reciprocal shifts preserve the cell gauge, and that a
  nonbijective atom mapping is rejected.
- Both final builds completed without warnings. make_pretty.sh: 42 files,
  37 cached, five checked, zero failures; git diff --check passed.

All final drivers also set OMP_STACKSIZE=128M. The primitive-diamond MPI failure
without it was an OpenMP worker stack overflow in grid_cpu_collint.h's local
exponential tables, before the property calculation; serial timed out there.
The same calculations complete with the larger stack. No grid source or numerical
reference was changed for this runtime configuration issue. Unused restart
output is disabled in the Bloch comparison inputs, removing the driver's
2-MiB output-limit failures without adding suppressions.
