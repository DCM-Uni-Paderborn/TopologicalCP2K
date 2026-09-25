# Floating-point pivot correction and unresolved material controls

This is a separate checkpoint, not a replacement for `../stanene-convergence/`.
The CP2K checkpoint is `b5c9c64ff6`; the Tacho skew extension is revision
`02ef047cb4e1e16959714ce50462bb09511d0e9a` plus the archived floating-pivot patch.
The production change uses typed magnitudes rather than integer `abs` in the
dependency. It leaves CP2K matching, local candidate rows, the 0.01 pivot
acceptance ratio and the 1e-10 solve-residual threshold unchanged.

## Evidence

- The instrumented original kernel reports `abs(-0.082575...) == 0` with a
  four-byte integer result, while `ArithTraits<double>::abs` retains the magnitude.
- A frozen 8x8 integer matrix has exact Pfaffian 60 by recursive expansion.
  The original dependency fails six scale/permutation variants; the corrected
  dependency resolves all six, with largest residual 1.86e-15.
- The native sparse unit program passes on two and four MPI ranks. The MPI
  localizer directory passes 45/45 dense and 29/29 MUMPS/Tacho checks. Formatting
  passes all three changed CP2K files. Serial builds do not link Tacho and were
  not rerun for this dependency-only change.
- All 33 original material queries were repeated with the rebuilt CP2K library
  and the process-local thread-safe BLAS. All spectral gaps are unchanged.
  All 27 cases resolved both before and after retain their index. N=3 DZVP,
  eta/Delta=0.5 changes from unresolved to nu=1, matching the independent dense
  oracle. Four previous failures remain, and one TZVP query becomes unresolved.

The five unresolved queries in this runtime are DZVP N=6 and N=12 at eta=0.5,
DZVP N=12 and N=18 at eta=0.75, and TZVP N=12 at eta=0.75. No reference integer
is substituted for a failed factorization. The first batch stopped at an
assertion requiring every previously resolved query to remain resolved. Its
complete failure log is retained. A separate resume reuses only completed,
hash-checked reports and finishes the missing cases while explicitly recording
the new unresolved outcome.

For the new TZVP case, the original dependency with the local BLAS resolves
nu=1 (8.02e-12 residual), as does the patched dependency with Homebrew BLAS
(7.81e-11). The patched dependency/local-BLAS combination fails all three
orderings, with smallest residual 1.03405e-10. Thus this correction is necessary
but does not establish platform-independent numerical robustness.

## Exploratory Policy, Not Adopted

An isolated all-row maximum-column pivot variant resolves that TZVP case with
residual 2.25153e-13. It still fails the tested DZVP N=6/eta=0.5 and N=12/eta=0.75
queries, and has not passed the full matrix of tests. The private experimental
header and candidate-matching diagnostic adapter in this archive are not the
production sources. The production dependency header is separately retained in
`dependency/`, and the patch is under `tools/toolchain/scripts/stage5/`.

Near-singular local fronts are also observed with strong initial matchings and
full local-column pivot searches. Delayed pivots between supernodes, accurate
factor solves and bounded fill/memory remain the next solver work. No pivot
perturbation, element dropping or relaxed solve threshold was introduced.

## Replay

The archive contains 53 files plus a per-file manifest, totals 97,410 bytes,
and has SHA-256 `6b74e1d0f2fc82b3e7f9a5518d0364b2ce3d111eef3a359afb64f882be9293a3`.
It reuses the unchanged earlier raw-export archives. Runs overlap and are not
performance measurements. Library hashes distinguish rebuilt binaries; runtime
commit banners alone are not used as provenance.

```bash
python validation/development/replay_pfaffian_pivot.py \
  validation/development/pfaffian-pivot /tmp/pfaffian-pivot-replay.json
```

With NumPy/SciPy and the corrected Tacho-enabled shared library:

```bash
OPENBLAS_NUM_THREADS=1 python validation/development/replay_pfaffian_pivot.py \
  validation/development/pfaffian-pivot /tmp/pfaffian-pivot-recomputed.json \
  --library /path/to/libcp2k.dylib --recompute mesh3
```

The retained numerical replay reruns all three N=3 queries and compares dense
eigenvalues and the independent Hessenberg Pfaffian. Match the archived BLAS
environment for reproducing larger near-threshold success/unresolved outcomes.
