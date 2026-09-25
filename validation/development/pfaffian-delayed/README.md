# Delayed skew pivots without symbolic restarts

This is a separate serial, skew-specific numerical extension using Tacho
symbolic analysis and BLAS. It is **not installed in the native CP2K path**.
The native tree remains at checkpoint `5218a42573`; no public CP2K push is
implied. The general indefinite solver is unchanged.

## Method

Unstable variables and their full Schur contributions move to the parent
front, including delayed-variable couplings. Original entries are inserted
once. Pivot selection examines both complete updated columns, including the
separator, at relative threshold 0.01 and normalized absolute cutoff 1e-10.
The root uses a rook search. Upper 2x2-block signs and the complete elimination
permutation determine the Pfaffian. Stored factors solve four independent
right-hand sides; the worst relative two-norm residual must be at most 1e-10.
The independent spectral-gap check is retained.

The first unblocked implementation is retained separately. The blocked version
collects 64 accepted columns for a BLAS-3 Schur update. The final public and
checked targets use public `GraphTools_Metis` and `SymbolicTools` interfaces,
without changes to Tacho's private driver. Earlier private read-only accessors
are preserved only as experimental provenance.

## Verified Results

- All **33/33** unchanged stanene queries resolve in the first external
  ordering. The maximum solve residual is 2.05625e-11. All previously resolved
  parities, including trivial small-torus counterexamples, are unchanged.
- All **15 dense comparisons** agree; maximum gap discrepancy is
  4.23273e-15 Ha. Larger cases use independent near-zero Ritz estimates and
  residuals, not rigorous shifted-inertia brackets or dense sign oracles.
- Both synthetic suites pass: **804 nonsingular queries** have the expected
  signs and **432 singular queries** are rejected. Sizes reach 258; tests
  include random permutations, three overall scales, disconnected, separator,
  banded and dense patterns, and analytically singular dyadic congruences.
- **834 component factors** are explicitly reconstructed, with maximum
  relative Frobenius error 1.03665e-13. This is a factor check, not an index
  fallback. The current debug limit is order 512; older binaries used 192.
- The new adapter/extension passes UndefinedBehaviorSanitizer and standard
  library assertions. Dependencies were not rebuilt with instrumentation.
  A separate standalone build reproduces the synthetic classifications.
- A 676-million-entry mock root front is rejected before dense allocation;
  the check's measured process footprint is only 12.2 MB.

Recovered cases relative to the floating-pivot-corrected library scan:

| Case | eta/Delta | Gap (Ha) | Z2 | Worst solve residual |
| --- | ---: | ---: | ---: | ---: |
| DZVP N=6 | 0.50 | 0.036234046932 | 1 | 1.14e-12 |
| DZVP N=12 | 0.50 | 0.016280349816 | 1 | 2.06e-11 |
| DZVP N=12 | 0.75 | 0.291265521469 | 1 | 2.18e-13 |
| DZVP N=18 | 0.75 | 0.319876256213 | 1 | 2.88e-13 |
| TZVP N=12 | 0.75 | 0.293046308528 | 1 | 5.82e-13 |

The complete rerun uses archived energies -0.1592 Ha (DZVP) and -0.16050731 Ha
(TZVP). One earlier standalone TZVP diagnostic uses -0.1592 Ha; it is not
substituted for the basis-specific-energy rerun.

## Memory and Scope

N=12/eta=0.5 has order 14976, 84 fronts, largest front 4582, 46,049,224 scalar
panel entries and 15,804 delayed transfers. Transfers count each ancestor
pass, not distinct variables. N=18/eta=0.75 has order 33696, largest front
5408 and 77,431,958 panel entries. The four-query N=18 process has a measured
peak footprint of 10,524,059,848 bytes, including Python/reference containers.

Panel entries are not condensed 2x2 blocks or whole-process memory. The
prototype's 650-million-entry cap counts front, update and factor arrays but
omits maps, RHS and small work arrays. A configurable workflow-level resource
policy, native integration and distributed factors remain required. No
controlled multi-node scaling or size-converged material classification is
claimed. Single-host research timings can overlap other work.

## Reproduction

The 1,268,493-byte `methods-results.tar.gz` archive has SHA-256
`8a3b851bc3063edee03d16e8eca2e4a970e9263a81fb2812789da5222acde62c`.
It contains per-file hashes, source snapshots, build/test logs and all material
results. Previous reference archives are not overwritten.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  validation/development/replay_pfaffian_delayed.py \
  validation/development/pfaffian-delayed /tmp/delayed-replay.json
```

The replay verifies 791 new archived files, the corrected 53-file baseline,
and unchanged raw exports in `stanene-convergence/` and
`stanene-flattening.tar.gz`. It compares all 33 material results and 15 dense
references. Add `--library /path/to/libpfaffian_delayed_checked.dylib` to
recompute the 1236 synthetic queries. NumPy and SciPy are required. Match the
recorded single-thread/process-local BLAS environment for numerical replay.
`replay-verified.log.gz` and `replay.json` retain the successful replay, including
recomputation. `replay.log.gz` retains an earlier replay-script path error; it
was a wrong archive-member path, not a numerical failure.

To rebuild after unpacking the archive, configure the **standalone** directory:

```bash
cmake -S build-serial/pfaffian-stability/standalone -B delayed-build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/tacho-install
cmake --build delayed-build --parallel 2
delayed-build/check_delayed_memory
```

Tacho (the recorded skew branch), its dependencies and Boost headers are
required. The checked target uses GCC-compatible sanitizer flags. The enclosing
historical CMake file is provenance, not the standalone reproduction entry point.
