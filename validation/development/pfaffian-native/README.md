# Native Delayed-Skew Integration

This record follows, without replacing, the separate prototype in
`../pfaffian-delayed/`. The native CP2K source uses public Tacho symbolic
analysis with its own real-skew numerical panels, delayed 2x2 pivots,
blocked updates and independent solve probes. Source and binary hashes,
build configuration, inputs, outputs and test logs are in `index.json`
and the checksummed `methods-results.tar.gz`.

## Results

- All 33 complete-band material queries are resolved in the first external
  ordering, with unchanged inputs, query energies, gaps and parities.
  Maximum worst-RHS relative residual: `4.0664823872484995e-11`, below
  the unchanged `1e-10` acceptance threshold.
- The previous independently verified archive contains 15 dense references.
  Its maximum dense/sparse gap difference is `4.232725281383409e-15` Ha.
  The material replay invokes native Pfaffian code on the complete-band
  matrices; it does not repeat the large-system SCF calculations.
- Runtime-instrumented native C++: 804 nonsingular synthetic queries have
  the independently expected sign; all 432 singular queries are rejected.
  The prebuilt dependencies are not instrumented.
- Native unit tests pass with two and four MPI ranks and two OpenMP threads.
  The finite/periodic dense, MUMPS and Tacho directories pass 74 checks
  with four MPI ranks; the build without optional sparse libraries passes
  45 dense checks. No regression reference values were changed.

## Diagnosed Integration Failure

The first integration run failed on bismuth SOC. Tolerated one-sided
roundoff gave an asymmetric graph to symbolic analysis, which omitted
an original entry from a front. The reduced 32-site chain reproduces this
failure at 435 of the 930 possible one-sided off-pattern placements.
After consistent real-skew projection and symmetric pattern construction,
all 930 variants give the correct sign. Native unit tests add transposition,
odd permutations, extreme scales and rejection above the skew tolerance.

For physical AII inputs, the discarded part's maximum row/column absolute
sum bounds its spectral norm and is checked relative to the independent
metric-aware gap. This is not permission to restore broken physical time
reversal. Bismuth's unchanged-reference gap is `0.0681101299687` Ha and
the corrected solve residual is `1.39425592196289e-15`.

## Resources and Limits

`PFAFFIAN_MEMORY` defaults to 4096 MiB for requested floating-point
factor, pivot, frontal, Schur-update and solve buffers. It is not a
whole-process cap: sparse input, matching, ordering, integer metadata,
allocator overhead, DBCSR and MUMPS are excluded. Boundary tests cover
an exact sufficient budget, one byte less, collective failure and
inconsistent limits across ranks. The native ordering uses public
defaults, not the prototype's explicit METIS connected-component option.
Factor counts and residuals therefore need not match bitwise.

Pfaffian factors remain on one MPI rank. The material gaps here are
independent near-zero Ritz estimates, not MUMPS inertia brackets.
These runs do not establish asymptotic material convergence or multi-node
scaling; overlapping research timings are not controlled benchmarks.

## Replay

From the paper repository, with Python, NumPy and SciPy available:

```sh
python validation/development/replay_pfaffian_native.py \
  validation/development/pfaffian-native /tmp/native-replay.json
```

The replay verifies every retained file and the preceding prototype's
independent-reference archive before comparing native material results.
Add `--library /path/to/libcp2k.dylib` to recompute the two synthetic
suites through a Tacho-enabled native library. The archived standalone
`native-check/CMakeLists.txt` builds the same C++ source with UBSan and
standard-library assertions. `replay.json` and `replay.log.gz` record a
successful replay including this instrumented synthetic recomputation.
