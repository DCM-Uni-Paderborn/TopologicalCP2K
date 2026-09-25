# Independent solve probes and symbolic-retry research

CP2K checkpoint: `5218a42573`. This checkpoint changes the deterministic
right-hand sides from `A*x` to vectors independent of `A`, at the same 1e-10
residual threshold. It does not replace the independent spectral-gap test.
The installed Tacho dependency remains the typed-magnitude correction documented
in `../pfaffian-pivot/`; no adaptive pivoting was installed in CP2K.

## Completed Checks

- 600 synthetic size/scale/permutation queries, including 132 singular queries.
  The former probes accepted 72 singular queries. The new probes reject all 132.
  Both native versions resolve the same 450 nonsingular queries correctly and
  leave 18 unresolved. No accepted nonsingular sign changes.
- A unit-triangular dyadic congruence constructs an exact rank-22 matrix of
  order 24 with no zero rows. Six scale/permutation variants are now native
  unit tests. Analytic nullity, not roundoff-sized dense pivots, is the oracle.
- Native units pass with 2 and 4 MPI ranks. Dense and optional sparse directory
  checks pass 45/45 and 29/29, respectively. All three changed CP2K files pass
  `make_pretty.sh --no-cache`. The serial build without Tacho was not rerun.
- Three N=3 stanene queries still agree with dense eigenvalues and the
  independent Hessenberg Pfaffian. The full 33-query material suite has not
  been repeated for this probe-only change.

## Uninstalled Prototype

`tacho_pair_adaptive.cpp` and the private headers implement all-row local
pivots, early weak-pivot exceptions and restart-based symbolic coarsening.
Only problematic pairs are merged with a later coupled group. This is not
true dynamic delayed elimination and is not a distributed factorization.

With independent right-hand sides, the prototype resolves all 468 nonsingular
synthetic queries and rejects the 132 singular queries. The N=6 DZVP/eta=0.5
material query resolves nu=1 after 17 retries, with residual 1.51858e-12 and
gap 0.03623404693181014 hartree. The original dense reference and export hashes
are independently checked through the immutable `stanene-flattening.tar.gz`.
The 10,951,808 scalar panel entries are not the production adapter's condensed
2x2-block count and are not peak process memory. The largest front is still
2950/3744, so this does not establish scalable sparse behavior.

The retained N=12 log used the earlier `A*x`-probe prototype. That exploratory
run was terminated after more than 100 retries without a completed query;
it is not a passing test and provides no index or performance claim.

## Reproduction

The 1,055,892-byte archive has SHA-256
`6d5f299cc104b387d90df23509926e8e1c69859f2748d1d1385f4da6c59d6e42`.
It includes source snapshots, the private header tree, scripts, results and
logs with per-file checksums. Previous archives remain unchanged.

```bash
OPENBLAS_NUM_THREADS=1 python validation/development/replay_pfaffian_probes.py \
  validation/development/pfaffian-probes /tmp/pfaffian-probes-replay.json \
  --library /path/to/libcp2k.dylib
```

NumPy and SciPy are required. Match the archived process-local BLAS environment
for strict numerical replay. Omitting `--library` verifies retained evidence
and the independent dense-reference archive without a CP2K installation.
To rebuild the separate prototype, unpack the archive, configure
`build-serial/pfaffian-stability` against the same Tacho/Boost installation,
and build only `pfaffian_pair_checked`. Other historical target names in its
CMake file do not reproduce the old binaries from current CP2K sources.

No CP2K branch or pull request was published for this checkpoint.
