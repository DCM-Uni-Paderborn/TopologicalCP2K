# Larger complete-band stanene localizer controls

This archive retains the original dependency build and its unresolved outcomes.
The subsequent floating-pivot correction and separate numerical rechecks are in
`../pfaffian-pivot/`; they do not replace the reports or hashes recorded here.
Use `replay_pfaffian_pivot.py` for the corrected-library rechecks. The original
replay intentionally requires the original success/unresolved classification.

These records distinguish three changes: native flattening in complete primitive
Bloch blocks, a quadratic-cost native AII basis conversion, and an independent
sparse full-band reference used to explore larger analysis tori. The last of
these is not a new sparse production representation in CP2K.

## Physical and numerical controls

The earlier stanene geometry, PBE, GTH-q4 SOC, full-band exports and unmodified
Gaussian coordinate links are retained. Unless noted otherwise, the scalar SCF
mesh is 8x8x1, with UZH DZVP/200 Ry. The analysis tori are independent of this
SCF integration mesh. Flattening uses all bands with energies +/-1 Ha.
The absolute query energy is -0.1592 Ha for DZVP and -0.16050731 Ha for TZVP.

| Torus/basis | Gap at eta/Delta=1 (Ha) | Z2 |
| --- | ---: | ---: |
| 12x12, DZVP | 0.229913902646 | 1 |
| 15x15, DZVP | 0.252124305564 | 1 |
| 18x18, DZVP | 0.268152343435 | 1 |
| 12x12, TZVP/400 Ry | 0.231852482381 | 1 |
| 12x12, DZVP, 12x12 SCF mesh | 0.229840135569 | 1 |

All three larger DZVP sizes also give Z2=1 at eta/Delta=1.25 and 1.5.
These are sampled scales, not proof of an uninterrupted plateau. The gaps
still change with size, and these numbers are not electronic band gaps.
The finer SCF mesh changes the sampled indirect electronic gap from
0.002708107 to 0.002707169 Ha and the eta=1 localizer gap by 7.38e-5 Ha.

The sparse reference reproduces 15 dense reference gaps at N=3,6,8,9 DZVP
and N=6 TZVP within 4.3e-15 Ha. Its complete MO frames are expressed in an
orthonormal AO basis, then k/-k cosine/sine pairs retain physical time reversal.
No coordinate link is replaced by its polar factor, and no band is discarded.
Sparse shift-invert eigenvalues are checked against the original skew matrix;
these are not rigorous inertia-based gap brackets.

Five sparse Pfaffians remain unresolved: eta=0.5 at N=3,6,12, and eta=0.75
at N=12,18. Their nonzero spectral-gap estimates do not rescue an inaccurate
factorization. Dense N=3 and N=6 references separately determine those indices;
they are not substituted for failed sparse values. The fixed solve-probe
tolerance stays 1e-10. At most two additional exact permutations are attempted,
with the corresponding Pfaffian permutation parity included. Every attempt,
including unsuccessful ones, is recorded. A resolved integer is never obtained
by relaxing the tolerance or dropping small matrix elements.

## Native code and runtime evidence

The complete-Bloch native flattening agrees with exact constructed pencils,
including degeneracy, nonorthogonal metrics and invalid-input rejection.
The native AII conversion uses
`Q^dagger L Q = (L + CLC + i(CL-LC))/2` for a signed permutation C, instead
of multiplying by a dense Q. Thirteen direct matrix comparisons reproduce the
original dense formula within 1.4e-15, including broken-structure rejection.
The serial and MPI localizer directories pass 45/45 and 74/74 assertions.

One initial native N=8 run crashed in the loaded Homebrew complex BLAS kernel.
The retained crash excerpt localizes the failure to complex Cholesky inside
the generalized eigensolver; it is not a fitted-reference discrepancy.
A standalone identity-matrix Cholesky probe passes with both BLAS builds,
so this is not evidence that every Homebrew Cholesky call is broken.
The native rerun selects the local ARMv8 BLAS consistently through its process
environment. No shell startup file or system-wide library installation changes.
The first local-BLAS rerun was intentionally terminated after profiling the
cubic AII transform; its replacement is retained as `native/mesh8-final`.
The interrupted run is not counted as successful validation.
The completed two-rank N=8 run gives gaps 0.0591899670928424,
0.0242159444548418 and 0.0410187646703921 Ha at eta=0.75,1,1.25,
with indices 1,1,0. The largest independent-reference difference is
5.8e-14 Ha. This retains the small-torus counterexample. Binary hashes
and the process-only BLAS environment are recorded separately: the native
run predates the additional nonfinite-output guard, while the final rebuilt
unit/regression tests include it. Incremental runtime banners and the Git
revision recorded by the launcher at completion are not binary provenance.

## Archives and replay

`index.json` records each archive hash and size, the CP2K source snapshot and
the Tacho dependency revision. New complete exports are split by case to avoid
oversized repository objects. `methods-results.tar.gz` contains method sources,
scripts, all sparse reports, native logs, runtime diagnostics and unit/regression
logs. Earlier exports are reused from the unchanged `../stanene-flattening.tar.gz`.
No earlier archive or counterexample is overwritten.
The retained replay validates 195 archived members across ten cases.
The additional numerical replay repeats all seven sparse N=3 and N=6 queries,
including both unresolved eta=0.5 factorizations. All gaps and resolved indices
are reproduced. `replay.json` and `recompute.json` retain the checks and the
rebuilt library hash; `recompute.log` retains the numerical rerun output.

Verify all members and numerical comparisons without running a new calculation:

```bash
python validation/development/replay_stanene_convergence.py \
  validation/development/stanene-convergence /tmp/stanene-convergence-replay.json
```

With NumPy/SciPy and a Tacho-enabled CP2K shared library, append for example:

```bash
  --library /path/to/libcp2k.dylib --recompute mesh3 mesh6
```

The replay rejects any changed input, missing report, unresolved-to-resolved
substitution or contradictory index. `--recompute all` repeats the larger
reference studies and needs substantially more memory and time. Elapsed times
in these overlapping single-machine runs are not scaling benchmarks.

This evidence improves finite-volume and basis/SCF controls. It does not prove
asymptotic material convergence, finite-range theorem hypotheses, a large-flake
classification, robust Pfaffian factorization for every query, or multi-node
scalability.
