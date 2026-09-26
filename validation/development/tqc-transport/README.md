# Little-group MPI transport validation

Native checkpoint: `ef37b012c7`; baseline: `7df186991f`.
This is a local arithmetic/transport validation, not a scaling benchmark or a material survey.

- Final `QS/regtest-little-group`: 139/139 checks in SSMP, three-rank PSMP and four-rank PSMP.
- Shared `QS/regtest-property-wilson`: 47/47 checks in SSMP and four-rank PSMP.
- Baseline little-group directory: 135/135 checks per build before the new unitary-only test.
- All runs use two OpenMP threads and one BLAS thread; native calculations run sequentially.
- Full ordered TQC records agree in 66 comparisons (451,330 records), with exact integer tokens.
  Floating tolerance is fixed at 1e-10; the largest difference is 2.49e-14.
- Nine production-graph controls distinguish valid data, incompatible counts and invalid factors.
  The invalid second edge belongs to rank one in the three-/four-rank tests and aborts collectively.
- All 530 Hall settings pass in the optimized mathematical driver and with `topology_symmetry.F`
  plus its driver instrumented. Instrumentation does not cover the whole CP2K library.

The archive retains 494 payload files: native input/output records, final modified sources,
build/formatting logs, comparison reports and methods, and the synthetic collective driver.
Binary/shared-library hashes are in `index.json`; preceding hashes are in the archive.
No restart wavefunctions, build objects, credentials or Git metadata are included.

From the manuscript root, replay without running CP2K:

```sh
python validation/development/archive_tqc_transport.py . validation/development/tqc-transport --replay
```

To reproduce the native directory runs from the matching CP2K source checkout, use the standard
`tests/do_regtest.py` driver with `--maxtasks 1 --timeout 180 --skip_unittests --ompthreads 2`,
the exact directory as `--restrictdir`, and `--mpiranks 1`, `3` or `4` as above.
Build-specific driver logs retain the full settings and executable directory. `make_pretty.sh`
passes with zero failed checks. Mathematical sweeps use
`topology_symmetry_unittest.psmp --all-settings` and the focused debug executable.

The archived `build-mpi/tqc-transport-probes/edge_collectives.f90` links against the matching
`libcp2k` using `mpifort`, `-Ibuild-mpi/src/mod_files`, `-Lbuild-mpi/src -lcp2k.2026.2` and the
corresponding runtime library path. It is compiled with `-O1 -g -fcheck=all`
`-ffpe-trap=invalid,zero,overflow -fopenmp`. Run each of `valid`, `incompatible`, `invalid`
as its first argument on one, three and four ranks, in separate working directories.
The first two return zero with respectively zero and one incompatible segment; the last returns
nonzero with `Little-group compatibility algebra failed`. Each is bounded by a 40-second timeout.

The focused mathematical debug build compiles `src/topology_symmetry.F` and
`src/topology_symmetry_unittest.F` with `-D__SPGLIB -O1 -g -fcheck=all`
`-ffpe-trap=invalid,zero,overflow -fsanitize=undefined -fno-sanitize-recover=all`, using its module
directory before the serial build's module directory and linking CP2K, SPGLIB and BLAS.
Both optimized/debug logs report 12,720 little-group cases and 12,720 segments,
6,360 projective character tables and 5,324 corepresentation tables. These are algebraic checks.

The AO metric, selected frames, small representation tables, character construction and integer
lattice calculations still remain replicated. Neither large SOC material convergence nor
multi-node speedup is established here.
