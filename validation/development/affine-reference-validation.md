# Affine-reference and site-resolved EBR validation

Date: 2026-09-24. CP2K source: local commit
`f8961d2c10474882aff067d3719891c8a0f20016`, branch
`feature/spectral-localizer`. The CP2K branch has not been pushed.

## Method

`topology_wyckoff::wyckoff_match_affine` matches a complete supplied affine
catalogue to the generated site families. It checks the integer rank of each
reference direction matrix, the whole tangent identity `B W D = 0`, the affine
offset modulo integers, and the conventional-cell orbit multiplicity. A
complete bijection is required. Output indices and operation witnesses remain
unallocated on invalid, missing, ambiguous or nonbijective input.

The tangent products use bounded 64-bit integer arithmetic; finite coordinates
and offsets use the requested fractional tolerance. The reference must be in
the same cell/origin convention as the group. Labels are supplied metadata,
not an inference from a single sampled coordinate. No reference table or new
Gaussian input keyword is installed at runtime.

The existing native site test accepts `--reference=FILE`. A text fixture starts
with a group count. Each group has a group-number/order/family-count header,
then full operation records (nine rotation integers in Fortran column order,
three translations), and family records (label, multiplicity, three affine
origin coordinates, nine integer direction components in column order).

## Independent references

All external data are pinned to Crystalline.jl revision
`9cfb644a1a4cc1c7baec457b321885459bb98bd1`:

- `test/data/xyzt/generators/sgs/3d/{1..230}.csv` supplies generators.
- `data/wyckpos/3d/{1..230}.csv` supplies named parametrized positions.
- `data/bandreps/3d/{elementary,elementaryTR}/maxpaths/{1..230}.csv`
  supplies the archived Bilbao EBR reference tables.

`prepare_affine_reference.py` accepts only arithmetic AST nodes and x/y/z,
parses coefficients with rational arithmetic, and closes the generators with
integer translations modulo 24. It does not evaluate source expressions or
import third-party code. The source data use translations whose denominators
divide 24; this is asserted, not rounded. Identity is emitted first. Generator
closure and the supplied positions use the same conventional setting.
Each of the 460 generator/coordinate downloads has its immutable source URL
and SHA-256 retained in `affine-reference-sources.jsonl.gz`. EBR table hashes
are retained in the preceding `site-induction-catalogue-comparison.log.gz`.
The reference CSV files are not redistributed in this manuscript repository.

## Results

- All 230 groups and all 1,731 named Wyckoff families match.
- Ordinary/grey and scalar/spinful variants give **6,924** correspondences.
- A separate simultaneous unimodular cell shear, arbitrary origin shift,
  integer coset change and reversal of the operation order preserves the
  complete named specialization graph in every group.
- The EBR dimension multiset agrees **at every individual Wyckoff position**,
  covering **10,398** unreduced maximal-site columns, with zero mismatches.
- Conventional-cell band ranks are divided by the independently counted
  number of pure centering translations before comparison with primitive-cell
  reference band dimensions. No column is filtered by the legacy CSV
  `Decomposable` field.
- The full **24,111** deterministic family/column/convention records agree
  between optimized SSMP, each rank of a two-rank MPI run, and a separately
  instrumented executable. Their common normalized-record SHA-256 is
  `44c09430972ab854506ded33a3af0c95fad7e88734e484cd0d95d2180f427d36`.
- The analytic matcher rejects a wrong tangent even when its origin lies on
  the correct plane, duplicate entries, wrong multiplicity, NaN and enormous
  coordinates. Rescaling the free parameter directions does not alter the
  represented family. These tests also pass under runtime instrumentation.
- The final official regression drivers pass **171/171 in SSMP and MPI**:
  six native unit programs, 129 Gaussian little-group assertions, and 36
  existing topology/export/k-point assertions. No physical reference value
  was changed. The preliminary runs also passed 171/171 in each build.
- `make_pretty.sh --no-cache` checks all four changed CP2K files: 4 checked,
  0 failed. `git diff --check` is clean.

The instrumented files are `topology_wyckoff.F` and the native site/Wyckoff
test programs, compiled with GNU Fortran 16.1.0 using `-O1 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fsanitize=undefined
-fno-sanitize-recover=all`. The rest of CP2K/dependencies is not instrumented.
The instrumented reference executable omits `__SPGLIB` and therefore has no
database interface; it still links the existing CP2K library. This is not a
complete no-SPGLIB CP2K build. Array-temporary diagnostics in the raw debug
log are nonfatal and excluded only from the deterministic-record comparison.

## Reproduction

The scripts use the Python standard library. From a configured source tree:

```text
python prepare_affine_reference.py CACHE reference.dat sources.jsonl
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 topology_site_unittest.ssmp --reference=reference.dat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 mpiexec -n 2 --tag-output topology_site_unittest.psmp --reference=reference.dat
python compare_affine_catalogue.py serial.log CACHE reference.dat
python verify_affine_outputs.py serial.log mpi.log debug.log
```

The EBR CSV cache is prepared with the previously retained
`compare_site_catalogue.py`. The final physical drivers use
`UNIT/topology_(integer|atomic|wyckoff|band|reciprocal|site)_unittest`,
`QS/regtest-little-group`, `QS/regtest-topology`, and `QS/regtest-kp-1`.
As in the preceding validation, Apple ARM runs use one BLAS thread, two
OpenMP threads, a 128 MiB OpenMP stack and the local thread-safe OpenBLAS
preload. The final work directories are
`build-serial/affine-verified-regtests/TEST-2026-09-24_21-42-46` and
`build-mpi/affine-verified-regtests/TEST-2026-09-24_21-42-47`.

Retained evidence includes the preparation/comparison scripts, compressed
final native and driver logs, source hashes, comparison JSON, instrumented
compile/link commands, rejection-test output and formatting log. The numeric
fixture is reproducible from the pinned sources, not a runtime dependency.

## Boundaries and next step

Matching geometry and band dimensions is stronger than aggregate counts,
but does not distinguish different onsite irreps of the same dimension or
prove their complete induced characters. It must not be advertised as
canonical EBR-by-EBR identification. The next step is character-level
correspondence, including double-group factor conventions and antiunitary
corepresentations, followed by attaching validated canonical metadata to
Gaussian-band analysis. A generator for which no composite witness was found
is not declared elementary merely because its dimension matches a table.
Material gap/isolation, complete reciprocal connectivity, transport and
large-system topological validation remain separate obligations.
