# Real-space site induction and composite band witnesses

Validated 24 September 2026. CP2K source checkpoint:
`a380ac91c17f7cb5e0ba181962c6c1da929862c5`, local branch
`feature/spectral-localizer`. This checkpoint was not pushed to CP2K or its fork.

## What changed

`topology_site_relations` constructs real-space specialization witnesses from
the generated connected Wyckoff closures and independently constructed onsite
irreps/corepresentations. It solves an integer endpoint lift, checks every
source stabilizer on the entire displacement and performs finite site-group
induction at the destination. Destination orbit images retain the original
seed's local column labels. The source and target onsite decompositions include
physical spin factors and the nonunit character norms of antiunitary
corepresentations.

Singleton inductions define equivalence classes. Direct sums provide composite
witnesses for those classes, including exceptions at maximal site-symmetry
groups. Recursive expansion terminates by strictly decreasing positive band
dimension. The production atomic catalogue retains its original columns and
independently verifies every induction matrix and the complete expansion
against its reciprocal signature matrix, together with band dimensions.

The Gaussian detail output now includes `SITE_INDUCTION_SIZES`,
`SITE_INDUCTION_EDGE`, `SITE_INDUCTION_ENTRY`, `SITE_INDUCTION_COLUMN` and
`SITE_INDUCTION_EXPANSION`. It supplies geometric endpoints, branching,
equivalent representatives and decomposition provenance separately from
sampled integer atomic membership. No physical reference value was changed.

## Verification

- All 530 Hall settings, ordinary/grey groups and scalar/spinful factors:
  2,120 variants, 22,592 site-specialization relations and 180,736 independent
  reciprocal character identities. Recursive expansions are also checked
  against whole-band characters at all eight sample points.
- Identical 4,774 result records in optimized SSMP, each MPI rank and the
  instrumented executable. Instrumentation recompiles the integer-lattice,
  band-induction, site-induction and unit-test sources with bounds checking,
  floating-point traps and undefined-behavior instrumentation; it does not
  instrument all of CP2K or the external libraries.
- Analytic inversion/mirror paths, intermediate positions, a nonorthogonal
  cell, shifted origin, integer coset changes, magnetic half translation,
  transitive direct sums and invalid witnesses.
- The ordinary spinful Fm-3m exception is detected. Spinful grey maximal-site
  representations are never classified composite in the complete sweep.
  The regular CI subset additionally checks the two spinless grey S4-site
  exceptions in space group 84.
- A separately compiled test without `__SPGLIB` passes the analytic cases;
  this is not a complete CP2K build without SPGLIB.
- Official regtest drivers: **171/171 in SSMP and 171/171 in two-rank MPI**.
  These comprise six native programs, 129 Gaussian little-group assertions
  and 36 existing topology-export/k-point assertions.
- Independent exact SymPy readers check 15 Gaussian site-induction files
  per build: nonnegative coefficients, `A E = A`, dimensions, lifted endpoint
  geometry, source stabilizers along the path, equivalence components and
  recursive witness provenance. The existing star-quotient oracle still
  validates all eleven exported compatibility quotients per build.
- `make_pretty.sh --no-cache` passes all eight changed files; `git diff --check`
  passes. The first unconstrained formatting sweep was stopped after touching
  unrelated baseline documentation/formatter issues; unrelated edits were
  undone. The retained cold targeted result is authoritative for this change.

An initial regtest attempt used three stale native test executables after
changing the shared Fortran derived-type/interface layout. Those executables
aborted; all Gaussian assertions passed. Rebuilding all six native targets,
then repeating the complete drivers, resolved the failures without source or
reference changes. Both initial logs are retained to document that distinction.

## Independent EBR references

`compare_site_exceptions.py` compares the *space-group coverage* of composite
maximal-site representations with Tables II-IV of Cano et al., Physical Review
B 97, 035139 (2018), https://arxiv.org/abs/1709.01935. All 2,120 variants agree.
This is deliberately not called a canonical irrep-by-irrep comparison.

`compare_site_catalogue.py` independently counts the scalar/spinful columns in
the Bilbao EBR tables archived by Crystalline.jl at revision
`9cfb644a1a4cc1c7baec457b321885459bb98bd1`. It checks both `elementary` and
`elementaryTR`, all 230 space groups and every native Hall-setting result.
All **2,120/2,120** per-setting counts match. The JSONL report records the
immutable source URL and SHA-256 of every reference CSV. Tables are used only
for validation, not copied into CP2K's runtime code or distributed here.

Counting each space group once yields 3,383 scalar / 2,258 spinful generators
without time reversal, and 3,141 scalar / 1,616 spinful generators with it.
The 2,258 value also occurs in the archived data and differs from the aggregate
2,263 printed in the 2018 paper. No native result was changed to target either
total. A discrepancy between historical totals is not a license to relabel
individual representations without a direct character/equivalence comparison.

## Reproduction and retained evidence

Build `cp2k-bin` and all six native test targets before running the drivers;
updating only the shared library is insufficient after an interface change.
The relevant target names are `topology_integer_unittest`,
`topology_band_unittest`, `topology_wyckoff_unittest`,
`topology_atomic_unittest`, `topology_reciprocal_unittest` and
`topology_site_unittest`.

Use `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` for algebra sweeps. For the
Gaussian drivers, use two OpenMP threads and `OMP_STACKSIZE=128M`; MPI uses
two ranks. The local GNU Fortran 16.1/OpenBLAS build uses the existing
thread-safe OpenBLAS preload documented in the preceding validation record.
The external operation database is SPGLIB 2.7.0.

```sh
topology_site_unittest.ssmp --all-settings
mpiexec -n 2 --tag-output topology_site_unittest.psmp --all-settings
python3 verify_site_sweeps.py serial.log debug.log mpi.log
python3 compare_site_exceptions.py serial.log /path/to/libsymspg.dylib
python3 compare_site_catalogue.py exception-comparison.jsonl reference-cache
python3 verify_site_induction.py /path/to/regtest-little-group/*.little_group
python3 verify_star_quotients.py /path/to/regtest-little-group/*.little_group
```

The two character/count comparison scripts use Python's standard library;
the exact Gaussian certificate check uses SymPy. None is a CP2K dependency.
The retained logs and Gaussian certificate archives use the `site-induction-`
prefix. The scripts exit nonzero on missing coverage or inconsistent results.

## Interpretation

The returned decompositions have real-space induction/homotopy witnesses, not
merely equal character vectors at selected momenta. The reference comparisons
test exception coverage and generator *counts*. They do not identify each
remaining generator with a canonical Wyckoff/irrep label, prove completeness
of an arbitrary supplied reciprocal graph, establish a global spectral gap,
or classify a real material from sampled atomic membership alone.
