# Recoverable onsite provenance in Gaussian atomic signatures

Validated on 25 September 2026. Local CP2K commit
`4fced6f5554b5c6535fc22934a05ef5579324fba` retains the onsite data used by
the actual Gaussian atomic-signature path. The CP2K commit is local only;
this record does not imply an upstream release or a new pull request.

## Production change

Every atomic signature column already has a site and a local representation
index. Previously, its complete onsite representation was discarded after
induction. The new compact snapshot retains the exact original operation
mapping, Cartesian polar rotations, site-group products, antiunitary flags,
physical spin lifts, local characters and integer atom-cell shifts. It does
not retain the larger orbit-sized induction cache. Full-group spin lifts
are obtained from the same construction, in original operation order.

The `ATOMIC_ONSITE_*` and `ATOMIC_REFERENCE_SPIN` records are written to the
existing `.little_group` file whenever `ATOMIC_SIGNATURES` is enabled.
There is no new input keyword or runtime character-table dependency.
The additional retained arrays are included in the memory estimate.
`docs/methods/properties/inversion_topology.md` specifies the record format.

The local action fixes the seed exactly. Its factors are the Gamma
double-group signs, not Bloch factors of the original representatives.
The writer checks that these factors are exactly real and equal to +1 or
-1 before storing a product row as signed operation indices. This is
lossless, including antiunitary products. Scalar onsite factors are one;
the spin-1/2 matrices are retained but must not act on scalar characters.

An initial full pair-by-pair text format exceeded the unchanged 2 MiB
regtest output limit for some cubic tests. The compact rows remove this
redundancy. The limit was not relaxed and no huge-output suppression was
added. No physical reference value or symmetry tolerance was changed.

## Independent reconstruction and reference identification

`verify_onsite_export.py` reads only the exported data and independently
enumerates the orbit of each seed with unitary orbit representatives.
Ordinary and grey groups have such representatives. For an operation g
fixing orbit point j modulo L_j, it forms h=t_j^-1 g t_j and evaluates

```text
chi_ind(g,k) = sum_j exp(-2*pi*i*k.L_j) alpha(g,t_j,h) chi_site(h).
U(t_j)^dagger U(g) U(t_j) = alpha(g,t_j,h) U(h).
```

The scalar calculation sets alpha=1. The spinor calculation checks the
intertwiner and its unit modulus explicitly. These sums are compared
entry by entry with the characters obtained from the exported atomic
integer matrix, little-group tables and corepresentation restrictions.
The script also checks complete stabilizer and little-group membership,
Cartesian axes, atom-cell shifts, product signs, column dimensions and
all declared table dimensions. It does not import CP2K's induction code.

The optional `--onsite-reference` audit uses the previously validated
irreptables 3.1.0 fixture and `prepare_wyckoff_characters.py`. It identifies
each complete local character table with a point-group reference by an
explicit orthogonal conjugacy and physical spin-lift conversion. The
source hash, reference group, Cartesian orientation and local-column
correspondence are retained. The reference bound is 2e-4, with characters
scaled by their own positive dimensions; numerical source rounding is
retained. The orbit reconstruction bound is independently 1e-8.

This is a named correspondence in declared point-group axes. It does not
assert that the first geometric conjugacy is Bilbao's canonical onsite
axis convention, equate different EBRs with coincident sampled signatures,
or classify the actual Gaussian occupied space as topologically trivial.

## Results

The final serial and two-rank MPI regression drivers each pass 179/179
checks. The drivers include eight native tests and QS/regtest-kp-1,
QS/regtest-topology and QS/regtest-little-group. Tests use two OpenMP
threads, OPENBLAS_NUM_THREADS=1 and OMP_STACKSIZE=64M.

Per build, the export audit covers:

- 16 Gaussian output files;
- 197 site tables and 697 local representation columns;
- 263 reciprocal queries and 31,537 induced character entries;
- GPW, GAPW, post-SCF GTH SOC, sheared cells, screw phases,
  generic reciprocal points, Wilson points and automatic sampling.

Across the two builds, the largest independent induction residual is
3.5759695934318345e-14. The largest combined reference spin/character
residual is 9.559579355022178e-6. Every onsite reference-column assignment
and its recorded axes agree between builds; the largest axis difference
is zero. These are 16 cases run twice, not 32 distinct materials.

Four deliberately corrupted record types are checked in every file:
onsite character, lattice shift, full-group spin entry and product sign.
All 128 negative checks reject the altered records. The native snapshot
audit separately passes all 2,120 combinations of the 530 Hall settings
with ordinary/grey and scalar/spinful factors. Its existing signature
checks cover 31,079 columns and 707,038 integer matrix entries per run.
Both MPI ranks run the complete audit; this is not MPI speedup evidence.

The focused instrumented build recompiles symmetry, little-group,
corepresentation, character-match, band-representation and atomic-signature
modules with bounds checks, floating-point traps and undefined-behavior
instrumentation. Its analytic test passes. Expected array-temporary
warnings are retained. Other modules are linked from the release library;
this is not a fully instrumented Gaussian/MPI or separate no-SPGLIB build.

## Local stack failure

The high-cutoff `neon-spin-lift.inp` initially crashed before the SCF with
the default OpenMP thread stack. It reproduced in a standalone run with
no concurrent build. The macOS report identifies an access to a thread
stack guard in `general_fill_exp_table` / `grid_cpu_collocate_pgf_product`.
Thus the earlier suspicion of a concurrent relink does not explain this
failure. With OMP_STACKSIZE=64M, the unchanged input succeeds in both full
drivers. The failed logs and relevant crash-report fields are retained.
No global shell configuration or grid-kernel code was changed.

## Evidence and replay

`onsite-export-records.tar.gz` retains the 32 final `.little_group`
files, Gaussian outputs, input directory, both final driver logs,
all-setting logs, focused instrumentation records, earlier failure logs,
the character fixture, environment, code patch and changed source files.
`onsite-export-files.sha256` records the retained file hashes and
`onsite-export-archive-check.json` records the archive verification.
The comparison summary is `onsite-provenance-comparison.json`.
Both export audits, including reference matching and all negative checks,
were repeated after extracting the archive into a fresh scratch directory.
The resulting `onsite-export-replay.json` is byte-identical to the original
comparison summary.

After extraction into a scratch directory, run the scripts from this
development directory against the extracted paths. The default reader
requires NumPy; the optional reference audit also requires SciPy and SymPy.
Exact versions are recorded in the archive's environment JSON.

```sh
tar -xzf onsite-export-records.tar.gz
OPENBLAS_NUM_THREADS=1 python3 verify_onsite_export.py --negative \
  --onsite-reference build-serial/character-convention-reference.dat \
  build-serial/onsite-provenance-final/TEST-2026-09-25_01-08-08/QS/regtest-little-group/*.little_group \
  > serial-export.jsonl
OPENBLAS_NUM_THREADS=1 python3 verify_onsite_export.py --negative \
  --onsite-reference build-serial/character-convention-reference.dat \
  build-mpi/onsite-provenance-verified/TEST-2026-09-25_01-08-59/QS/regtest-little-group/*.little_group \
  > mpi-export.jsonl
python3 compare_onsite_exports.py serial-export.jsonl mpi-export.jsonl
```

The instrumented build recipe is `build_onsite_provenance_instrumented.sh`;
it runs from a matching CP2K source checkout with the serial release
library already built. The final native audit command is
`topology_atomic_unittest.ssmp --all-settings`, or the `.psmp` executable
with `mpiexec --tag-output -n 2`. No canonical onsite/EBR catalogue,
complete connectivity proof or new material convergence result is claimed.
