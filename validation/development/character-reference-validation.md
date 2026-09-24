# Character-level reference correspondence

Validated on 24 September 2026 in the local spectral-localizer research tree.
The CP2K change is committed locally as `590d9ee878`; it has not been pushed
or presented as an upstream release. The manuscript repository remains private.

## Results

- 230 ordinary space groups, both scalar and spinful character tables.
- 2,700 tabulated reciprocal-point tables and 8,907 irreducible columns.
- Every complete column table has a unique reference correspondence.
- Largest external residual: 1.2092068868892655e-5.
- Native geometry/factor/character tolerance: 1e-9. Printed external tables
  and reference spin matrices: 2e-4, without rounding their characters.
- 11,607 normalized case/label records agree in SSMP, separately on each
  of two MPI processes, and in the focused instrumented executable.
- Normalized-record SHA-256:
  `5032f7500f36f119fe9d9edb088b493eca7a771f691e984fef945ff0253a9a87`.
- Official local drivers: SSMP 173/173 and MPI 173/173, including eight
  native unit programs, 129 Gaussian little-group assertions and 36 existing
  topology/export/k-point assertions. No physical references were changed.
- `./make_pretty.sh --no-cache` passed on all five changed CP2K files.

These are character/algebra and existing Gaussian integration tests, not
230 converged DFT material calculations or a multi-node scaling study.
The two-process reference run repeats the complete audit on each process;
it is not a distributed character-generation benchmark.

## Independent source and conventions

The fixtures read only the published data files from `irreptables==3.1.0`
(GPLv3), part of IrRep. No IrRep implementation is imported and no table
is linked into CP2K. The 460 source-file hashes are retained separately.
The source is described by Iraola et al., Comput. Phys. Commun. 272,
108226 (2022), DOI 10.1016/j.cpc.2021.108226, and
<https://pypi.org/project/irreptables/3.1.0/>.

The previously validated affine fixture supplies conventional centering
translations. Its primitive basis P is computed by an integer Hermite
normal form. Operations become P^-1 W P and P^-1 tau. The fixture uses
k_native = -P^T k_table for CP2K's exp(-2*pi*i*k*T) translation convention.
The global sign is declared before any table matching. No column-dependent
phase fitting, label permutation to conceal errors, or selective character
conjugation is performed.

An independent central-projector audit rejects 25 boundary-point tables
when the sign conversion is omitted, and accepts all 2,700 when it is
included. For the spinful I212121 W example, the supplied one-dimensional
characters have chi(z)*chi(y)/chi(x)=+i, whereas the unconverted native
factor is -i. This diagnoses a convention mismatch, not a claim that IrRep
is wrong. The raw and convention-converted audit logs are retained.

Printed coordinates are reconstructed with denominator at most 48 and
absolute changes below 5.1e-6, accounting for five-decimal coordinate data.
Character values are never rationalized. The reference Cartesian spin
frame is found by an orthogonal intertwiner. Polar-unitarized reference
spin matrices are used only to determine that frame; all actual comparisons
use the original finite-precision matrices.

An independent recomputation of a pi-rotation spin lift after a different
floating-point cell normalization changed its SU(2) representative in the
I23 reference case. `little_group_factor_system` now optionally returns
the exact spin lifts used to generate its factors, including the physical
Theta matrix in antiunitary operations. The matcher uses those lifts,
instead of guessing their signs from a second Cartesian reconstruction.

## Native checks

`topology_character_match.F` validates an explicit operation isomorphism,
unit-modulus gauge phases, and the twisted factor-coboundary identity.
It checks dimensions, orthogonality, Wigner norms, completeness, and
Hermitian central idempotents in the unitary twisted group algebra.
The supplied character columns must then match bijectively. A
one-dimensional character twist can preserve factors but change names;
the API therefore never infers canonical gauge phases from factors alone.

The standard unit test requires no reference data or SPGLIB database.
It covers ordinary/projective C3, complex rephasing, label-changing twists,
paired/doubled grey corepresentations, and reordered antiunitary elements.
A body-centered three-screw boundary model supplies analytic four-column
characters in rotated Cartesian frames, reordered operation lists, and
cell scales 1e-100, 1, and 1e100. It also rejects duplicate columns, wrong
gauge factors, nonfinite/huge values, and dimension-preserving orthogonal
column mixtures which are not representation characters.

The unitary restriction checks do not replace the separate Wigner test
for extension to antiunitary corepresentations. Reference correspondences
do not yet assign standard labels to arbitrary Gaussian input cells,
identify every onsite EBR column, or complete global connectivity.

## Reproduction

Run from the CP2K research source root, with the affine fixture generated
as described in `affine-reference-validation.md`. The test-only Python
environment uses NumPy, SymPy and `irreptables==3.1.0`.

```sh
python prepare_character_reference.py affine-reference.dat character-reference.dat character-sources.jsonl
OPENBLAS_NUM_THREADS=1 python audit_character_reference.py character-reference.dat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 build-serial/bin/topology_character_unittest.ssmp --reference=character-reference.dat
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 mpiexec -n 2 --tag-output build-mpi/bin/topology_character_unittest.psmp --reference=character-reference.dat
bash build_character_instrumented.sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 build-serial/character-debug/character-reference --reference=character-reference.dat
python verify_character_outputs.py serial.log mpi.log debug.log
```

To reproduce the unconverted-convention audit, generate a separate fixture
with `prepare_character_reference.py --native-k-sign 1` and otherwise the
same positional arguments. The raw audit reports 25 failed tables;
conjugating all columns is an algebra diagnostic only, not the chosen
reference adapter.

The focused instrumented build uses GNU Fortran 16.1 on Apple ARM with
`-fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined` on the
little-group, corepresentation and character-matching modules plus the
test program. It still links the existing CP2K library and OpenBLAS.
It is not a separate full no-SPGLIB CP2K build.

The official driver selection is:

```text
UNIT/topology_(integer|atomic|wyckoff|band|reciprocal|site|character|symmetry)_unittest|QS/regtest-little-group|QS/regtest-topology$|QS/regtest-kp-1$
```

Serial uses one MPI rank, two OpenMP threads, three concurrent tasks;
MPI uses two ranks, two OpenMP threads, two concurrent tasks. Both set
`OPENBLAS_NUM_THREADS=1`, `OMP_STACKSIZE=128M`, `CP2K_DATA_DIR` to the source
data directory, and the validated thread-safe OpenBLAS insertion library
recorded in the logs. All logs are in `character-reference-logs.tar.gz`.
