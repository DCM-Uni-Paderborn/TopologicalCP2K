# Full invariant-torus character validation

Validated on 25 September 2026 with local CP2K commit
`77e51c818026749b0bb5a779a2ed80fecce5ddf0`. No CP2K push was made.
This extends the site-resolved EBR comparison with a production diagnostic,
not with a canonical EBR-label or Bloch-bundle classifier. It changes no
Hamiltonian, input keyword, energy reference or test tolerance.

## Mathematics and implementation

For a unitary direct-space rotation W, the reciprocal fixed group is
K_W={k mod Z^3: (transpose(W)-I)k in Z^3}. Its discrete character group
is Z^3/(W-I)Z^3. The checked Smith decomposition U(W-I)V=D supplies
canonical keys U*L for each integer Bloch shift L: reduce rows with
nonzero D modulo the positive divisor and retain free rows unchanged.
Collect the complex onsite/spin coefficients of identical keys. Distinct
keys are independent characters on the complete compact group K_W,
including disconnected torsion components.

`atomic_band_character_polynomial` builds these finite sums for every
unitary operation in the common group, cell, origin, coset and spin-lift
convention. Antiunitary orbit representatives conjugate onsite characters;
antiunitary operators themselves are not assigned traces. Integer products
and Smith factorizations use the existing checked 64-bit implementation.
`atomic_character_compare` checks finite coefficients, storage shapes,
unique sorted keys and common quotient conventions before comparing them.
Tolerance is dimension-scaled, not a floating-point rank decision.

The atomic catalogue first uses its integer sampled signatures to discard
pairs that already differ. Remaining pairs undergo the full-function test.
Nontransitive tolerance neighborhoods are rejected. The detailed Gaussian
output contains `ATOMIC_CHARACTER_CLASS column representative` and
`ATOMIC_CHARACTER_CLASS_RESIDUAL`. Representatives are local column IDs,
not canonical labels. All columns, integer membership decisions and
real-space induction witnesses remain unchanged. Memory estimates include
the transient Fourier cache and pairwise comparison storage.

## Results

The analytic inversion test distinguishes the two parity pairs at centers
(0,0,0) and (1/2,0,0), although their Gamma characters coincide. A
Gamma-only catalogue retains 17 distinct full-torus classes: sixteen
center/parity columns and one generic orbit, also with Kramers pairing.
Absent keys, nonfinite coefficients and empty input are rejected.

All 378 multi-candidate ordinary EBR references from the preceding
230-space-group audit have equal unitary-character functions within their
candidate sets. This holds in both original and combined lattice-shear,
Cartesian-rotation, origin-shift, operation-order and coset conventions.
Denser reciprocal sampling cannot separate these cases. The result does
not prove bundle equivalence or identify canonical onsite axes.

| Run | Prepared bands | Fourier/raw-induction queries | Maximum dimension-scaled residual |
| --- | ---: | ---: | ---: |
| Optimized serial, EBR fixture and bounded controls | 3,248 | 213,238 | 3.9740e-15 |
| MPI rank 0, same audit | 3,248 | 213,238 | 4.0279e-15 |
| MPI rank 1, same audit | 3,248 | 213,238 | 4.0279e-15 |
| Focused instrumentation, EBR and analytic controls | 3,200 | 205,614 | 4.0279e-15 |
| Serial, all 530 Hall settings and analytic controls | 6,362 | 536,410 | 3.7741e-16 |

Raw induction is evaluated at every torsion component and at three free
coordinate samples y, using k=transpose(U)*y. Those numerical samples
validate the Fourier implementation; equality is decided from the finite
coefficient lists, not from the numerical samples. The full Hall sweep
includes ordinary/grey and scalar/spinful groups at three seeds and
114,480 independent Cartesian s/p/d orbital-spin comparisons. It is a
symmetry-algebra audit, not 530 DFT material calculations.

The atomic-signature all-setting audit passes in serial and on both MPI
processes: 2,120 matrices, 31,079 columns, 707,038 entries, 217,553
column-edge checks and 33,199 nonnegative witnesses. MPI processes repeat
the algebra independently; no distributed scaling result is claimed.

The relevant native, Gaussian, topology and k-point driver passes 179/179
checks in each optimized build. Both EBR negative-test runs accept the
baseline and reject all six deliberate corruptions. Independent Gaussian
export validation checks 16 files per build, 197 sites, 697 columns,
263 queries and 31,537 induced character entries. Maximum induction
residual is 3.57597e-14; named point-group reference residual is
9.55958e-6 at the unchanged rounded-table tolerance. All 128 corrupted
onsite/spin/product records are rejected across the two builds.

The new exported class bookkeeping passes independently. Examples in
`fourier-character-classes.json` include He GPW (125 columns, 18 sampled
classes, 125 full-torus classes), He GAPW (11,5,11), Ne SOC (33,17,29)
and the Ne spin-lift case (26,23,23). The largest accepted class
coefficient residual is 1.85399e-17. These count reference columns, not
occupied bands, and are not permission to merge Bloch bundles.

## Reproduction and scope

Use the fixture and original reference data archived in
`ebr-character-records.tar.gz`. The companion
`fourier-character-records.tar.gz` retains the current fixture, changed
source files and patch, build/format/test logs, Gaussian detail outputs,
negative fixtures, audit scripts and a SHA-256 manifest. It is an evidence
archive, not a standalone CP2K source distribution. Restore it under a
matching source checkout with the dependencies recorded in the environment
JSON. The CP2K branch has not been published.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  build-serial/bin/topology_band_unittest.ssmp \
  --ebr-reference=build-serial/ebr-character-full.dat > serial.log
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  mpiexec --tag-output -n 2 build-mpi/bin/topology_band_unittest.psmp \
  --ebr-reference=build-serial/ebr-character-full.dat > mpi.log
python3 audit_fourier_characters.py serial.log mpi.log
python3 audit_ebr_characters.py build-serial/ebr-character-full.json serial.log mpi.log
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  build-serial/bin/topology_band_unittest.ssmp --all-settings
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 OMP_STACKSIZE=64M \
  build-serial/bin/topology_atomic_unittest.ssmp --all-settings
```

`build_fourier_instrumented.sh` rebuilds the test, symmetry, little-group,
corepresentation, character-matching, integer-lattice and induction
kernels with bounds checking, floating-point traps and undefined-behavior
instrumentation. Other dependencies remain serial release-library code.
This executable does not compile the optional SPGLIB database driver and
is not a fully instrumented Gaussian or MPI build. Its output nevertheless
checks the supplied EBR fixtures for all 230 ordinary groups.

The evidence includes ordinary EBR matching, not grey EBR-reference
certification, canonical onsite labels, elementarity of the 63 unmatched
native columns, global connectivity, converged material topology, or
multi-node performance. `make_pretty.sh --no-cache` passed all six
changed CP2K files. No regression reference values were edited.

The evidence archive contains 232 files (39,216,301 uncompressed bytes;
2,459,290 compressed bytes). Its SHA-256 is
`5ab7973e87ae4af0b20bd1285b161619816e5a5ffe43042206f9d8680bc6f35a`.
A fresh temporary extraction verifies all 232 file hashes and reruns the
native fixture and the independent log audit successfully. The replay
again finds 378 equal-function candidate sets and the serial residual
3.97399e-15; the owned temporary directory is removed afterward.
`fourier-character-replay.json` records the result and
`replay_fourier_archive.py` supplies the recipe. Original external-table
sources and their provenance remain in the earlier EBR evidence archive.
