# Stanene SOC Torus: Corrected Operator Comparison

The 25 September full-band cross-check found an operator-assembly error,
not merely an unconverged invariant. `build_core_ppnl` already reverses
the angular-momentum integrals for swapped atom blocks; the periodic
localizer applied an extra `is_complex` sign during image folding.
The correction uses the same transform as CP2K's existing post-SCF SOC
band path. No existing regression reference or numerical tolerance was
changed. The Wilson calculations themselves are unchanged.

## Independent Matrix Check

All 52 spinor bands on the 3x3 Gamma-centred property mesh are exported
at the same 8x8x1 converged scalar potential, with complete neighbour
overlap matrices. The polar factors of the square AO coefficient
matrices form the orthonormal Bloch frames. A discrete Fourier transform
of the Hamiltonian and coordinate-link matrices reconstructs the full
936-order localizer, compared against the metric-normalized direct AO
matrix. This is independent assembly/eigensolver/Pfaffian code, not an
independent pseudopotential or SCF implementation.

| Check | Before fix | After fix |
| --- | ---: | ---: |
| Maximum diagonal-block discrepancy (Ha) | 3.09907e-2 | 1.84364e-12 |
| Maximum spin-averaged scalar discrepancy (Ha) | 3.87512e-12 | 1.82532e-12 |
| Maximum coordinate-block discrepancy (Ha) | 1.31252e-15 | 1.31252e-15 |
| Native gap (Ha) | 1.91347978651e-4 | 6.26487247807e-4 |
| Native index | 1 | 0 |

The independent gap is 6.26487247768e-4 Ha with index zero. Its Pfaffian
uses orthogonal Hessenberg reduction, cross-checked against recursive
Pfaffians on 32 random skew matrices of orders 2, 4, 6 and 8. No band
flattening, coordinate-link unitarization or occupied-space projection
is used. Research switches for such transformations are not activated.

## Corrected Controls

At E=-0.1592 Ha and ETA=0.03 Ha, DZVP/200 Ry gives indices 0, 1, 0, 0
for torus sides 3, 4, 5, 6. N=4 is zero at ETA=0.025 and 0.035 Ha.
The TZVP/400 Ry N=4 control is zero at E=-0.16050731 Ha, ETA=0.03 Ha.
These results remain incompatible with a converged material localizer
plateau, despite the stable occupied-band Wilson index of one.

For corrected N=3, the dense gap lies inside the independent MUMPS
bracket [6.264872243636e-4, 6.264873096930e-4] Ha. Tacho gives the same
zero index, with solve residual 2.11e-14; the serial/MPI dense gap
difference is 1.40e-13 Ha. The run uses 66 shifted factorizations and
97 MB process-summed solver memory, not whole-process memory.

The small 2x2-SCF/2x2-torus fixture has independent and native gaps
differing by 2.19e-12 Ha. Two new ordinary inputs cover dense and Tacho
execution. After make_pretty, the localizer driver passes 37/37 serial
and 62/62 two-rank MPI assertions (including MUMPS/Tacho). The earlier
two-rank localizer/quadratic run passes all 62 assertions, with no
changes to the old references. These counts overlap and must not be
summed as independent tests.

## Archive and Replay

`stanene-soc-correction.tar.gz` contains 91 members plus manifest,
20,976,628 bytes, SHA-256:
`04054fd5d8e296106466b69271d7e3d8cb65d7cb52ad2f168928718d01c7e0a4`.
It contains literal inputs, all corrected outputs, both complete Bloch
datasets, before/after dense matrix dumps, basis/potential data,
source snapshots, build caches, regression inputs and research scripts.
The source snapshot is local CP2K commit `9891c230a0`; no CP2K push is
implied. Temporary matrix-output instrumentation was removed.

With Python, NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python3 replay_stanene_soc.py stanene-soc-correction.tar.gz stanene-soc-summary.json
```

The replay verifies each member hash, reconstructs both complete Bloch
localizers, compares before/after AO matrices, checks the independent
Pfaffian, and validates corrected serial/dense/sparse agreement. The
summary CSV distinguishes corrected queries from the retained negative
control. Earlier `stanene-localizer-{dzvp,tzvp}` archives remain useful
for the unchanged Wilson references and explicitly historical output.
