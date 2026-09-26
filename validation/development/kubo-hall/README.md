# Antisymmetric charge DC response, 26 September 2026

This checkpoint adds opt-in `HALL_RESPONSE T` to `PROJECTED_AO` and `BLOCH`.
It does not add magnetism, an external magnetic field, spin Hall currents,
extrinsic skew/side-jump scattering or a microscopic collision operator.
Native CP2K changes remain local; no public CP2K branch or PR was pushed.

## Convention and mathematical checks

The first index denotes current and the second electric field. The resolvent
is `(gamma + i[H,.])^-1`, and CP2K's current is anti-Hermitian `J=-i*hbar*v`.
Pairing two states gives
`2*f[en,em]*(en-em)*Im[Ja(m,n)*Jb(n,m)]/(gamma^2+(en-em)^2)`.
The decreasing Fermi divided difference is negative. Exactly degenerate
pairs contribute zero; no inverse band splitting is evaluated.

The native mathematical test checks a two-level analytic response,
near-degeneracy down to 1e-14 model energy, band permutations, energy-zero
shifts, unitary mixing inside two degenerate doublets, Cartesian rotations
and reflections, and time reversal. A full 16x16 complex Liouvillian is
independently assembled and inverted after mixing the four-level Hamiltonian
into a nondiagonal frame. Maximum discrepancy is 7.55e-15 over optimized
and instrumented executions.

For `H=(sin(kx),sin(ky),m+cos(kx)+cos(ky)).sigma`, a 48x48 mesh with
`kT=0.01`, `gamma=0.001` gives:

| m | Lower-band Chern integral | sigma_xy/(e^2/h) |
|---|---:|---:|
| -1 | +1 | -0.99999976562 |
| +1 | -1 | +0.99999976562 |
| +3 | 0 within 5.3e-15 | 9.65149189e-9 |

Every k point is also checked against the analytic finite-temperature,
finite-relaxation Berry response. The maximum discrepancy is 7.22e-16.
Repeating all three meshes after k-dependent, nonunitary AO basis changes
with transformed H/S derivatives and connection changes the response by at
most 2.34e-15. The model Berry convention is `A=i<u|d u>`; here the cold
insulating limit is `sigma_xy=-(e^2/h) C`. These are model results, not DFT
material Hall conductivities.

## Native validation

- SSMP: 67/67 assertions, including the Hall unit and three Kubo directories.
- Four MPI ranks: the same 67/67 assertions.
- Two MPI ranks: 14/14 assertions in the final Hall-only directory/unit subset.
- Two OpenMP threads, one BLAS thread, thread-safe OpenBLAS; 64 MiB OpenMP
  stack for SSMP/two ranks and 128 MiB for four ranks.
- Both existing projected/Bloch operator unit executables also pass.
- The new unit passes a standalone Fortran 2008 build with `-O0 -g`,
  `-fcheck=all` and invalid/division-by-zero/overflow traps.
- `make_pretty.sh`: 20 files, zero failures. Both full CP2K builds succeed.

New Ne/GTH-SOC controls cover finite GPW/GAPW, periodic GPW/GAPW, full-grid,
K290 and SPGLIB property reconstruction, and 1D/2D normalization. UZH MOLOPT
bases are used throughout. Time reversal requires zero charge Hall response:
the largest 3D/finite residual is below 1e-11 S/m; the sheet residual is below
1e-21 S; the strictly 1D tensor is exactly zero. This does not test a nonzero
magnetic DFT Hall response. Existing Bi/diamond and lower-dimensional Kubo
references are unchanged.

Seven additional native calculations compare Hall on/off at identical input:
3x3x3 scalar Bloch, 3x3x3 SOC Bloch, and finite projected SOC, plus rejection
of Hall with ATOM_EMBEDDING. All three pairs give identical printed isotropic
conductivity, and full-tensor differences are at most 3.913e-12 S/m.

The initial draft of the new regression fixture accidentally compared a
3x3x3 grid against a pre-existing 2x2x2 reference. Its initial 1e-12 S/m Hall
zero tolerance was also tighter than the finite-SOC eigensolver cancellation
noise. That failing log is retained explicitly. The fixture now uses the
matching 2x2x2 grid; existing physical references were not adjusted. Its new
absolute 3D Hall tolerance is 1e-9 S/m, below 1e-10 of the dissipative scale.
Separate 3x3x3 on/off controls verify the original mesh independently.

## Frozen evidence and replay

`methods-results.tar.gz` contains numerical sources, native inputs/outputs,
all final driver logs, the initial fixture-failure log, controls, formatting
and build evidence, and the native patch. `manifest.json` records every
SHA-256 and the local native checkpoint. The GPL source license is included.

Verify frozen evidence:

```bash
python3 replay.py
```

Rebuild the mathematical unit directly from the archived sources, with
bounds and floating-point checks, against a locally installed BLAS/LAPACK:

```bash
python3 replay.py --compiler gfortran --blas-library /absolute/path/to/libopenblas.a
```

The replay verifies recorded native results; it does not rebuild all of
CP2K. Re-executing the native fixtures requires the recorded local source
checkpoint, CP2K data files and corresponding build dependencies. Full
material convergence, microscopic scattering and noncollinear magnetic
DFT Hall benchmarks remain separate work.
