# Site-induced atomic band validation

Verified 24 September 2026. These are native algebra/reference tests, not
new material calculations or a completed general EBR classifier.

## Full sweep

- `atomic-band-serial-all-settings.log`: serial release executable.
- `atomic-band-mpi-all-settings.log`: two MPI ranks, each exercising the full kernel.
- `atomic-band-debug-all-settings.log`: bounds checking and invalid/zero/overflow traps.

Every sweep covers all 530 Hall settings, three fractional site seeds, six
reciprocal points, scalar/spinful factors, and ordinary/grey groups. It reports
38,160 site/k-point cases, 67,824 induced-representation decompositions and the
same number of directed compatibility checks. Independent Cartesian s/p/d
orbital-spin characters give 114,480 matching comparisons.

The seeds are (0,0,0), (1/2,0,1/4), and (0.173,0.217,0.319). Reciprocal points are
Gamma, (1/2,1/2,1/2), (1/2,0,0), (1/3,1/3,0), (0,0,1/4), and
(0.173,0.217,0.319). Connections start at Gamma and include its (0,0,1) image.
This samples special and general site orbits; it is not exhaustive over Wyckoff
positions or all Brillouin-zone manifolds.

## Additional checks

Analytic checks cover inversion-center parities, generic inversion orbits,
nonsymmorphic screws, undoubled quaternion spin irreps, and C3 orbital partners
whose orbit representatives include an antiunitary half translation. Permuted
input ordering, invalid/empty inputs and nonfinite coordinates are checked.
Further tests verify passive origin changes, sheared unimodular cells and
the Bloch rephasing under integer changes of coset representatives.

The official focused directory driver checks the new native unit program and
the existing Gaussian little-group directory. Its 48 assertions are distinct
from the database counts; passing those assertions does not constitute a new
set of 48 material calculations. The regular CI database sample is deliberately
bounded; `--all-settings` retains the complete sweep above.

## Reproduction and scope

Build `topology_band_unittest` with the CP2K CMake build. Run its versioned
executable with `--all-settings` for the full sweep; without that switch it
runs the analytic tests and bounded CI sample. SPGLIB is required only for
the database sweep, not for construction from supplied operations. A separately
compiled test without the SPGLIB preprocessor option passes the analytic suite.

Local Apple ARM runs use one OpenBLAS thread and the retained thread-safe
OpenBLAS installation. Gaussian MPI regtests additionally use two OpenMP
threads per rank and OMP_STACKSIZE=128M. The debug executable is compiled with
`-O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`.

The implementation generates characters induced from supplied positions and
onsite irreps/corepresentations. It has no new production input keyword, does
not determine elementary/composite status, and does not supply an exhaustive
Wyckoff catalogue or a general topological classification. Integer onsite
content and compatibility cannot by themselves prove a trivial Bloch bundle.
