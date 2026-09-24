# Shared spin lifts in Gaussian SOC sewing

Date: 2026-09-24

Implementation: local CP2K commit `1a84c86d65753221fda8c18f9c54cc8044bf141b`,
parent `590d9ee878dae35a33497263a55f821421582715`.
The CP2K branch has not been pushed. This report and the manuscript
are synchronized only to the private paper repository.

## Defect and correction

`little_group_factor_system` constructs the spin double-group factors
using a normalized cell. Previously `write_little_group` independently
computed a Cartesian rotation using `cell%hmat` and `cell%h_inv`, and then
called `spin_rotation` again. At a tied quaternion branch, rounding can
choose a different global SU(2) sign. Individually valid spin matrices
then fail the product law of the independently generated factor table.
This is not an SCF energy or EBR-reference-value error.

The Gaussian path now requests the actual `spin_lifts` from the factor
constructor. CP2K's conjugate physical-spin representation is retained.
For antiunitary operations, the stored lift already contains
`J = i sigma_y`. The columns are converted back from `U J` to `U`
before `transform_spinor_frame` applies `J K` to the source. Thus Theta
is applied once, not twice. The extra four complex entries per operation
are included in the existing memory estimate. No tolerance, input
default, factor construction or old physical reference was changed.

## Gaussian before/after reproduction

The new directory-regtest input `neon-spin-lift.inp` uses:

- One Ne atom, primitive BCC cell with conventional lattice constant
  4.2 Angstrom and roundoff-sized input-cell perturbations.
- PBE, `BASIS_MOLOPT_UZH`, `DZVP-MOLOPT-GGA-GTH-q8`,
  `GTH_SOC_POTENTIALS`, and post-SCF GTH SOC.
- Gamma SCF and Gamma NNKP property point, complex MOs.
- All 13 scalar states for second variation; the lowest two spinor
  states form the isolated represented Kramers doublet.
- One 1200 Ry grid; `EPS_SCF 1e-10`; unchanged default little-group
  tolerance 1e-6 and energy-commutator tolerance 1e-6 Hartree.

The parent executable fails **at the product-law check**, after passing
metric and subspace closure. The before/after input differs only by
prettify whitespace and keyword ordering, not physical settings. Direct-run results:

| Quantity | Parent | Shared lifts |
| --- | ---: | ---: |
| SCF energy / Hartree | -34.90553744771596 | -34.90553744771598 |
| Metric residual | 1.21322e-15 | 1.21322e-15 |
| Maximum frame residual | 2.89849e-7 | 2.89849e-7 |
| Projective product residual | 2.00000 | 8.49321e-14 |
| Energy commutator / Hartree | 1.77636e-15 | 8.88178e-16 |

The fixed path identifies one 2D irrep, one undoubled 2D corepresentation,
and both signed and nonnegative atomic membership. The official two-rank
MPI run gives product residual 8.48210e-14 and energy
-34.90553744771597 Hartree. The official SSMP run gives 8.39331e-14
and -34.90553744771599 Hartree. These small differences are numerical
reduction/roundoff changes, not adjusted reference values.

## Distinct convergence effects

At 400 Ry the Gamma primitive-cell input fails subspace closure at
9.19484e-5, before the product-law check. The 1200 Ry single-grid run
reduces this residual to 2.89849e-7. No closure guard was weakened.

The initially attempted primitive-cell 2x2x2 shifted MP mesh preserves
only 8/48 cubic spatial operations. Its occupied subspace is not a
full-cubic reference even if the atomic geometry is cubic. The grid
count can be reproduced with the standalone script archived alongside
this report. With that mesh, raising the real-space cutoff to 3600 Ry
does not remove the remaining closure defect (7.24275e-6 in the
retained two-band diagnostic). This is a separate integration-set
issue and was not hidden by fitting band transformations.

A nearly diagonal diagnostic cell is automatically made exactly
orthorhombic by CP2K. It consequently does not reproduce the branch
tie. The final regression intentionally keeps the nonorthogonal
primitive basis. The complete 26-band case correctly encounters the
existing requirement for an isolated selected subspace and computed
excluded states; that guard was not changed either.

## Verification

- Both serial and MPI executables rebuilt.
- `make_pretty.sh --no-cache` on the five changed files: 5 checked,
  0 failed.
- Standard drivers: **179/179** assertions in SSMP (one rank, two OMP
  threads) and PSMP (two ranks, two OMP threads each).
- Selection: eight topology unit programs, 135 Gaussian little-group
  assertions, 36 existing topology/export/k-point assertions. This
  includes the existing nonsymmorphic, antiunitary, automatic geometry,
  GAPW, and Wannier90 checks.
- New bounded test applies the 48 cubic polar operations and their
  grey extension to complex frames in a primitive basis; the cell is
  scaled by 4.2 times 1e-100, 1 and 1e100.
- Focused GNU Fortran build of the spin/factor modules and symmetry
  unittest passes with `-fcheck=all`, traps for invalid/zero/overflow,
  and undefined-behavior instrumentation. This focused program does
  not enable the optional SPGLIB database sweep and links the existing
  CP2K library. It is not a fully instrumented Gaussian/MPI build.
- The new Gaussian test takes about 42 s in the concurrently running
  local serial/MPI regression jobs; this is not a scaling benchmark.

The log archive includes the parent failure, direct fixed run, final
serial/MPI driver logs and new-case outputs, the focused unit result,
prettify/build logs and the independent lift diagnostic. It deliberately
retains unsuccessful convergence probes for traceability.

## Remaining scope

This fixes production Gaussian SOC sewing in the current topology
branch. It does not supply automatic canonical EBR naming, a full
magnetic SCF treatment or general band connectivity. The independent
all-space-group character correspondence and the real-space site
catalogue remain the basis for the next EBR-identification step.
