# Spectral localizers

`FORCE_EVAL/PROPERTIES/SPECTRAL_LOCALIZER` evaluates a two-dimensional localizer from the converged
**full AO Hamiltonian**, overlap and analytic Cartesian position integrals. Assembly uses
distributed DBCSR matrices, represented by real/imaginary components in CP2K's matrix interface.
`SOLVER DENSE` provides the LAPACK reference; optional `SOLVER MUMPS` supplies sparse class A
inertia and metric-aware gap brackets without dense AO matrices. Optional `SOLVER TACHO` adds sparse
class AII/Z2 Pfaffians with independent MUMPS gap validation. It supports isolated GPW and GAPW
calculations, with either a scalar collinear Hamiltonian or a restricted Hamiltonian augmented by
GTH spin-orbit coupling (SOC). With the default `FORMULATION FINITE`, both `SUBSYS/CELL/PERIODIC`
and `DFT/POISSON/PERIODIC` must be `NONE`. The optional `FORMULATION PERIODIC` supports
two-dimensionally periodic Gamma supercells and explicit full analysis meshes using Berry integrals,
as described below. The periodic analysis may start from a symmetry-reduced k-point SCF, but builds
the complete analysis torus rather than reducing its property operators. Unconverged SCF
calculations are rejected.

## Input and related analyses

For an isolated system, add to `FORCE_EVAL`:

```text
&PROPERTIES
  &SPECTRAL_LOCALIZER ON
    INVARIANT CHERN
    SOLVER DENSE
    PLANE XY
    POSITION [angstrom] 0.0 0.0 0.0
    ENERGY [eV] -0.1 0.0 0.1
    KAPPA [hartree*bohr^-1] 0.01 0.02
    EPS_GAP [hartree] 1.0E-8
    EPS_METRIC 1.0E-10
    MAX_AO 512
  &END SPECTRAL_LOCALIZER
&END PROPERTIES
```

`POSITION` is repeatable; all combinations with the `ENERGY` and `KAPPA` lists are evaluated.
Coordinates are Cartesian in the input structure's coordinate frame, not fractional. The component
perpendicular to `PLANE` is ignored. Energies are absolute Hamiltonian energies, **not relative to
an automatically selected Fermi energy**. Choose them using the spectrum and converge the positive
energy/length scale `KAPPA`. `XY`, `XZ`, and `YZ` specify oriented planes; changing orientation
reverses a Chern sign.

This analysis is independent of the sibling `KUBO_TRANSPORT` section. Its analytic position
integrals do not change Kubo's default atom-embedded current approximation. The optional finite
`KUBO_TRANSPORT/CURRENT_OPERATOR PROJECTED_AO` uses these integrals with a metric-aware commutator;
see the limitations below. The periodic Wilson loops, Chern numbers and inversion-subgroup
indicators under `DFT/PRINT/WANNIER90` are also unchanged. The localizer needs neither Wannier90 nor
SPGLIB, does not require crystal symmetry, and does not identify elementary band representations. A
resolved finite localizer index is not a certificate that a bulk DFT material or a metallic occupied
manifold has a particular band topology.

The separate [quadratic pseudospectrum](quadratic_pseudospectrum.md) shares the finite AO operators
and finds joint position/energy approximate states with separated residuals. It reports
localization, not a Chern or Z2 index, and can be enabled alongside the spectral localizer.

For physical spinful time reversal, use instead:

```text
INVARIANT Z2
SOC T
TIME_REVERSAL T
```

This requires a restricted SCF and **GTH potentials containing SOC parameters for all atomic
kinds**, for example from `GTH_SOC_POTENTIALS`, with valence-compatible bases. Ordinary
scalar-relativistic GTH potentials are insufficient. The existing CP2K SOC operator is added in the
complete AO spinor space; this is **not self-consistent noncollinear SOC DFT**. The time-reversal
relations of the Hamiltonian, overlap and positions are checked, not merely assumed from the
keyword. `Z2` reports 0 for the trivial and 1 for the nontrivial class AII local index. This is not
an inversion-parity calculation. A real scalar Hamiltonian without time-reversal breaking cannot
serve as a nonzero class A Chern benchmark.

## Nonorthogonal metric and index

Let $H_{ab}=\langle a|\hat H|b\rangle$, $S_{ab}=\langle a|b\rangle$ and
$X_{ab}=\langle a|\hat x|b\rangle$, $Y_{ab}=\langle a|\hat y|b\rangle$ in the same finite basis. The
covariant localizer at $(x,y,E)$ is

$$
L=\begin{pmatrix}
H-ES & \kappa[(X-xS)-i(Y-yS)]\\
\kappa[(X-xS)+i(Y-yS)] & -(H-ES)
\end{pmatrix}.
$$

The protection gap is $\mu=\min_j|\lambda_j|$ for $Lv=\lambda\,\mathrm{diag}(S,S)v$. These
generalized eigenvalues are invariant under a nonsingular AO basis change. Using raw AO eigenvalues
would give a basis-dependent protection gap. No occupied-space truncation, energy-window projection
or overlap-eigenvector removal is performed. An indefinite or poorly conditioned metric is rejected
according to `EPS_METRIC`, the minimum allowed smallest/largest overlap eigenvalue ratio. The sparse
solver uses a conservative numerical lower bound on this ratio, rather than its exact value.

For class A, the index is half the signature of $L$. If the gap is below `EPS_GAP` or a metric- and
spectrum-dependent roundoff threshold, output is **UNRESOLVED**, not a reliable zero index. All
indices require convergence with system size, boundary geometry, basis and a suitable interval of
`KAPPA`. The projected AO position operators need not commute exactly in a finite Gaussian basis.
The reported gap is for that finite projected model; basis incompleteness is not covered by its
numerical tolerance.

For class AII, physical spinors are ordered with all up AOs followed by all down AOs, and
$\Theta=JK$, $J=(0,I;-I,0)$. With the auxiliary localizer index outermost, $C=J_{\rm aux}\otimes J$
and $Q=(I-iC)/\sqrt{2}$ transform $L$ into a purely imaginary skew matrix. Dense pivoted skew
elimination evaluates the sign of the real Pfaffian of $iQ^\dagger LQ$. Its orientation is
calibrated against $L_0=\mathrm{diag}(S,-S)$ in the same basis. Only signs are accumulated, avoiding
Pfaffian-product overflow. This dense routine is independently implemented; it is not an integration
of the external sparse Pfaffian package.

### Dense and sparse factorization

Both paths assemble the same covariant operator in DBCSR, retaining the union of Hamiltonian,
overlap, position and SOC sparsity patterns without additional dropping. The localizer dimension is
twice the scalar AO count without SOC and four times it with SOC. `MAX_AO` limits the scalar AO
count for both paths; increase it explicitly for larger sparse calculations.

`DENSE` gathers the assembled localizer and metric on each MPI rank and evaluates on the source
rank. Three mathematically distinct operations are used:

| Quantity             | Dense operation                                                                   | Sparse operation                                                             |
| -------------------- | --------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| Chern half-signature | Pivoted `ZHETRF` ($LDL^\dagger$), or `DSYTRF` ($LDL^T$) for exactly real matrices | MUMPS symmetric-indefinite inertia after realification                       |
| Protection gap       | Generalized Hermitian eigenvalues of $(L,B)$                                      | Inertia counts for shifted pencils $L\pm\sigma B$                            |
| AII/Z2 index         | Pivoted real skew elimination and reference-calibrated Pfaffian sign              | Tacho skew-LDL with permutation signs, plus independent MUMPS gap validation |

The dense Chern index comes from the LDL inertia, cross-checked against the generalized spectrum.
Both [LAPACK factorizations](https://www.netlib.org/lapack/complex16/zhetrf.f) have 1-by-1 and
2-by-2 diagonal blocks. Shared block evaluation counts their eigenvalue signs, including
zero-diagonal 2-by-2 blocks. The [real path](https://www.netlib.org/lapack/double/dsytrf.f) never
drops a small nonzero imaginary part to select real arithmetic. Pivot magnitudes are **not**
protection gaps. Ordinary inertia and determinant signs are **not** Pfaffian signs: an odd
simultaneous row/column permutation preserves the determinant but reverses the Pfaffian. No
additional keyword is needed for these distinctions.

Dense memory is quadratic and time cubic; this path is the small-system oracle, not a scalable
solver.

`MUMPS` requires an MPI build with [MUMPS](../../technologies/mumps.md). It checks Hermiticity
before converting the distributed blocks to lower-triangular coordinates. Realification of complex
Hermitian matrices doubles their dimension and duplicates their inertia; the index compensates for
this duplication. There is no dense gather and no $S^{-1/2}$.

The metric is checked separately with the non-pivoting positive-definite MUMPS mode. Internally
MUMPS uses $LDL^T$ here, so a successful factorization is accepted only when there are also **no
negative pivots**. Positive shifted metrics provide a lower bound on their smallest eigenvalue; a
row-sum norm bounds the largest eigenvalue. This deliberately conservative test can reject a metric
that the exact dense condition ratio would admit.

The indefinite localizer uses numerical pivoting, with static pivot perturbation and compression
disabled. Inertia counts of $L\pm\sigma B$ bracket the nearest generalized eigenvalue. The printed
`Gap` is the lower end, and `Gap bracket` gives both ends. These are finite-precision numerical
bounds, not interval-arithmetic certificates. Symbolic analysis is reused within the metric stage
and within the indefinite stage. A too-small or unresolved gap never produces an index.

Output reports factorization counts, input entries (including additive coordinate duplicates),
factor entries, factor workspace, and MUMPS's process-summed solver memory. The latter excludes
CP2K/DBCSR and application buffers and is not whole-process RSS. Sparse fill-in and a sequential
root front can limit MPI scaling; sparse is not automatically faster for small systems.

The `MUMPS` backend covers class A, including complex SOC matrices. For `INVARIANT Z2`, use
`SOLVER DENSE` or optional [`SOLVER TACHO`](../../technologies/tacho.md). Ordinary symmetric inertia
does not determine a Pfaffian sign. Tacho applies the AII basis transformation and skew-LDL
factorization sparsely, with explicit permutation signs and atomic-reference orientation. Pfaffian
factors reside on the source rank, not across MPI ranks. Factorization solve-probe residuals and
symbolic 2x2-block counts are reported separately from MUMPS diagnostics. Unresolved factors yield
no index; there is no automatic dense fallback.

## Export and validation

For external finite-Hamiltonian analysis, `DFT/PRINT/AO_MATRICES` now also accepts `POSITION T`.
With `OVERLAP T`, `KOHN_SHAM_MATRIX T` and `NDIGITS 16`, it prints the three analytic position
matrices in Bohr about the coordinate origin, in the same AO ordering as the other matrices. This
export does not require activating the localizer. The exported KS matrices are scalar collinear
matrices; `SOC T` in the localizer does not silently change their meaning or export the assembled
spinor Hamiltonian.

The standalone matrix kernel also accepts other finite Hermitian Hamiltonians if a consistent
metric, position operators and, for AII, physical time-reversal representation are supplied.
Automatic DFTB/xTB adapters and external-Hamiltonian file input are not implemented. Neither a
velocity matrix nor $[H,r]$ alone determines all position matrix elements, especially diagonal and
degenerate-subspace components.

Tests include nontrivial/trivial open Qi-Wu-Zhang models, spin-mixed time-reversed copies,
orientation reversal, complex nonorthogonal basis changes, origin/energy shifts, singular metrics,
gap closure, broken time reversal and Pfaffian permutations. Additional LDL tests cover
zero-diagonal 2-by-2 pivots, every inertia signature for orders 1 through 12, real/complex
congruences, rescaling from 1e-150 to 1e150, and tiny sign changes. Independent recursive Pfaffian
polynomials for orders 2 through 8 check negative scaling, odd permutations and negative-determinant
basis changes. An analytic generalized-gap example deliberately separates eigenvalues from pivot
sizes, before and after a nonorthogonal basis change, with separate tests on both sides of zero. The
sparse unit test compares indices and generalized gaps against the dense oracle, including complex
metrics, singular/indefinite metrics and gap closure. Small molecular GPW/GAPW and GTH-SOC regtests
exercise the CP2K integration; these are not nontrivial DFT-material benchmarks.

MPI-distributed Pfaffian factorization, large SOC flakes, and basis/size-converged material
benchmarks remain separate developments. Periodic coordinates must not be substituted into the
open-boundary formula; use the torus construction below. A suitable material validation sequence is
crystalline bismuthene with SOC, followed by disordered Bi flakes; it requires independent
band-topology references and systematic finite-size convergence.

See [Dixon, Loring and Cerjan (2023)](https://doi.org/10.1103/PhysRevLett.131.213801) for local
classification in gapless surroundings and [Wong et al. (2026)](https://doi.org/10.1103/6hj9-tgct)
for parity-valued local classification and sparse Pfaffian methods.

## Periodic torus formulation

For a two-dimensional periodic calculation, select:

```text
&SPECTRAL_LOCALIZER ON
  FORMULATION PERIODIC
  MP_GRID 3 1 2
  PLANE XZ
  ENERGY [hartree] -0.1
  ETA [hartree] 0.5 1.0
  POSITION [angstrom] 3.0 3.0 3.0
&END SPECTRAL_LOCALIZER
```

Here both `CELL/PERIODIC` and `POISSON/PERIODIC` must be `XZ`. `XY` and `YZ` are also supported. For
this formulation, `PLANE` selects **lattice-vector indices**, and reciprocal vectors are obtained
from the actual inverse cell, not from Cartesian box lengths. This includes oblique cells supported
by the underlying electrostatics. `POSITION` remains Cartesian; its origin phases are
$\phi_j=\mathbf G_j\cdot\mathbf r_0$, where $\mathbf G_j$ belongs to the enlarged torus cell. Moving
the query by a periodic torus vector leaves the operator unchanged. The nonperiodic cell direction
must remain sufficiently large for the electronic-structure calculation.

The implementation uses the periodic construction of
[Doll, Loring and Schulz-Baldes (2025)](https://doi.org/10.1007/s11040-025-09508-0), Eq. (1). CP2K
integrates $C_j=\langle a|\cos(\mathbf G_j\cdot\mathbf r-\phi_j)|b\rangle$ and
$T_j=\langle a|\sin(\mathbf G_j\cdot\mathbf r-\phi_j)|b\rangle$ analytically. These are **not**
matrix functions of the finite projected Cartesian AO positions. In CP2K's existing $(H,X,Y)$
orientation convention the covariant matrix is

$$
L_{\rm per}=\begin{pmatrix}
H-ES-\eta(2S-C_1-C_2) & \eta(T_1-iT_2)\\
\eta(T_1+iT_2) & -H+ES+\eta(2S-C_1-C_2)
\end{pmatrix}.
$$

This is $-\eta$ times the paper's Eq. (1), conjugated by $\mathrm{diag}(I,-I)$. That convention
retains the Chern orientation of CP2K's finite localizer. `ETA` has energy units; the generalized
gap of $(L_{\rm per},\mathrm{diag}(S,S))$ is consequently still in Hartree. Explicit `KAPPA` is
rejected in this mode, and explicit `ETA` is rejected in `FINITE` mode. Both lists can be scanned
along with `ENERGY` and `POSITION`. The overlap, Hamiltonian and all trigonometric operators use the
same periodic AO basis and boundary-crossing matrix elements.

The dense, MUMPS and Tacho paths share the DBCSR assembly. Post-SCF GTH SOC and class AII are
available under the same restricted-SCF and time-reversal conditions as for the finite localizer.
The periodic mass term retains the required skew structure. Pfaffian orientation is calibrated to
the trivial atomic limit in the same basis, and is tested against spin-mixed pairs of periodic Chern
models on both sides of the transition. This does not implement self-consistent SOC DFT.

**Interpretation matters:** The paper's bulk index theorem assumes a gapped short-range lattice
Hamiltonian with appropriate energy and volume bounds. Projection into a finite Gaussian basis does
not automatically establish those assumptions; projected exponential operators need not be unitary
or commute exactly. CP2K therefore reports a **finite AO torus diagnostic**, not a theorem
certificate. Converge the basis, supercell and `ETA`, with independent band-topology checks for a
material claim. The reported localizer gap is not the electronic band gap. No automatic occupied
space projection, spectral flattening or adjustment of reference values is performed.

The nontrivial matrix tests cover opposite QWZ Chern phases, a trivial phase, spin-mixed AII,
integer-cell and lattice-site shifts, reversed plane orientation, complex nonorthogonal basis
changes and energy rescaling. A full Fourier transform verifies the same gap and index in real space
and k space, including the torus-coordinate couplings between different k sectors. Periodic GPW/GAPW
SOC regression inputs compare all three solvers and exercise `XY`, `XZ` and `YZ`; their neon system
is a trivial integration test, not a topological material benchmark.

`MP_GRID` sets a **full Gamma-centered property mesh**, independent of the SCF scheme, shift or
symmetry reduction. The default `1 1 1` retains the Gamma-supercell analysis. Nonperiodic directions
must have size one. Explicit k-point SCF uses its image-resolved, converged Hamiltonian directly; a
Gamma-only SCF is converted to image-resolved operators at its frozen potential. Neither path reuses
occupied MOs or replaces the SCF density by a differently sampled density.

The implementation transforms $H(k)$ and $S(k)$ back to the associated Born-von-Karman real-cell
basis using CP2K's existing image Fourier convention. Ordered analytic Berry integrals are lifted
with both their ket-cell displacement and bra-cell translation phase. All full-mesh sectors are
retained. This is essential: multiplication by $\exp(i\mathbf G\cdot\mathbf r)$ connects different k
sectors. Independent localizers at irreducible k points with SCF weights are not equivalent. The
real-cell basis also retains the standard physical time-reversal representation for the Pfaffian
path; no unaccounted k-to-minus-k permutation is hidden in that calculation.

The construction is distributed over AO blocks in DBCSR and works with the dense, MUMPS and Tacho
solvers. The current reference Fourier assembly takes $O(N_k^2)$ primitive-block operations.
`MAX_AO` limits the **whole torus**: scalar AOs per primitive cell times the product of `MP_GRID`.
This is not yet a performance claim for dense k meshes. A symmetry-reduced SCF is accepted, but the
property calculation itself is not symmetry reduced.

Wilson/TRIM property symmetry reduction, general space-group TQC/EBR classification, and
three-dimensional AII localizers remain separate blocks.

## Finite projected-AO Kubo current

The sibling property section can use the shared analytic positions and post-SCF GTH SOC:

```text
&KUBO_TRANSPORT ON
  CURRENT_OPERATOR PROJECTED_AO
  SOC T
  TEMPERATURE 300
  DISSIPATION [K] 1000
&END KUBO_TRANSPORT
```

This requires a converged finite GPW/GAPW calculation (`PERIODIC NONE`, no explicit k-points). SOC
additionally requires restricted SCF and SOC parameters in every GTH potential. Complex spinors span
the full AO space, with occupation one per spinor, rather than two as in restricted scalar
calculations. `MAX_AO` bounds the scalar dimension before allocating dense matrices.

For the generalized eigensystem $HC=SCE$ and analytic covariant position matrices $X_\alpha$, the
anti-Hermitian commutators are evaluated without forming an inverse overlap:

$$
J^\alpha_{mn}
= [C^\dagger(HS^{-1}X_\alpha-X_\alpha S^{-1}H)C]_{mn}
= (\epsilon_m-\epsilon_n)[C^\dagger X_\alpha C]_{mn}.
$$

The velocity of this projected model is $iJ$ in atomic units. The conductivity contraction retains
complex products, taking their real part for the **dissipative symmetric response**. No anomalous
Hall conductivity is reported. The unchanged default is `CURRENT_OPERATOR ATOM_EMBEDDING`.

This commutator includes the projected nonlocal and SOC Hamiltonian contributions, but a commutator
of projected operators is **not**, in a finite basis, the exact projection of the continuum current.
Basis convergence is necessary. It is not the periodic Bloch velocity either; that extension needs k
derivatives and AO connection terms, including nonlocal/SOC contributions. The finite model must not
be advertised as fully periodic SOC Kubo transport.

The dense kernel tests an analytic Kramers-degenerate model, complex nonorthogonal changes mixing
spin and orbital components, origin shifts, anti-Hermiticity, invariant transition-strength
contractions, and rejection of invalid metrics and position matrices. Neon and Bi2 regression
calculations exercise actual GTH-SOC operators together with the localizer.

## Periodic Bloch Kubo current

For periodic GPW/GAPW, the distinct `CURRENT_OPERATOR BLOCH` evaluates full-mesh currents at the
frozen converged SCF potential. For example, for `CELL/PERIODIC XZ`:

```text
&KUBO_TRANSPORT ON
  CURRENT_OPERATOR BLOCH
  MP_GRID 3 1 2
  SOC T
  TEMPERATURE 300
  DISSIPATION [K] 1000
&END KUBO_TRANSPORT
```

The full Gamma-centered property mesh is independent of SCF sampling or reduction. Nonperiodic mesh
directions must have size one. A Gamma SCF is converted to image operators without changing its
potential. SOC requires a restricted SCF and SOC parameters for every kind, and is added in the
complete AO spinor space, not an occupied-state truncation. The default current operator is
unchanged; explicit SCF k-points require `BLOCH`.

Cartesian Fourier derivatives and analytic image first moments give

$$
v^a_{nm}= [C^\dagger\partial_a H C]_{nm}
-\epsilon_n[C^\dagger\partial_a S C]_{nm}
+i(\epsilon_n-\epsilon_m)[C^\dagger D_a C]_{nm}.
$$

This uses the cell-gauge form of
[Esteve-Paredes and Palacios](https://doi.org/10.21468/SciPostPhysCore.6.1.002), Eqs. (21) and (25).
The identity $D_a-D_a^\dagger=-i\partial_a S$ is checked before use. Including both derivative and
connection terms retains gauge covariance and incorporates the image Hamiltonian's nonlocal/SOC
contributions. It is not the bare Peierls derivative or the canonical momentum. Finite AO basis
errors still require convergence.

The response sums normalized full-mesh weights and projects the tensor onto the periodic lattice
subspace. It includes the Fermi divided-difference limit at equal energies, including diagonal band
velocities and entire degenerate blocks. No division by a band gap or eigenvector derivative is
needed. Only symmetric dissipative DC conductivity is evaluated, not anomalous/spin Hall response.
One-, two- and three-dimensional normalization follows the existing Kubo units.

This is a bounded dense reference implementation. `MAX_AO` limits primitive-cell scalar AOs;
`MAX_MEMORY_MB` bounds estimated replicated array storage per rank (MiB), excluding CP2K/DBCSR and
library workspace. Image operators and final currents are replicated; k-point eigensystems are
distributed over MPI ranks. This is not yet a sparse property solver; optional eigenframe symmetry
reconstruction is described below.

Matrix tests cover analytic band velocities, k-dependent nonunitary AO changes, energy and origin
shifts, the finite-system limit, connection rejection, and degenerate Drude weights. CP2K tests
include scalar and GTH-SOC Ne, Bi2, GPW/GAPW, reduced/full SCF, a Gamma-potential conversion, and
one-, two- and three-dimensional periodicity. A converged Ne/SOC comparison gives the same
conductivity for a primitive 3x1x2 mesh and its explicitly Gamma-sampled 3x1x2 supercell. The
inexpensive regression inputs have separate finite-cutoff references; equality between independently
converged SCF potentials is a convergence test, not imposed by changing references. Wilson/TRIM
property symmetry and material-converged transport benchmarks remain separate validation targets.

### Property eigenframe symmetry

The optional Kubo setting

```text
CURRENT_OPERATOR BLOCH
SYMMETRY T
SYMMETRY_BACKEND SPGLIB
```

reduces the **property** eigensystems, independently of the SCF mesh. `K290` remains the default
backend and `SYMMETRY` defaults to false. Both scalar and GTH-SOC spinor frames are supported. Only
the irreducible points are diagonalized. Full-rank eigenframes at equivalent points are
reconstructed with atom permutations, Gaussian AO rotations, cell-image Bloch phases, physical spin
rotations, and time reversal. Degenerate spaces retain the transformed source gauge; they are not
aligned by assigning individual eigenvectors.

In the cell Fourier convention used here, an atom mapped into image $L$ contributes
$\exp(-2\pi i k'\cdot L)$. Reciprocal-vector phases from an atomic Bloch convention must not be
added. CP2K's SOC block convention uses the conjugate spin representation, so the spatial SU(2) lift
is conjugated consistently. This changes the representation, not the SOC Hamiltonian.

Every reconstructed frame is checked against the target overlap and Hamiltonian. Currents are
evaluated with the target derivative/connection operators and checked against the transformed source
current, including the antiunitary action. Other operations reaching the same target may be tried,
but are accepted only after the same checks. Failure aborts; no unchecked state or silently relaxed
tolerance enters the conductivity. `EPS_SYMMETRY` controls geometry detection;
`EPS_SYMMETRY_OPERATORS` controls the matrix checks. The independent frozen potential must actually
have the requested symmetry. A reduced SCF calculation alone is not evidence of that covariance.

The final response still sums the full mesh with explicit weights and rotated tensors, not an
irreducible scalar-weight approximation. Dense target matrix work and replicated arrays remain, so
the smaller diagonalization count is not a claim of proportional time or memory savings. This path
is not yet connected to Wilson links, general TQC irreps/EBRs or localizer point-orbits.
