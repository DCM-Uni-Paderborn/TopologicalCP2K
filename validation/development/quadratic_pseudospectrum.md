# Quadratic pseudospectrum

`FORCE_EVAL/PROPERTIES/QUADRATIC_PSEUDOSPECTRUM` finds states that approximately localize
simultaneously near a position and an energy in an isolated GPW/GAPW calculation. It is independent
of the [spectral localizer](spectral_localizer.md), although both analyses share the analytic AO
position matrices and optional GTH spin-orbit operators. The quadratic gap is a localization
diagnostic, **not a topological index**. The construction follows
[Cerjan, Loring and Vides, J. Math. Phys. 64, 023501 (2023)](https://doi.org/10.1063/5.0098336). No
MATLAB or external research-script dependency is introduced.

## AO metric and observables

Let $S$ be the positive overlap, $H$ the converged AO Hamiltonian, and $X_j$ the three analytic
covariant first-moment matrices. At query $(E,\boldsymbol{x})$, define $A_0=H-ES$ and
$A_j=X_j-x_jS$, with $w_0=1$ and $w_j=\kappa$ for $j=1,2,3$. The AO problem is

$$
 Q_{\mathrm{AO}}c=\lambda Sc,\qquad
 Q_{\mathrm{AO}}=\sum_{j=0}^3 w_j^2 A_j^\dagger S^{-1} A_j,\qquad c^\dagger Sc=1.
$$

The lowest $\sqrt{\lambda}$ is the quadratic gap. Higher requested eigenpairs give the next joint
approximate states. CP2K prints each eigenpair's square root, energy residual $\delta E$, three
position residuals $\delta x_j$, expectation values, and generalized-eigenpair residual:

$$
 \delta E^2=(A_0c)^\dagger S^{-1}(A_0c),\quad
 \delta x_j^2=(A_jc)^\dagger S^{-1}(A_jc),\quad
 \lambda=\delta E^2+\kappa^2\sum_j\delta x_j^2.
$$

These formulas are invariant under nonsingular AO basis changes. Multiplying covariant AO matrices
directly without the intervening metric solution would not have this property. Likewise, translating
both the coordinate operators and the query, or shifting both the Hamiltonian and the query energy,
does not change the result.

The position residual concerns the **position operator projected into the finite AO space**. It is
not the exact continuum variance involving an independently evaluated $\langle r_j^2\rangle$. In
particular, extremely small AO spaces can underestimate localization residuals. Converge the basis,
finite-system extent and $\kappa$ before interpreting results. The complete AO space is used, not an
occupied-state projection.

## Input and states

Use `CELL PERIODIC NONE`, `POISSON PERIODIC NONE`, a suitable isolated Poisson solver and no
`KPOINTS` section. After a converged SCF calculation:

```text
&PROPERTIES
  &QUADRATIC_PSEUDOSPECTRUM ON
    SOLVER DENSE
    POSITION [angstrom] 3.0 3.0 3.0
    ENERGY [hartree] -0.1 0.0 0.1
    KAPPA [hartree*bohr^-1] 0.1
    NSTATES 2
    EPS_EIGEN 1.0E-8
    &STATES ON
      FILENAME localized
    &END STATES
  &END QUADRATIC_PSEUDOSPECTRUM
&END PROPERTIES
```

Repeat `POSITION` and supply lists of energies/scales to scan their Cartesian product. Coordinates
refer to the `SUBSYS/COORD` frame. The `.states` print key contains complex S-normalized AO
coefficients, with query and eigenvalue headers. For SOC the order is all spin-up AOs followed by
all spin-down AOs. These are eigenvectors of $Q$, not eigenvectors of $H$. Within degenerate
eigenspaces coefficients are gauge-dependent; compare projectors or complete-subspace observables.

`SPIN_CHANNEL` selects a scalar collinear Hamiltonian. `SOC T` instead requires a restricted SCF
calculation and GTH pseudopotentials containing SOC parameters. It adds CP2K's complex GTH SOC
operator after SCF, not a self-consistent noncollinear density. GPW, GAPW and all-electron GAPW are
supported for scalar analysis; the SOC path specifically requires GTH potentials.

## Dense and iterative paths

- `SOLVER DENSE` is the default finite-system reference. It gathers AO matrices and uses a full
  S-orthonormal trial space, reusable Cholesky metric factors and LAPACK diagonalization. Memory is
  quadratic, runtime cubic. It requires neither MUMPS nor Tacho.
- `SOLVER ITERATIVE` uses DBCSR block products to apply $A_jS^{-1}A_j$ without forming $Q$, $S^{-1}$
  or $S^{-1/2}$. An MPI build with [MUMPS](../../technologies/mumps.md) is required. Sparse positive
  metric factors are validated once and reused across all products and queries. No PARDISO is
  needed. A complex restarted block Rayleigh-Ritz iteration expands preconditioned residuals and
  performs double S-metric reorthogonalization. Thick restarts retain additional low Ritz vectors to
  improve convergence of clustered spectra. `MAX_SUBSPACE` bounds the replicated vector space;
  sparse AO matrices and their products remain distributed. MUMPS currently accepts centralized
  multiple RHS and broadcasts the solutions. This is not a fully distributed eigensolver.

`MAX_AO` defaults to 512 scalar functions as a safety limit for either solver; SOC doubles the
physical dimension. Raise it deliberately. Set `MAX_SUBSPACE` to at least twice `NSTATES`, unless
the whole AO/spinor space fits. `MAX_ITER` bounds iteration and `EPS_EIGEN` checks
$\|Qc-\lambda Sc\|_{S^{-1}}\leq\epsilon\max(1,|\lambda|)$ in atomic units. Unconverged solves abort
rather than presenting partial results as converged. Increase the trial space for clusters of nearby
eigenvalues and cross-check against dense results whenever feasible. A small residual is an
eigenpair check, not an independent certificate that no lower eigenvalue was missed.

`EPS_METRIC` rejects nearly dependent bases. Dense mode checks the actual overlap eigenvalue ratio;
the sparse mode uses a conservative lower bound and can therefore reject additional borderline
cases. Neither mode silently truncates basis directions.

## Validation and scope

The unit test covers an analytic complex Pauli example, nonsingular complex basis changes,
coordinate/energy-origin invariance, a localized 40-site defect, degenerate subspaces, dense versus
restarted iteration, nonconvergence, and invalid overlap metrics. MUMPS tests verify repeated
multi-RHS solutions without refactorization. Molecular regression directories compare identical
dense/iterative references for GPW, GAPW, all-electron H2 and GTH-SOC Ne using MOLOPT UZH bases.

Local validation on 2026-09-24 passed 37 MPI and 19 serial assertions across the quadratic and
existing localizer directories, using GNU Fortran 16.1, OpenMPI 5.0.9, OpenBLAS and MUMPS 5.9.1. An
independent Cholesky-based calculation from the printed H2 AO matrices reproduced all scanned gaps
within 3.2e-11 hartree and checked S-normalization and dense/iterative state-subspace overlap.
Additional Bi2 GPW/GAPW runs used TZVP-MOLOPT-GGA-GTH-q5, GTH SOC, 68 spinor coefficients, four
states and a 32-dimensional trial space. On four MPI ranks, dense and iterative gaps agreed within
6e-15 hartree and all four metric-subspace singular values agreed with unity within 1e-10. These
larger molecular runs are additional local checks, not extra mandatory regression inputs.

This is an isolated-system implementation. Periodic coordinates, k-points, symmetry-reduced
k-meshes, Berry-connection formulations and band unfolding are not enabled. Large-system and
multi-node scaling are not established by the small validation cases. The low-level Hermitian matrix
interface can accept model Hamiltonians, but other Quickstep Hamiltonian families require their own
validated AO-operator adapters.
