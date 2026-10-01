96-atom Bi DZVP/TZVP control and bounded Pfaffian roundoff correction
=================================================================

Scope
-----
This is a fixed-geometry basis-sensitivity test of a bare finite Bi patch,
not a complete-basis extrapolation, converged bulk classification, new
material campaign or multi-node benchmark. The preceding DZVP size-control
archive is required. The source/executable of the TZVP spectrum run matches
that archive. The corrected native control additionally loads the isolated,
fingerprinted module explicitly recorded in tzvp-roundoff-state.json.

Completed checks
----------------
The TZVP scalar space has 1,632 AOs and all 3,264 post-SCF spinors are kept.
The neutral SOC gap changes from 17.8745 to 21.8655 meV relative to DZVP.
Nine independently evaluated scales share eight indices; kappa=0.002
hartree/bohr changes from zero to one. Each basis uses its own neutral
midpoint, so this is not a common-absolute-energy comparison.

Four native TZVP queries at kappa=(0.0003, 0.0015, 0.002, 0.003) have
indices (1, 1, 1, 0). Independent Gaussian reintegration of their own full
export agrees within 3.01693e-11 hartree at the unchanged 2e-8 tolerance.
The original native run failed before assigning an index. The failure and
diagnostic are preserved in failures/, with their unsuccessful run-record.json
files. They must not be counted as successful native calculations.

The correction projects a validated real skew representation explicitly
onto its antisymmetric part and includes the discarded component in the
metric gap perturbation bound. It does not relax the symmetry, Pfaffian or
gap tolerances. localizer_roundoff.F is the exact isolated module source;
the compiled shared object is deliberately not archived. Its fingerprint,
compiler command, preloading environment and unit-test outcomes are retained.

The controlled physical cross-basis comparison uses an analytic rectangular
Gaussian overlap. The occupied-space projector distance is 0.15107169 and
mean retained weight 0.99702886. Degenerate individual eigenvectors are not
matched. The common-space Hamiltonians differ by 0.00534356 hartree in
spectral norm, or 0.00203258 after removing the mean energy shift.

Archive and replay
------------------
bismuth-basis96-roundoff-evidence/index.json describes 78 files compressed
into 13 parts (504,846,313 bytes total). Both individual and concatenated
SHA-256 hashes are checked, together with all extracted member fingerprints.
The original failed workflow state is retained rather than overwritten.

Required base bundle, relative to this directory:
  ../bismuth-size96-validation/
Its concatenated archive hash is
  4fad0127a2f5bdd54819f70b8f952b1eca36aeb1ecb06c5526e5c0f069b12582
This extension's concatenated archive hash is
  2044f7413ca9543e34d10050e67e480bcd297baa2900e8de31001af635efbe12

Check the downloaded extension without unpacking large files:
  python verify_size_checks.py bismuth-basis96-roundoff-evidence

With NumPy and SciPy installed, replay from the two evidence bundles:
  OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=1 python replay_basis_roundoff.py \
    ../bismuth-size96-validation bismuth-basis96-roundoff-evidence NEW_OUTPUT

NEW_OUTPUT must not exist. The replay reintegrates Gaussian AO moments and
repeats 9 DZVP, 9 TZVP-spectrum and 4 corrected-native queries. It recomputes
the cross-basis comparison and runs 15 provenance/covariance/archive tests.
The retained replay-report.json and roundoff-archive-replay.log record zero
replay differences on Python 3.13.5, NumPy 2.4.2 and SciPy 1.17.0. The gap
replay tolerance is 1e-10 hartree; native/reference comparison is 2e-8 hartree.
This is archive-only reanalysis, not a fresh CP2K SCF calculation.

Exclusions
----------
Unchanged DZVP exports and basis/potential data are in the required base.
Build products, runtime binaries, moment caches, restart backups and unused
Wannier files are omitted. This is research-branch evidence and does not
represent an upstream CP2K release or establish publication rights for all
dependencies. Historical source notices are preserved, not reassigned.
