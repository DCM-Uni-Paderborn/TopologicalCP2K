Completed Linux build and finite-Bi controls, 1 October 2026

Scope
  Each host uses the same archived research source. Five unit programs
  pass at 1 MPI rank/2 OpenMP threads and at 4 ranks/1 thread. Both
  selected regression suites pass 81/81 numerical checks. These builds
  omit MUMPS and Tacho and do not establish multi-node scaling.

Library negative control
  Spark's installed OpenBLAS 0.3.26 automatic armv8sve path computes
  (1+2i,3+4i)^dagger (5+6i,7+8i) as 46-8i, instead of 70-8i.
  The ctypes CBLAS reproducer uses scalar Python references, not NumPy.
  OPENBLAS_CORETYPE=ARMV8 fixes the tested kernels and CP2K tests.
  This environment setting is scoped to the individual processes.
  No global installation or CP2K reference values were changed.
  The archived initial build script and failed logs deliberately retain
  automatic dispatch. To reproduce passing tests after building, set
  OPENBLAS_CORETYPE=ARMV8 and run qualify_build.py in a fresh output tree.

Material control
  The fresh Terok 48-atom Bi restart is compared to the retained macOS
  calculation with the identical initial restart fingerprint. The
  comparator also checks geometry, basis/potential data and numerical
  settings. Total energy is identical at printed precision. The maximum
  spinor-energy difference is 3.42e-14 hartree, and the complete-AO
  physical-Hamiltonian difference has spectral norm 1.14e-13 hartree.
  The original comparison to a different restart stage failed the
  unchanged tolerance and is retained as a rejected control.
  Two resource-limited Spark material attempts are retained as aborted,
  not presented as valid energies or scientific validation.

Reproduction
  python verify_platform_archives.py .
  This checks archive hashes, member hashes and recorded test consistency.
  It does not rerun CP2K or certify untested platforms.

  For independent numerical replay, extract terok-platform-checks.tar.gz
  into a new empty directory, install NumPy and SciPy, and run there:
  OPENBLAS_NUM_THREADS=1 python compare_platform.py reference48 \
      --reference reference48-restarted --output replay.json
  This reintegrates Gaussian overlap/position matrices and reconstructs
  both complete physical Hamiltonians from the archived SOC eigenstates.
  It requires no CP2K executable or new SCF calculation.

Provenance
  Source commit c4fa12b8bd2001b390a92fdfcd12e56a25ba999f.
  Source archive SHA256:
  13286f5dacc4f69265d333f0dd55ac4bc8ffbfbeb2937b6231ccb681056b466e
  Both hosts verify 9,787 archived regular files, executable bits, and
  symlink targets. Owner/group IDs and umasks are platform dependent.
  Individual archive indexes contain binary/library/source fingerprints.

Not established
  The larger 96-atom workflow is ongoing and is not included in these
  completed results. Finite-size, basis, cutoff, temperature and bulk
  convergence are distinct from this matched-restart comparison.
