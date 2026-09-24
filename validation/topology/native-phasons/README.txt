Native Chern and cross-geometry validation, 18 September 2026

Scope
-----
This is a separate validation set for the added native Chern/phason methods.
The manuscript's Neon/Stanene inputs, spectra, figure and reported values
are unchanged. This archive does not relabel those earlier calculations
as results from the newer executable.

Sources
-------
Native CP2K revision: b717751cccfb8a3d98c9cfd0a40df65c3726b581
Base revision: 15e70cf (after merging the earlier Wilson/C1/state-export work)
Source: https://github.com/DCM-Uni-Paderborn/cp2k/tree/b717751cccfb8a3d98c9cfd0a40df65c3726b581
Upstream proposal: https://github.com/cp2k/cp2k/pull/6052
Companion tools revision: 812a517f76f5d1104e2300dd3e3ca53f406c4fb2
Companion version: 0.1.0.dev0
Optional adapter: Z2Pack 2.2.1

topology-tools-812a517.tar.gz contains the committed numerical modules,
tests, validators, examples, user documentation and GPL license. Internal
status notes and build products are omitted. The native executable needs
neither these Python tools nor Z2Pack nor Wannier90 at runtime.

Evidence
--------
overlap-and-pump-results.json: full-precision physical overlap differences,
  C1 pumps at 4x4, 8x8 and 12x12, the separable C2 pump at 4^4, and diagnostics.
pump-refinement.log: native physical 4x4 to 8x8 check at tolerance 1e-8.
regtests.log: four existing topology regression inputs, MPI2/OpenMP2.
tests-without-z2pack.log: 34 passed, 5 optional comparisons skipped.
tests-with-z2pack.log: 39 passed.
physical-state-fixtures.tar.gz: retained CP2K version-1 state snapshots and
  associated input/output records used by validate_native_snapshots.py.
  The snapshots were produced in the earlier state-export validation work;
  b717751 identifies the native evaluator used for the reported comparison.
  Historical absolute paths in raw logs are provenance, not required paths.

The native unit tests topology_wilson_unittest, topology_curvature_unittest
and topology_snapshot_unittest are part of the pinned CP2K source. They
passed with two OpenMP threads. The curvature test includes the occupied
doublet of a 4D Dirac model at mass -3, with mesh refinement from 6^4 to 8^4.
The 8^4 value -0.742412846261164 is deliberately retained as a finite-mesh
regression value, not rounded or presented as a converged integer.

Reproduction
------------
Build the pinned CP2K source using the normal CP2K CMake dependencies and
the targets topology_phasons, topology_wilson_unittest,
topology_curvature_unittest, and topology_snapshot_unittest. The retained
native build used GNU Fortran 16, OpenBLAS and the MPI/OpenMP variant.

In an empty scratch directory, extract both archives. The following
commands assume the resulting paths topology-tools-812a517/ and fixtures/
and a CP2K build at /absolute/path/to/build. Create fixtures/ before
extracting physical-state-fixtures.tar.gz into it with tar -xzf ... -C fixtures.

python3 -m venv validation-env
validation-env/bin/python -m pip install numpy scipy pytest
cd topology-tools-812a517
CP2K_SOURCE_DIR=/absolute/path/to/pinned-cp2k-source \
  CP2K_PHASON_BINARY=/absolute/path/to/build/bin/topology_phasons.psmp \
  ../validation-env/bin/python -m pytest -q

No-Z2Pack run: expect 34 passed and 5 optional tests skipped.
Then install z2pack==2.2.1 into the same environment and rerun for 39 passed.
The independent Wilson-kernel tests compile a small Fortran driver, so a
Fortran compiler and LAPACK/BLAS are required. FC and WILSON_LAPACK_FLAGS
can override gfortran and '-llapack -lblas', respectively. The macOS run
used WILSON_LAPACK_FLAGS='-L/opt/homebrew/opt/openblas/lib -lopenblas'.
To recompute the physical comparison from the supplied snapshots:

../validation-env/bin/python validate_native_snapshots.py \
  /absolute/path/to/build/bin/topology_phasons.psmp ../fixtures ../recomputed

The output directory must not already exist. No new SCF calculation is
needed for this comparison. Mesh manifests are regenerated with the local
fixture paths. The native auxiliary program runs as one process, optionally
with OpenMP, and is not an MPI-distributed parameter-mesh calculation.

Native definitions and mesh format are documented at:
https://github.com/DCM-Uni-Paderborn/cp2k/blob/b717751cccfb8a3d98c9cfd0a40df65c3726b581/docs/methods/electronic_structure/phason-topology.md
