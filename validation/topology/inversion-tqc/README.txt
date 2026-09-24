Native inversion representations, parity indices and atomic signatures
====================================================================

Source commit: d97839c4bc47c7d921137641d3a06086886a6056
Base commit: 15e70cf947456415d363dbaab9c9280995de5207
Development branch: topology-native-tqc (separate from the phason branch)
Date: 18 September 2026

Scope: inversion subgroup, spinful time reversal for the indicator/EBR path.
The generic character-decomposition kernel is native, but automatic irreps,
compatibility relations and EBR catalogues for all space groups are not supplied.
There is no runtime Z2Pack, Wannier90 or Python dependency for this analysis.

source.patch is the complete difference from the public CP2K base commit.
Apply it to a checkout of that base with git apply source.patch. It includes
implementation, CMake registration, unit tests, regression inputs and manual.
Published in https://github.com/cp2k/cp2k/pull/6054 (not yet merged).
Rebased PR commit: d82311e028d800d12db797e6467baf410433f8c2.
This uses mainline 251e95f, including the merged phason extension,
and combines the build/test lists. The TQC source and regression inputs
are unchanged. CMake configuration and a freshly compiled standalone
native unit suite passed after integration. The full physical regression
record here remains pinned to the implementation commit d97839c above.

Build used GNU Fortran/C/C++ 16.1.0, Open MPI 5.0.9, OpenBLAS 0.3.33,
DBCSR 2.10.0, FFTW 3.3.11 and Libxc 7.0.0. CMake Release, MPI and OpenMP enabled,
Fortran bounds checking enabled, LIBXS and LIBINT2 disabled, external DBCSR.
Representative configuration:

  cmake -S CP2K_SOURCE -B BUILD -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc-16 \
    -DCMAKE_CXX_COMPILER=g++-16 -DCMAKE_Fortran_FLAGS=-fcheck=bounds \
    -DCP2K_USE_MPI=ON -DCP2K_USE_LIBXS=OFF -DCP2K_USE_LIBINT2=OFF \
    -DCP2K_USE_FFTW3=ON -DCP2K_USE_LIBXC=ON -DCP2K_USE_DBCSR_CONFIG=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS -DBUILD_SHARED_LIBS=ON
  cmake --build BUILD --target cp2k-bin topology_symmetry_unittest -j 8

Set DBCSR_DIR/CMAKE_PREFIX_PATH to the installed dependencies when needed.
Fortran MPI_F08 was not enabled. No source change to the existing local
k-point diagonalizer is included in this extension.

regression/ contains all seven inputs and completed outputs, run with two
MPI ranks and two OpenMP threads per rank, OPENBLAS_NUM_THREADS=1. For each:

  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 CP2K_DATA_DIR=CP2K_SOURCE/data \
    mpiexec -n 2 BUILD/bin/cp2k.psmp -i INPUT -o INPUT.out

stanene-serial.out repeats stanene-tqc.inp with one rank and one OpenMP
thread. The unit output is from the built native executable. bad-center.inp
deliberately places the inversion center at (1/4,1/4,1/4) for the neon test
and must abort without printing an index. The common include is supplied.

Recheck the frozen successful results (Python >=3.11, standard library only):

  python3 verify_results.py CP2K_SOURCE

The source tree must include source.patch because its matcher definitions
are used. SHA256SUMS identifies the archived files. This reproduces checks
of recorded output; use the commands above to rerun electronic structure.

Known independent limitation: single-rank multi-thread local k-point
diagonalization can produce non-normalized states on the tested build.
The parity checks reject such states rather than reporting an invariant.
Use one OpenMP thread for single-rank runs. A preliminary sparse-access
serialization did not pass repeated tests and was discarded, not included
as a purported fix. The earlier k-point-paper threading diagnosis remains
open. The successful two-rank tests do not establish scaling performance.

The earlier Si/Al calculations, retained neon/stanene Wilson spectra and
PRL model have not been changed. These are complementary parity checks at
the same fixed material settings, not new converged material predictions.
