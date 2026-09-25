# Stanene: Band Reference and Finite-Torus Controls

**The localizer results below are superseded.** A subsequent full-band
matrix comparison exposed a redundant atom-block sign in the SOC torus
assembly. The archived Wilson results are unchanged, but the original
localizer scans must not be used as current material-validation results.
See `stanene-soc-validation.md` and its separate corrected archive/replay.
The original archives and summaries remain byte-for-byte unchanged.

Fresh calculations from 25 September 2026. These use true XY periodicity
for both cell and Poisson equations, not the earlier XYZ slab data.
The original development scan retained counterexamples as well as
nontrivial queries. Numerical acceptance thresholds were not relaxed.

## Original Results (Localizer Superseded)

- DZVP/200 Ry and TZVP/400 Ry use the UZH q4 bases and GTH-PBE-q4 SOC,
  restricted PBE SCF, an 8x8x1 full mesh, EPS_SCF=1e-9 and all scalar AOs.
- Both fresh Wilson surfaces converge to Z2=1 on 192 points in 193 loops.
  The sampled indirect gaps are 0.073691334871 and 0.076568514545 eV.
- Independent NumPy SVD, matrix products and eigenvalues replay all final
  overlaps. Maximum circular WCC errors are 7.33e-15 and 7.55e-15.
- All 35 localizer queries are retained in the CSV/JSON summaries, not
  just queries agreeing with the Wilson reference.
- At DZVP N=3, E=-0.1592 Ha and ETA=0.03 Ha, serial dense, MPI dense and
  MPI Tacho all give Z2=1. The dense gap lies within the MUMPS bracket.
  Serial/MPI gap difference: 2.86e-13 Ha; sparse solve residual: 1.26e-14.
- At the same energy and ETA, N=6 instead gives Z2=0 with gap 1.14e-6 Ha.
  Nearby scales and the TZVP control also show parameter dependence.

These are finite-AO solver/integration controls, **not a size-converged
nontrivial material localizer benchmark**. No occupied-band truncation,
flattening, relaxed tolerance or substituted reference value was used.
The sufficient energy/volume conditions of the periodic-localizer theorem
have not been established for these Gaussian operators. The localizer
gap is not the electronic band gap. The basis and cutoff were changed
together, not independently extrapolated.

## Archives and Replay

Each archive contains inputs, all final Wilson overlaps/eigenvalues,
native Wilson spectra, completed output logs, run metadata, basis and
potential data, operator source snapshots, CMake caches, two research
scripts and a SHA-256 member manifest. Executables are not included;
their hashes are recorded. Runtime git banners can predate incremental
rebuilds, so the matching source snapshot is retained explicitly.

- `stanene-localizer-dzvp.tar.gz`: 69 members plus manifest,
  46,318,095 bytes, SHA-256
  `3141de6b174ed9bc376c19bc7723002f5f0ba0c342f0e9db24d2513b8efa097f`.
- `stanene-localizer-tzvp.tar.gz`: 34 members plus manifest,
  46,310,467 bytes, SHA-256
  `06fb1d70da74f58ba34ece85ebad10fda249fcff0280c749051778f65784d195`.

From this directory, with Python and NumPy:

```sh
python3 replay_stanene_localizer.py stanene-localizer-dzvp.tar.gz stanene-localizer-dzvp-summary.json
python3 replay_stanene_localizer.py stanene-localizer-tzvp.tar.gz stanene-localizer-tzvp-summary.json
```

The replay checks every member hash, SCF and program completion, common
SCF energies, query energies inside the sampled gap, native/NumPy Wilson
spectra and parity, and the serial/dense/sparse gap comparison. It writes
all localizer queries to CSV as well as JSON. It does not rerun DFT or
pretend that re-parsing a printed finite index proves bulk convergence.

Full DFT reruns can use the literal archived inputs and supplied basis/
potential files with the matching CP2K implementation. The research runner
documents generation and execution, including MPI and thread settings.
These calculations are not added as costly routine CI regressions.
