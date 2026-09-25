# Separated Stanene Numerical Controls

These are fresh restricted scalar SCF calculations followed by GTH SOC in
the complete scalar band space. They are not self-consistent noncollinear
SOC DFT. The fixed two-atom stanene geometry, PBE and q4 potential are
unchanged. UZH MOLOPT DZVP and TZVP retain 26 and 34 scalar AOs, respectively.

## Controlled Changes

- At a fixed 12x12 analysis torus and 8x8 scalar SCF mesh: DZVP at
  200/400/600 Ry, and TZVP at 400/600 Ry. This separates basis from cutoff.
- With TZVP/600 Ry: SCF meshes 8x8, 12x12 and 18x18, independently of
  the 12x12 property mesh.
- With TZVP/600 Ry and 12x12 SCF: relative cutoff 40 versus 60 Ry, and
  open-axis cell height 20 versus 25 Angstrom. The Cartesian buckling
  remains 0.852 Angstrom; the fractional coordinates change accordingly.
- A separate DZVP/200 Ry, 8x8-SCF export enlarges the analysis torus to
  21x21. It is a volume control, not another SCF-integration refinement.
- A fresh TZVP/600 Ry, 12x12-SCF Wilson surface supplies an independent
  converged band-invariant reference.

All complete-band localizer queries retain the original Gaussian coordinate
links and all empty states. The energy is the midpoint of the sampled global
electronic gap; sign flattening gives +/-1 Ha. Each case uses eta/Delta
0.75, 1 and 1.25. A changed midpoint within the same electronic gap does not
change this flattened operator. Electronic gaps and rescaled localizer gaps
are distinct observables and must not be compared as the same energy gap.

## Evidence and Resources

Each case archive retains input, output, commands, source/binary/data hashes,
the complete state and neighbour exports (or Wilson surface), and every
localizer diagnostic. Existing archives are reused only after verifying
their complete content against the raw sources. The method archive includes
the helper scripts, basis/potential files, source snapshots, build cache,
Python/NumPy/SciPy versions, and hashes of loaded native libraries.

CP2K exports run with two MPI ranks and two OpenMP threads per rank. Sparse
reference analyses run sequentially in fresh processes with one BLAS thread.
Their recorded peak RSS is the whole analysis-process high-water mark, not
MPI-summed SCF memory or the Pfaffian numerical-buffer allocation count.
Numerical Pfaffian factors remain serial. These runs are not a scaling study.

Independent near-zero Ritz values are tested against the original skew
matrix. They remain estimates, not certified shifted-inertia gap brackets.
No tolerance is relaxed, small element dropped, or coordinate link
unitarized to obtain an expected integer. Native Pfaffian solve probes keep
the existing 1e-10 threshold. The shared-library analysis is not a rerun of
the entire large-torus native real-cell property path.

## Reproduction

`run_stanene_controls.py` regenerates the complete matrix on the configured
macOS CP2K research build. It requires the four retained helper scripts in
`build-serial`, the matching data, and the MPI executable and shared library.
The process-only local OpenBLAS override recorded in each run must be set;
no shell startup file is modified. `--cases` selects individual controls.
Existing calculations are accepted only when all inputs and provenance agree.
Insufficient disk space prevents starting another case.

The platform-independent archive check uses Python:

```sh
python validation/development/replay_stanene_controls.py \
  validation/development/stanene-separated-controls /tmp/stanene-controls-replay.json
```

With NumPy and SciPy installed, append
`--library /path/to/libcp2k --recompute dzvp-400-scf8` to repeat
three complete-band gap/Pfaffian calculations through a native Tacho-enabled
library. The replay verifies every archive member, fresh-versus-prior controls,
constant buckling, gap/residual acceptance and the converged Wilson output.
It also regenerates `table.tex` and `comparisons.json` byte for byte from
`summary.json`. The table is included directly by the Supporting Information.

The existing small-volume counterexamples remain in the preceding archives.
These numerical controls do not establish a complete-basis limit, asymptotic
localizer-gap convergence, the finite-range theorem's hypotheses for
flattened Gaussian operators, or a converged large SOC flake classification.

## Results

All 30 complete-band localizer queries give Z2=1, with native Pfaffians
resolved in the first external ordering. The new Wilson surface converges
to Z2=1 after four refinements; its two coarsest candidates were zero.
It uses 193 loops with 192 points each and has a sampled electronic gap
of 76.5410 meV. The localizer and Wilson comparisons are independent
analyses of the same underlying post-SCF band Hamiltonian.

At eta/Delta=1, DZVP cutoff variation changes the localizer gap by less
than 4.4e-7 Ha, whereas the common-cutoff DZVP/TZVP difference is about
1.94e-3 Ha. The TZVP/600 Ry scalar SCF refinement from 12 to 18 changes
it by 6.91e-6 Ha; the relative-cutoff and cell-height controls change it
by 1.20e-9 and -2.24e-8 Ha, respectively.

The 21x21 localizer has order 45,864 and gaps 0.331864624081,
0.282051686763 and 0.240476515333 Ha at the three scales. The eta/Delta=1
gap is still 5.2% above the earlier 18x18 result. Peak analysis RSS is
9.05 GiB. The index is consistent across these tested controls; the
gap is not demonstrated to be asymptotically size-converged.

The retained replay verifies 126 archive members and the generated table.
Fresh baseline gaps agree with the earlier native-library controls within
1.22e-13 Ha. A numerical replay of all three DZVP/400 Ry queries from the
archive reproduces the gaps exactly in the retained runtime and agrees in
parity. Across the 30 new queries, the largest independent-solve and
eigenpair residuals are 2.51e-12 and 2.50e-11, respectively.
