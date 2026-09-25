"""Native CP2K Bi(111) periodic reference and isolated, unrelaxed flakes.

Geometry: Li et al., New J. Phys. 23, 063042 (2021),
doi:10.1088/1367-2630/ac04c9, section 4.1. These are restricted PBE
potentials with post-SCF GTH SOC, not self-consistent spinor DFT.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def geometry(size, vacuum):
    a, buckling = 4.33, 1.74
    height = math.sqrt(3) * a / 2
    if size == 0:
        cell = [[a, 0, 0], [a / 2, height, 0], [0, 0, vacuum]]
        atoms = [[a / 2, height / 3, vacuum / 2 + buckling / 2],
                 [a, 2 * height / 3, vacuum / 2 - buckling / 2]]
    else:
        atoms = [[a * (ix + u) + a / 2 * (iy + u), height * (iy + u), z]
                 for iy in range(size) for ix in range(size)
                 for u, z in ((1 / 3, buckling / 2), (2 / 3, -buckling / 2))]
        lo = [min(p[d] for p in atoms) for d in range(3)]
        hi = [max(p[d] for p in atoms) for d in range(3)]
        lengths = [hi[d] - lo[d] + vacuum for d in range(3)]
        cell = [[lengths[i] if i == j else 0 for j in range(3)] for i in range(3)]
        atoms = [[p[d] - lo[d] + vacuum / 2 for d in range(3)] for p in atoms]
    return cell, atoms


def nnkp(cell):
    # Finite-cell spectrum: one Gamma point, self link, no reciprocal translation.
    import numpy as np

    reciprocal = 2 * np.pi * np.linalg.inv(cell).T
    lines = ["begin real_lattice"]
    lines += [" ".join(f"{x:.16g}" for x in row) for row in cell]
    lines += ["end real_lattice", "begin recip_lattice"]
    lines += [" ".join(f"{x:.16g}" for x in row) for row in reciprocal]
    lines += ["end recip_lattice", "begin kpoints", "1", "0 0 0", "end kpoints",
              "begin nnkpts", "1", "1 1 0 0 0", "end nnkpts", ""]
    return "\n".join(lines)


def inputs(args):
    finite = all(math.isfinite(x) for x in (args.size, args.vacuum, args.temperature, args.cutoff, args.scf_mesh))
    if not finite or args.size < 0 or args.vacuum <= 0 or args.temperature < 0 or args.cutoff <= 0 or args.scf_mesh < 1:
        raise ValueError("Invalid geometry or numerical setting")
    cell, atoms = geometry(args.size, args.vacuum)
    nao = len(atoms) * {"DZVP": 13, "TZVP": 17, "TZV2P": 29}[args.basis]
    if args.mode == "wilson" and args.size != 0:
        raise ValueError("The Wilson reference must be periodic")
    if args.mode != "wilson" and args.size < 1:
        raise ValueError("Finite analysis requires a nonzero flake size")
    periodic = "XY" if args.size == 0 else "NONE"
    solver = "ANALYTIC" if args.size == 0 else "MT"
    kpoints = ""
    smear = ""
    if args.size == 0:
        kpoints = f"""    &KPOINTS
      SCHEME MONKHORST-PACK {args.scf_mesh} {args.scf_mesh} 1
      FULL_GRID T
      SYMMETRY F
      WAVEFUNCTIONS COMPLEX
    &END
"""
    elif args.temperature > 0:
        smear = f"""      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE {args.temperature}
      &END
"""
    analysis, properties = "", ""
    if args.mode == "wilson":
        analysis = f"""    &PRINT
      &WANNIER90
        KPOINTS_SOURCE WILSON
        SEED_NAME bismuth
        SOC T
        Z2 T
        TIME_REVERSAL T
        REQUIRE_GLOBAL_GAP T
        EXCLUDE_BANDS {' '.join(str(i) for i in range(11, 2 * nao + 1))}
        WILSON_MESH 12 13
        WILSON_MAX_REFINEMENT 4
        WILSON_TOL 1.0E-3
      &END
    &END
"""
    elif args.mode == "spectrum" or getattr(args, "export_spectrum", False):
        analysis = """    &PRINT
      &WANNIER90
        KPOINTS_SOURCE NNKP
        NNKP_FILE gamma.nnkp
        SEED_NAME bismuth
        SOC T
        TIME_REVERSAL T
        STATE_EXPORT T
      &END
    &END
"""
    if args.mode == "localizer":
        if not args.energy:
            raise ValueError("Localizer queries require explicit absolute energies")
        if (not all(math.isfinite(x) for x in (*args.energy, *args.kappa, *args.offset))
                or not args.kappa or min(args.kappa) <= 0 or not args.offset):
            raise ValueError("Nonfinite query or nonpositive localizer scale")
        center = [sum(p[d] for p in atoms) / len(atoms) for d in range(3)]
        positions = ["      POSITION [angstrom] " + " ".join(
            f"{center[d] + (x if d == 0 else 0):.16g}" for d in range(3))
                     for x in args.offset]
        properties = f"""  &PROPERTIES
    &SPECTRAL_LOCALIZER
      FORMULATION FINITE
      INVARIANT Z2
      SOC T
      TIME_REVERSAL T
      SOLVER {args.solver}
      MAX_AO {nao}
      PLANE XY
{chr(10).join(positions)}
      ENERGY [hartree] {' '.join(map(str, args.energy))}
      KAPPA [hartree/bohr] {' '.join(map(str, args.kappa))}
    &END
  &END
"""
    text = f"""! Fixed bare Bi(111) patch, not an edge-relaxed materials prediction.
! Geometry: doi:10.1088/1367-2630/ac04c9. All {nao} scalar AOs retained for SOC.
&GLOBAL
  PROJECT bismuth
  PRINT_LEVEL LOW
  RUN_TYPE ENERGY
&END
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME GTH_SOC_POTENTIALS
{kpoints}    &POISSON
      PERIODIC {periodic}
      POISSON_SOLVER {solver}
    &END
    &MGRID
      CUTOFF {args.cutoff}
      REL_CUTOFF 40
    &END
    &QS
      EPS_DEFAULT 1e-12
      EPS_PGF_ORB 1e-14
      EPS_PPNL 1e-14
      METHOD GPW
    &END
    &SCF
      ADDED_MOS -1
      EPS_SCF 1e-9
      MAX_SCF 250
      SCF_GUESS ATOMIC
      &DIAGONALIZATION
      &END
      &MIXING
        ALPHA 0.10
        METHOD BROYDEN_MIXING
      &END
{smear}      &PRINT
        &RESTART OFF
        &END
        &RESTART_HISTORY OFF
        &END
      &END
    &END
    &XC
      &XC_FUNCTIONAL PBE
      &END
    &END
{analysis}  &END
  &SUBSYS
    &CELL
{chr(10).join('      ' + axis + ' ' + ' '.join(f'{v:.16g}' for v in row) for axis, row in zip('ABC', cell))}
      PERIODIC {periodic}
    &END
    &COORD
{chr(10).join('      Bi ' + ' '.join(f'{v:.16g}' for v in p) for p in atoms)}
    &END
    &KIND Bi
      BASIS_SET {args.basis}-MOLOPT-GGA-GTH-q5
      POTENTIAL GTH-PBE-q5
    &END
  &END
{properties}&END
"""
    if getattr(args, "ao_matrices", False):
        if not analysis:
            raise ValueError("AO diagnostics require a spectrum export")
        text = text.replace("\n    &PRINT\n", "\n    &PRINT\n      &AO_MATRICES ON\n"
                            "        OVERLAP T\n        KOHN_SHAM_MATRIX T\n"
                            "        POSITION T\n        SOC T\n        NDIGITS 16\n        FILENAME ao\n"
                            "      &END\n", 1)
    result = {"input.inp": text}
    if getattr(args, "checkpoint", False):
        result["input.inp"] = result["input.inp"].replace("&RESTART OFF", "&RESTART ON\n"
            "          BACKUP_COPIES 1\n          &EACH\n            QS_SCF 5\n          &END")
    if getattr(args, "restart", None):
        result["input.inp"] = result["input.inp"].replace("SCF_GUESS ATOMIC", "SCF_GUESS RESTART")
        result["input.inp"] = result["input.inp"].replace("    BASIS_SET_FILE_NAME", "    WFN_RESTART_FILE_NAME initial.wfn\n    BASIS_SET_FILE_NAME")
    if args.mode == "spectrum" or getattr(args, "export_spectrum", False):
        result["gamma.nnkp"] = nnkp(cell)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=["wilson", "spectrum", "localizer"], required=True)
    parser.add_argument("--size", type=int, default=0)
    parser.add_argument("--vacuum", type=float, default=20.)
    parser.add_argument("--basis", choices=["DZVP", "TZVP", "TZV2P"], default="DZVP")
    parser.add_argument("--cutoff", type=int, default=400)
    parser.add_argument("--scf-mesh", type=int, default=8)
    parser.add_argument("--temperature", type=float, default=300.)
    parser.add_argument("--solver", choices=["DENSE", "TACHO"], default="TACHO")
    parser.add_argument("--energy", type=float, nargs="+")
    parser.add_argument("--kappa", type=float, nargs="+", default=[.001, .003, .01])
    parser.add_argument("--offset", type=float, nargs="+", default=[0.])
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--checkpoint", action="store_true")
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--export-spectrum", action="store_true")
    parser.add_argument("--ao-matrices", action="store_true")
    args = parser.parse_args()
    root, work = args.root.resolve(), args.output.resolve()
    if shutil.disk_usage(root).free < 2 * 1024**3:
        raise RuntimeError("Less than 2 GiB free; no new CP2K calculation started")
    expected = inputs(args)
    work.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, work / "runner.py")
    if args.restart:
        shutil.copyfile(args.restart.resolve(), work / "initial.wfn")
    for name, text in expected.items():
        with (work / name).open("x") as stream:
            stream.write(text)
    executable = root / ("build-mpi/bin/cp2k.psmp" if args.ranks else "build-serial/bin/cp2k.ssmp")
    command = (["mpiexec", "-n", str(args.ranks)] if args.ranks else []) + [str(executable),
              "-i", "input.inp", "-o", "output.out"]
    env = dict(os.environ, CP2K_DATA_DIR=str(root / "data"), OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="2", OMP_STACKSIZE="64M")
    stamp = dict(source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"],
                 cwd=root, text=True).strip(), executable_sha256=digest(executable),
                 runner_sha256=digest(Path(__file__)), sources={name: digest(root / name) for name in
                 ("data/BASIS_MOLOPT_UZH", "data/GTH_SOC_POTENTIALS")},
                 runtime={name: env.get(name) for name in
                          ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "DYLD_LIBRARY_PATH")})
    library = root / ("build-mpi/src/libcp2k.2026.2.dylib" if args.ranks else "build-serial/src/libcp2k.2026.2.dylib")
    stamp["library_sha256"] = digest(library)
    source_diff = subprocess.check_output(["git", "diff", "HEAD", "--", "src", "CMakeLists.txt"], cwd=root)
    if source_diff:
        (work / "source.patch").write_bytes(source_diff)
        stamp["source_diff_sha256"] = digest(work / "source.patch")
    if args.restart:
        stamp["restart_sha256"] = digest(work / "initial.wfn")
    started = time.monotonic()
    interrupted_reason = None
    with (work / "stdout.log").open("x") as stream:
        run = subprocess.Popen(command, cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT)
        with (work / "process.json").open("x") as state:
            json.dump(dict(runner_pid=os.getpid(), child_pid=run.pid, command=command), state, indent=2)
        while run.poll() is None:
            try:
                run.wait(timeout=15)
            except subprocess.TimeoutExpired:
                if shutil.disk_usage(root).free < 512 * 1024**2:
                    interrupted_reason = "Less than 512 MiB free during execution"
                    run.terminate()
                    run.wait()
    text = (work / "output.out").read_text() if (work / "output.out").exists() else ""
    record = dict(options={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  provenance=stamp, command=command, returncode=run.returncode,
                  elapsed_seconds=time.monotonic() - started,
                  interrupted_reason=interrupted_reason,
                  completed="PROGRAM ENDED AT" in text,
                  scf_converged="*** SCF run converged" in text,
                  energy=re.findall(r"ENERGY\|.*?([-\d.]+)\s*$", text, re.M),
                  files={p.name: digest(p) for p in sorted(work.iterdir()) if p.is_file()},
                  diagnostics=[line.strip() for line in text.splitlines()
                               if re.search(r"(?:TOPOLOGY|WILSON|SPECTRAL_LOCALIZER)\|", line)])
    with (work / "run.json").open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    print(json.dumps(record, indent=2), flush=True)
    if not record["completed"] or not record["scf_converged"] or run.returncode:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
