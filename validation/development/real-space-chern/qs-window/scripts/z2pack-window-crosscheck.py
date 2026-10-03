"""Use this branch's CP2K export with the existing optional Z2Pack adapter."""

import argparse
from pathlib import Path
import re
import sys

import numpy as np
import z2pack


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, "/Users/tkuehne/work/cp2k-topology-tools-phasons")
    from cp2k_z2pack import CP2KSystem

    text = (root / "tests/QS/regtest-topology/helium-nnkp.inp").read_text()
    template = work / "helium.inp"
    template.write_text(text.replace("NNKP_FILE helium.nnkp", "NNKP_FILE loop.nnkp").replace("He 0.5 0.5 0.5", "He 0.3 0.5 0.5"))
    system = CP2KSystem(
        input_file=template,
        lattice=np.eye(3) * 4,
        command=[str(args.binary.resolve())],
        workdir=work / "loops",
        num_bands=1,
        timeout=120,
        env={"CP2K_DATA_DIR": str(root / "data"), "OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "1"},
    )
    result = z2pack.line.run(system=system, line=lambda t: [t, 0, 0], iterator=[8, 16], pos_tol=1e-6)
    assert result.convergence_report["PosCheck"]
    assert abs(np.exp(2j * np.pi * result.wcc[0]) - np.exp(-0.6j * np.pi)) < 1e-6
    for directory in sorted((work / "loops").glob("loop-*")):
        output = (directory / "run.log").read_text()
        phase = float(re.search(r"Loop 1 Berry phase \[rad\]:\s+(\S+)", output)[1])
        assert abs(np.exp(1j * phase) - np.exp(2j * np.pi * result.wcc[0])) < 1e-6, phase
    print("External Z2Pack WCC:", result.wcc)
    print("Native CP2K Berry phases agree modulo 2*pi for both directed loops.")
    print("Scope: one-band periodic compatibility control, not native linked Z2Pack or a SOC material sweep.")


if __name__ == "__main__":
    main()
