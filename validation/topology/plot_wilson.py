"""Regenerate the manuscript figure from retained native CP2K Wilson spectra."""

import hashlib
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
PAPER = ROOT.parent.parent


def read_case(name, expected_loops, expected_parity):
    folder = ROOT / name
    values = np.loadtxt(folder / f"{name}.wilson")
    assert values.shape == (expected_loops, 10)
    assert np.all(np.isfinite(values))
    assert np.array_equal(values[:, 0], np.arange(1, expected_loops + 1))
    assert np.all((values[:, 2:] >= 0) & (values[:, 2:] < 1))
    assert np.all((values[:, 1] > 0) & (values[:, 1] <= 1 + 1e-10))
    log = (folder / "run.log").read_text()
    assert "PROGRAM ENDED AT" in log and "*** SCF run converged" in log
    parity = re.findall(r"Converged Z2 invariant:\s+(\d+)", log)
    assert parity == [str(expected_parity)]
    total_points = int(re.findall(r"Explicit points:\s+(\d+)", log)[-1])
    assert total_points % expected_loops == 0
    # Defaults in the archived input schema: origin 0, winding (1,0,0),
    # transverse vector (0,1/2,0); the endpoint loops are both retained.
    transverse = np.linspace(0.0, 0.5, expected_loops)
    np.savetxt(
        folder / "figure-data.csv",
        np.column_stack([transverse, values[:, 1:]]),
        delimiter=",",
        header="k2,min_singular_value," + ",".join(f"center_{i}" for i in range(1, 9)),
        comments="",
        fmt="%.17g",
    )
    return transverse, values[:, 2:], {
        "loops": expected_loops,
        "points_per_loop": total_points // expected_loops,
        "z2": expected_parity,
        "minimum_link_singular_value": float(values[:, 1].min()),
    }


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.7,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharey=True)
    summary = {}
    cases = [
        ("neon", 5, 0, "(a) Neon reference", "#17688b"),
        ("stanene", 193, 1, "(b) Buckled stanene", "#9a3a50"),
    ]
    for axis, (name, loops, parity, title, color) in zip(axes, cases):
        transverse, centers, summary[name] = read_case(name, loops, parity)
        axis.scatter(
            np.repeat(transverse, centers.shape[1]), centers.ravel(),
            s=15 if name == "neon" else 2.6, color=color, linewidths=0,
        )
        axis.set(xlim=(-0.012, 0.512), ylim=(-0.025, 1.025))
        axis.set_xticks([0, 0.25, 0.5], ["0", "1/4", "1/2"])
        axis.set_yticks([0, 0.25, 0.5, 0.75, 1])
        axis.set_xlabel(r"Transverse coordinate $k_2$")
        axis.set_title(title, loc="left", pad=10)
        axis.text(0.97, 0.96, rf"$\nu={parity}$", ha="right", va="top",
                  transform=axis.transAxes, fontsize=13)
        axis.grid(axis="y", color="0.9", linewidth=0.5)
        axis.set_axisbelow(True)
        axis.tick_params(direction="out", length=3)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"Hybrid Wannier center $\bar{x}_n$")
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.20, top=0.87, wspace=0.15)
    output = PAPER / "figures"
    output.mkdir(exist_ok=True)
    fig.savefig(output / "wilson-centers.pdf", metadata={"Creator": "Matplotlib"})
    fig.savefig(output / "wilson-centers.png", dpi=240)
    plt.close(fig)
    summary["inputs_sha256"] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for name, *_ in cases
        for path in sorted((ROOT / name).iterdir())
        if path.suffix in {".inp", ".log", ".wilson", ".json"}
    }
    (ROOT / "figure-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
