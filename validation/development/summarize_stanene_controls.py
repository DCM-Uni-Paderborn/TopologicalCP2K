"""Generate the controlled-comparison table directly from archived results."""

import argparse
import json
from pathlib import Path


EV_PER_HARTREE = 27.211386245988


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    args = parser.parse_args()
    summary = json.loads((args.bundle / "summary.json").read_text())
    data = {entry["case"]["label"]: entry for entry in summary}
    pairs = [("DZVP cutoff 200 to 400", "dzvp-200-scf8", "dzvp-400-scf8"),
             ("DZVP cutoff 400 to 600", "dzvp-400-scf8", "dzvp-600-scf8"),
             ("TZVP cutoff 400 to 600", "tzvp-400-scf8", "tzvp-600-scf8"),
             ("Basis at 400 Ry", "dzvp-400-scf8", "tzvp-400-scf8"),
             ("Basis at 600 Ry", "dzvp-600-scf8", "tzvp-600-scf8"),
             ("SCF mesh 8 to 12", "tzvp-600-scf8", "tzvp-600-scf12"),
             ("SCF mesh 12 to 18", "tzvp-600-scf12", "tzvp-600-scf18"),
             ("Relative cutoff 40 to 60", "tzvp-600-scf12", "tzvp-600-rel60"),
             ("Open-axis height 20 to 25", "tzvp-600-scf12", "tzvp-600-height25")]
    comparisons = []
    for name, before, after in pairs:
        old, new = (data[label]["analysis"] for label in (before, after))
        changes = [b["gap"] - a["gap"] for a, b in zip(old["queries"], new["queries"], strict=True)]
        comparisons.append(dict(control=name, before=before, after=after,
            electronic_gap_change_ev=(new["sampled_indirect_gap_hartree"] - old["sampled_indirect_gap_hartree"]) * EV_PER_HARTREE,
            gap_changes_hartree=changes, largest_gap_change_hartree=max(map(abs, changes)),
            unchanged_indices=all(a["nu"] == b["nu"] for a, b in zip(old["queries"], new["queries"], strict=True))))
    lines = [r"\begin{table}[htbp]", r"\centering\small",
             r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Basis/cutoff/SCF & $E_g$ (meV) & $g_{0.75}$ & $g_1$ & $g_{1.25}$\\", r"\midrule"]
    for entry in summary:
        c = entry["case"]
        if c["mode"] != "bloch" or c["size"] != 12:
            continue
        label = f"{c['basis']}/{c['cutoff']}/{c['scf_mesh']}"
        if c["rel_cutoff"] != 40:
            label += r", $R=60$"
        if c["height"] != 20:
            label += r", $h=25$"
        a = entry["analysis"]
        eg = a["sampled_indirect_gap_hartree"] * EV_PER_HARTREE * 1000
        gaps = " & ".join(f"{q['gap']:.6f}" for q in a["queries"])
        lines.append(f"{label} & {eg:.6f} & {gaps}" + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}",
              r"\caption{Separately varied numerical controls for complete-band stanene",
              r"on a fixed $12\times12$ analysis torus. Cutoffs are in Ry; the final",
              r"number specifies the square scalar SCF mesh. Unless marked otherwise,",
              r"the relative cutoff $R=40$ Ry and open-axis cell height $h=20$ \AA.",
              r"The buckling is fixed at $0.852$ \AA. $E_g$ is the sampled electronic",
              r"indirect gap; $g_r$ is the localizer gap in hartree at",
              r"$\eta/\Delta=r$, with $\Delta=1$ hartree. All listed indices are",
              r"$\nu=1$. The rescaled localizer gaps are not electronic band gaps.}",
              r"\label{tab:si-stanene-separated}", r"\end{table}", ""]
    assert all(q["nu"] == 1 for entry in summary if "analysis" in entry for q in entry["analysis"]["queries"])
    for filename, content in (("table.tex", "\n".join(lines)),
                              ("comparisons.json", json.dumps(comparisons, indent=2) + "\n")):
        with (args.bundle / filename).open("x") as stream:
            stream.write(content)
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
