#!/usr/bin/env python3
"""
Figures for Example 2 (the set of 70 proteins), written to figures/.

    protein_energy.pdf       Fig. 2: REG, C-DelPhi and C-APBS for each protein
    de12_histogram.pdf       Fig. S1: distribution of Delta E_12
    de12_rmsd_scatter.pdf    Fig. S2: |Delta E_12| against the structure-pair RMSD
                             (needs the .xyzqr files of both states)

Usage:
    python figures.py [--coords-solvent DIR] [--coords-vacuum DIR]
"""
import argparse
import csv
import math
import os
import statistics as st

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.ticker as ticker       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "figures")


def load_energies():
    with open(os.path.join(HERE, "results", "energies.csv")) as fh:
        return {r["pdb_id"]: {k: float(v) for k, v in r.items() if k != "pdb_id"} for r in csv.DictReader(fh)}


def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


def protein_energy(D):
    ids = list(D)
    plt.rcParams.update({"font.family": "sans-serif", "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(14, 5))
    for col, marker in zip(["REG", "C-DelPhi", "C-APBS"], ["o", "s", "^"]):
        ax.plot(range(len(ids)), [D[k][col] for k in ids], marker=marker, markersize=5, linewidth=1.6, label=col)
    ax.set_xlabel("PDB ID", fontsize=11)
    ax.set_ylabel("Energy (kcal/mol)", fontsize=11)
    ax.set_xticks(range(len(ids)))
    ax.set_xticklabels(ids, rotation=90, ha="right", fontsize=9)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.legend(title="Method", loc="lower right", fontsize=11)
    ax.grid(axis="y", color="grey", linestyle="--", linewidth=0.4, alpha=0.5)
    fig.tight_layout()
    path = os.path.join(OUT, "protein_energy.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    plt.rcdefaults()
    print(f"wrote {path}")


def serif_style():
    plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.linewidth": 0.8})


def de12_histogram(D):
    vals = [D[k]["dE12"] for k in D]
    mean, median = st.mean(vals), st.median(vals)
    serif_style()
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    width = 40.0
    lo, hi = width * (min(vals) // width), width * (max(vals) // width + 1)
    ax.hist(vals, bins=np.arange(lo, hi + width / 2, width), color="#6a8caf", edgecolor="white", linewidth=0.6)
    ax.axvline(mean, color="#b23a2f", linestyle="--", linewidth=1.2, label=f"mean = {mean:.1f}")
    ax.axvline(median, color="#2f6b3a", linestyle=":", linewidth=1.4, label=f"median = {median:.1f}")
    ax.set_xlabel(r"$\Delta E_{12}$ (kcal/mol)")
    ax.set_ylabel("number of proteins")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path = os.path.join(OUT, "de12_histogram.pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def de12_rmsd_scatter(D, solvent_dir, vacuum_dir):
    ids, rmsd, absde = [], [], []
    for k in D:
        fs, fv = os.path.join(solvent_dir, k + ".xyzqr"), os.path.join(vacuum_dir, k + ".xyzqr")
        if not (os.path.exists(fs) and os.path.exists(fv)):
            continue
        X, Y = np.loadtxt(fs, ndmin=2)[:, :3], np.loadtxt(fv, ndmin=2)[:, :3]
        ids.append(k)
        rmsd.append(math.sqrt(((X - Y) ** 2).sum(axis=1).mean()))      # no superposition
        absde.append(abs(D[k]["dE12"]))
    if len(ids) < len(D):
        print(f"de12_rmsd_scatter.pdf skipped: coordinates found for {len(ids)} of {len(D)} proteins "
              f"(pass --coords-solvent and --coords-vacuum)")
        return
    i_out = max(range(len(ids)), key=lambda i: absde[i])
    keep = [i for i in range(len(ids)) if i != i_out]
    r_all = pearson(absde, rmsd)
    r_rest = pearson([absde[i] for i in keep], [rmsd[i] for i in keep])
    serif_style()
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    ax.scatter(rmsd, absde, s=22, color="#6a8caf", edgecolor="#2f4a63", linewidth=0.5, alpha=0.85, zorder=3)
    ax.scatter([rmsd[i_out]], [absde[i_out]], s=30, color="#b23a2f", edgecolor="#6e211a", linewidth=0.6, zorder=4)
    ax.annotate(ids[i_out], (rmsd[i_out], absde[i_out]), textcoords="offset points", xytext=(-6, 8),
                fontsize=9, color="#6e211a")
    ax.set_xlabel("structure-pair RMSD (Å)")
    ax.set_ylabel(r"$|\Delta E_{12}|$ (kcal/mol)")
    ax.text(0.03, 0.95, f"$r = {r_all:.2f}$ (all {len(ids)})\n$r = {r_rest:.2f}$ (excl. {ids[i_out]})",
            transform=ax.transAxes, fontsize=9.5, va="top")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path = os.path.join(OUT, "de12_rmsd_scatter.pdf")
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coords-solvent", default=os.path.join(REPO, "data-TIP3P"))
    ap.add_argument("--coords-vacuum", default=os.path.join(REPO, "data-Xtal"))
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    D = load_energies()
    protein_energy(D)
    de12_histogram(D)
    de12_rmsd_scatter(D, args.coords_solvent, args.coords_vacuum)


if __name__ == "__main__":
    main()
