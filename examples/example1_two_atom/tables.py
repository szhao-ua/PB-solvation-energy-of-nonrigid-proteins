#!/usr/bin/env python3
"""
Tables for Example 1 (the perturbed two-atom system), built from
results/dipb_reg.csv, results/apbs.csv and the structure pairs in data/solvent/
and data/vacuum/ (files of the same name form one solvent/vacuum pair).

    Table 1   mesh refinement at Cp = 0.1 and the Cp sweep at h = 0.1
    Table S1  DIPB/REG domain study, h = 0.1
    Table S2  APBS domain study, h = 0.1
    Table S3  averages over five realizations of the perturbation, h = 0.1

The corrected energies are formed as in the paper:
    C-DIPB = Delta E_23(DIPB) + Delta E_12
    REG    = Delta E_23(REG)  + Delta E_12
    C-APBS = Delta E_23(APBS) + Delta E_12
with Delta E_12 computed by correction.py from the solvent and vacuum structures.

Usage:
    python tables.py            print the tables
    python tables.py --check    also compare every entry with the published value
"""
import argparse
import contextlib
import csv
import io
import os
import statistics as st
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
from correction import nonrigid_energy  # noqa: E402

SEEDS = [12345, 2026, 99999, 123456789, 987654321]
CPS = [0.10, 0.20, 0.30, 0.40, 0.50]
DIME = {0.8: 25, 0.6: 33, 0.4: 49, 0.2: 97, 0.1: 193}
EC2 = 332.06364       # kcal A / (mol e^2)
EPSP = 1.0

# ---------------------------------------------------------------- published values
# Table 1: DIPB, APBS, REG, C-DIPB, C-APBS
PUB_T1_LEFT = {0.8: (-60.73, -18.94, -52.83, -51.60, -101.11),
               0.6: (-79.86, 46.57, -66.48, -66.06, -90.42),
               0.4: (233.70, 212.80, -68.82, -68.28, -89.19),
               0.2: (1971.66, 1954.91, -72.42, -72.70, -89.47),
               0.1: (5715.48, 5701.88, -74.47, -74.94, -88.62)}
PUB_T1_RIGHT = {0.10: (5715.48, 5701.88, -74.47, -74.94, -88.62),
                0.20: (5618.11, 5604.51, -75.00, -75.48, -89.15),
                0.30: (3001.20, 2987.56, -75.50, -75.98, -89.65),
                0.40: (5668.06, 5654.47, -76.02, -76.50, -90.18),
                0.50: (5203.88, 5190.27, -76.51, -76.99, -90.66)}
# Table S1: C-DIPB, REG
PUB_S1 = {6: (-74.9416, -74.4646), 7: (-74.9433, -74.4654), 8: (-74.9422, -74.4650), 9: (-74.9433, -74.4663)}
# Table S2: APBS, C-APBS
PUB_S2 = {193: (5701.88, -88.617), 225: (5701.85, -88.645), 257: (5701.83, -88.658), 289: (5701.83, -88.663)}
# Table S3: mean and sd of DIPB, APBS, C-DIPB, REG, C-APBS
PUB_S3 = {0.10: ((5385.11, 869.60), (5371.50, 869.62), (-74.31, 0.58), (-73.83, 0.58), (-87.99, 0.58)),
          0.20: ((5150.60, 676.18), (5136.99, 676.19), (-74.19, 1.17), (-73.71, 1.17), (-87.86, 1.17)),
          0.30: ((4472.33, 926.01), (4458.71, 926.02), (-74.07, 1.77), (-73.59, 1.77), (-87.74, 1.77)),
          0.40: ((5168.68, 756.86), (5155.07, 756.87), (-73.95, 2.37), (-73.47, 2.37), (-87.62, 2.37)),
          0.50: ((4844.41, 578.61), (4830.80, 578.62), (-73.83, 2.98), (-73.35, 2.98), (-87.51, 2.98))}


# ---------------------------------------------------------------- inputs
def load_results():
    F, A = {}, {}
    with open(os.path.join(HERE, "results", "dipb_reg.csv")) as fh:
        for r in csv.DictReader(fh):
            F[(r["config"], float(r["h"]), int(r["edge"]), round(float(r["cp"]), 2), int(r["seed"]))] = float(r["energy_kcal"])
    with open(os.path.join(HERE, "results", "apbs.csv")) as fh:
        for r in csv.DictReader(fh):
            A[(r["kind"], float(r["h"]), int(r["dime"]), int(r["seed"]), round(float(r["cp"]), 2))] = float(r["energy_kcal"])
    return F, A


def load_xyzqr(path):
    a = np.loadtxt(path, ndmin=2)
    return a[:, 0:3].T, a[:, 3]


def dE12(seed, cp):
    """Delta E_12 (kcal/mol) for one realization, computed by correction.py."""
    name = f"2ATO_s{seed}_cp{cp:.2f}.xyzqr"
    pos, charges = load_xyzqr(os.path.join(HERE, "data", "solvent", name))
    vac_pos, _ = load_xyzqr(os.path.join(HERE, "data", "vacuum", name))
    with contextlib.redirect_stdout(io.StringIO()):
        return nonrigid_energy(pos, vac_pos, charges, EC2, EPSP)


# ---------------------------------------------------------------- output helpers
class Checker:
    def __init__(self, enabled):
        self.enabled, self.n, self.bad = enabled, 0, []

    def __call__(self, label, printed, computed, decimals):
        """Compare computed values, rounded as printed, with the published ones."""
        if not self.enabled:
            return ""
        tol = 0.5 * 10 ** -decimals + 1e-9
        self.n += len(printed)
        bad = [(p, c) for p, c in zip(printed, computed) if abs(c - p) > tol]
        if bad:
            self.bad.append(f"{label}: " + ", ".join(f"published {p} vs {c:.{decimals + 2}f}" for p, c in bad))
            return "  <-- differs from published"
        return "  ok"


def header(title):
    print("\n" + "=" * 100 + "\n" + title + "\n" + "=" * 100)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare with the published values")
    args = ap.parse_args()
    check = Checker(args.check)
    F, A = load_results()

    def row(h, cp, seed=12345, edge=8):
        e = dE12(seed, cp)
        return (F[("DIPB_dual", h, edge, cp, seed)],
                A[("dual", h, DIME[h], seed, cp)],
                F[("REG_rigid", h, edge, 0.10, 12345)] + e,
                F[("DIPB_rigid", h, edge, 0.10, 12345)] + e,
                A[("rigid", h, DIME[h], 0, 0.0)] + e)

    cols = f"{'DIPB':>10} {'APBS':>10} {'REG':>10} {'C-DIPB':>10} {'C-APBS':>10}"

    header("Table 1 (left): mesh refinement, Cp = 0.1, seed 12345")
    print(f"{'h':>6} {cols}")
    for h in DIME:
        v = row(h, 0.10)
        print(f"{h:6.1f} " + " ".join(f"{x:10.2f}" for x in v) + check(f"Table 1 h={h}", PUB_T1_LEFT[h], v, 2))

    header("Table 1 (right): varying Cp, h = 0.1, seed 12345")
    print(f"{'Cp':>6} {cols}")
    for cp in CPS:
        v = row(0.1, cp)
        print(f"{cp:6.1f} " + " ".join(f"{x:10.2f}" for x in v) + check(f"Table 1 Cp={cp}", PUB_T1_RIGHT[cp], v, 2))

    header("Table S1: DIPB/REG domain study, h = 0.1, Cp = 0.1")
    e = dE12(12345, 0.10)
    print(f"{'edge':>6} {'domain in x':>12} {'C-DIPB':>11} {'REG':>11} {'successive diff.':>17}")
    prev, cd_all, rg_all = None, [], []
    for edge in PUB_S1:
        cd = F[("DIPB_rigid", 0.1, edge, 0.10, 12345)] + e
        rg = F[("REG_rigid", 0.1, edge, 0.10, 12345)] + e
        half = int(4.2 + edge)      # the DIPB/REG domain is the atomic extent padded by the edge value
        diff = "" if prev is None else f"{cd - prev:+17.4f}"
        print(f"{edge:6d} {f'[-{half},{half}]':>12} {cd:11.4f} {rg:11.4f} {diff:>17}"
              + check(f"Table S1 edge={edge}", PUB_S1[edge], (cd, rg), 4))
        prev = cd
        cd_all.append(cd)
        rg_all.append(rg)
    print(f"spread: C-DIPB {max(cd_all) - min(cd_all):.4f}, REG {max(rg_all) - min(rg_all):.4f} kcal/mol")

    header("Table S2: APBS domain study, h = 0.1, Cp = 0.1")
    print(f"{'box (A)':>8} {'padding':>8} {'grid':>6} {'APBS':>10} {'C-APBS':>10}")
    ca_all = []
    for dime in PUB_S2:
        box = (dime - 1) * 0.1
        da = A[("dual", 0.1, dime, 12345, 0.10)]
        ca = A[("rigid", 0.1, dime, 0, 0.0)] + e
        ca_all.append(ca)
        # APBS is published to 2 decimals and C-APBS to 3
        status = (check(f"Table S2 box={box:.1f} APBS", PUB_S2[dime][:1], (da,), 2)
                  + check(f"Table S2 box={box:.1f} C-APBS", PUB_S2[dime][1:], (ca,), 3))
        print(f"{box:8.1f} {(box - 8.4) / 2:8.1f} {f'{dime}^3':>6} {da:10.2f} {ca:10.3f}{status}")
    print(f"spread: C-APBS {max(ca_all) - min(ca_all):.3f} kcal/mol")

    header("Table S3: five realizations of the perturbation, h = 0.1, mean +/- sample sd")
    print(f"{'Cp':>6} " + " ".join(f"{c:>19}" for c in ("DIPB", "APBS", "C-DIPB", "REG", "C-APBS")))
    for cp in CPS:
        es = [dE12(s, cp) for s in SEEDS]
        series = ([F[("DIPB_dual", 0.1, 8, cp, s)] for s in SEEDS],
                  [A[("dual", 0.1, 193, s, cp)] for s in SEEDS],
                  [F[("DIPB_rigid", 0.1, 8, 0.10, 12345)] + x for x in es],
                  [F[("REG_rigid", 0.1, 8, 0.10, 12345)] + x for x in es],
                  [A[("rigid", 0.1, 193, 0, 0.0)] + x for x in es])
        ms = [(st.mean(v), st.stdev(v)) for v in series]
        flat = [x for m in ms for x in m]
        pub = [x for m in PUB_S3[cp] for x in m]
        print(f"{cp:6.1f} " + " ".join(f"{m:10.2f} +/-{s:6.2f}" for m, s in ms)
              + check(f"Table S3 Cp={cp}", pub, flat, 2))

    if args.check:
        print("\n" + "=" * 100)
        print(f"{check.n} published values checked, {len(check.bad)} rows differ")
        for b in check.bad:
            print("  - " + b)
        sys.exit(1 if check.bad else 0)


if __name__ == "__main__":
    main()
