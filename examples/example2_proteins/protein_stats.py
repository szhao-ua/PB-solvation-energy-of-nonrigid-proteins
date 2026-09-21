#!/usr/bin/env python3
"""
Statistics for Example 2 (the set of 70 proteins), from results/energies.csv.

    Table 2        agreement between the corrected schemes: MAE, RMSE, bias,
                   mean relative deviation, least-squares slope, Pearson r
    Table S5       energies of ten representative proteins
    main text      sign of the uncorrected and corrected energies
    supplement     distribution of Delta E_12, and its size relative to |REG|

If the coordinate files are available (by default data-TIP3P/ and data-Xtal/ of this
repository), the script also recomputes Delta E_12 from the structures with
correction.py and gives the correlation of |Delta E_12| with the structure-pair RMSD,
the number of atoms and the net charge. With --structure-sets DIR, where DIR holds
data-TIP3P/, data-GBIS/, data-Xtal/ and data-noW/, it also gives Table S4.

Columns of results/energies.csv (kcal/mol):
    dE12                         Delta E_12 = E_C(solvent structure) - E_C(vacuum structure)
    DIPB, APBS                   uncorrected dual-structure energies
    REG                          regularized diffuse-interface energy (includes Delta E_12)
    C-DIPB, C-APBS, C-DelPhi     corrected energies, Delta E_23 + Delta E_12
    dE23_DIPB, dE23_APBS, dE23_DelPhi   rigid energies Delta E_23 on the solvent structure

Usage:
    python protein_stats.py [--check] [--coords-solvent DIR] [--coords-vacuum DIR] [--structure-sets DIR]
"""
import argparse
import contextlib
import csv
import io
import math
import os
import statistics as st
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
from correction import nonrigid_energy  # noqa: E402

EC2 = 332.06364       # kcal A / (mol e^2)
EPSP = 1.0

TABLE2 = [("C-DIPB", "REG"), ("C-DelPhi", "REG"), ("C-APBS", "REG"), ("C-APBS", "C-DelPhi")]
PUB_TABLE2 = {("C-DIPB", "REG"): (41.0, 42.8, -41.0, 3.2, 1.012, 0.9999),
              ("C-DelPhi", "REG"): (348.4, 365.4, -348.4, 28.0, 1.086, 0.9911),
              ("C-APBS", "REG"): (428.8, 448.3, -428.8, 34.3, 1.105, 0.9882),
              ("C-APBS", "C-DelPhi"): (81.4, 88.1, -80.4, 5.0, 1.020, 0.9990)}
TABLE_S5_IDS = ["1TG0", "1TQG", "1VBW", "1W0N", "1IQZ", "3LZT", "3O5Q", "3PUC", "3VOR", "4A02"]
TABLE_S5_COLS = ["DIPB", "APBS", "REG", "C-DIPB", "C-APBS", "C-DelPhi"]
PUB_S4 = {"TIP3P / restrained-min. crystal": (18.4, 8.9, 96.6),
          "TIP3P / vacuum minimized": (267.2, 246.4, 122.2),
          "GBIS / restrained-min. crystal": (137.4, 126.0, 103.5),
          "GBIS / vacuum minimized": (386.2, 354.8, 171.8)}
STRUCTURE_SETS = {"TIP3P / restrained-min. crystal": ("data-TIP3P", "data-Xtal"),
                  "TIP3P / vacuum minimized": ("data-TIP3P", "data-noW"),
                  "GBIS / restrained-min. crystal": ("data-GBIS", "data-Xtal"),
                  "GBIS / vacuum minimized": ("data-GBIS", "data-noW")}


def load_energies():
    with open(os.path.join(HERE, "results", "energies.csv")) as fh:
        return {r["pdb_id"]: {k: float(v) for k, v in r.items() if k != "pdb_id"} for r in csv.DictReader(fh)}


def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


def agreement(a, b):
    """MAE, RMSE, bias, mean relative deviation (%), slope and r of a against reference b."""
    e = [p - q for p, q in zip(a, b)]
    mb, ma = st.mean(b), st.mean(a)
    slope = sum((q - mb) * (p - ma) for p, q in zip(a, b)) / sum((q - mb) ** 2 for q in b)
    return (st.mean(abs(v) for v in e),
            math.sqrt(st.mean(v * v for v in e)),
            st.mean(e),
            100 * st.mean(abs(v) / abs(q) for v, q in zip(e, b)),
            slope,
            pearson(a, b))


def load_xyzqr(path):
    a = np.loadtxt(path, ndmin=2)
    return a[:, 0:3].T, a[:, 3]


def dE12_from_files(solvent, vacuum):
    pos, q = load_xyzqr(solvent)
    vpos, _ = load_xyzqr(vacuum)
    with contextlib.redirect_stdout(io.StringIO()):
        return nonrigid_energy(pos, vpos, q, EC2, EPSP)


def header(title):
    print("\n" + "=" * 96 + "\n" + title + "\n" + "=" * 96)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare Table 2 and Table S4 with the published values")
    ap.add_argument("--coords-solvent", default=os.path.join(REPO, "data-TIP3P"))
    ap.add_argument("--coords-vacuum", default=os.path.join(REPO, "data-Xtal"))
    ap.add_argument("--structure-sets", help="directory holding data-TIP3P, data-GBIS, data-Xtal and data-noW")
    args = ap.parse_args()

    D = load_energies()
    ids = list(D)
    bad = []
    print(f"{len(ids)} proteins in results/energies.csv")

    # ------------------------------------------------------------ internal consistency
    worst = max(abs(D[k][f"C-{s}"] - D[k][f"dE23_{s}"] - D[k]["dE12"]) for k in ids for s in ("DIPB", "APBS", "DelPhi"))
    print(f"check C-X = dE23_X + dE12 for X = DIPB, APBS, DelPhi: max |difference| {worst:.1e} kcal/mol")

    # ------------------------------------------------------------ Table 2
    header("Table 2: agreement between the corrected schemes (kcal/mol)")
    print(f"{'comparison':22s} {'MAE':>7} {'RMSE':>7} {'bias':>8} {'mean rel.':>9} {'slope':>6} {'r':>7}")
    for a, b in TABLE2:
        s = agreement([D[k][a] for k in ids], [D[k][b] for k in ids])
        line = f"{a + ' vs ' + b:22s} {s[0]:7.1f} {s[1]:7.1f} {s[2]:+8.1f} {s[3]:8.1f}% {s[4]:6.3f} {s[5]:7.4f}"
        if args.check:
            rounded = (round(s[0], 1), round(s[1], 1), round(s[2], 1), round(s[3], 1), round(s[4], 3), round(s[5], 4))
            ok = all(abs(x - p) < 1e-9 for x, p in zip(rounded, PUB_TABLE2[(a, b)]))
            line += "  ok" if ok else "  <-- differs from published"
            if not ok:
                bad.append(f"Table 2 {a} vs {b}: published {PUB_TABLE2[(a, b)]}, computed {rounded}")
        print(line)

    # ------------------------------------------------------------ Table S5
    header("Table S5: energies of ten representative proteins (kcal/mol)")
    print(f"{'PDB':6s}" + "".join(f"{c:>11}" for c in TABLE_S5_COLS))
    for k in TABLE_S5_IDS:
        print(f"{k:6s}" + "".join(f"{D[k][c]:11.2f}" for c in TABLE_S5_COLS))

    # ------------------------------------------------------------ signs
    header("Sign of the energies")
    for c in ["DIPB", "APBS", "REG", "C-DIPB", "C-APBS", "C-DelPhi"]:
        pos = [k for k in ids if D[k][c] > 0]
        print(f"{c:9s} positive for {len(pos):2d} of {len(ids)} proteins  {' '.join(pos)}")
    for c in ["C-APBS", "C-DelPhi"]:
        below = sum(D[k][c] < D[k]["REG"] for k in ids)
        print(f"{c} below REG for {below} of {len(ids)} proteins")

    # ------------------------------------------------------------ Delta E_12 distribution
    header("Distribution of Delta E_12 (kcal/mol)")
    de = [D[k]["dE12"] for k in ids]
    print(f"median {st.median(de):.1f}   mean {st.mean(de):.1f}   sd {st.stdev(de):.1f}   "
          f"min {min(de):.1f}   max {max(de):.1f}")
    rel = {k: 100 * abs(D[k]["dE12"]) / abs(D[k]["REG"]) for k in ids}
    kmax = max(rel, key=rel.get)
    print(f"|Delta E_12| / |REG|: median {st.median(rel.values()):.1f}%   "
          f"90th percentile {np.percentile(list(rel.values()), 90):.1f}%   max {rel[kmax]:.1f}% ({kmax})")

    # ------------------------------------------------------------ from coordinates
    have = [k for k in ids if os.path.exists(os.path.join(args.coords_solvent, k + ".xyzqr"))
            and os.path.exists(os.path.join(args.coords_vacuum, k + ".xyzqr"))]
    header(f"From the coordinates ({len(have)} of {len(ids)} proteins found in "
           f"{os.path.relpath(args.coords_solvent)} and {os.path.relpath(args.coords_vacuum)})")
    if have:
        diff, rmsd, nat, qnet = 0.0, {}, {}, {}
        for k in have:
            fs, fv = (os.path.join(args.coords_solvent, k + ".xyzqr"), os.path.join(args.coords_vacuum, k + ".xyzqr"))
            diff = max(diff, abs(dE12_from_files(fs, fv) - D[k]["dE12"]))
            (X, q), (Y, _) = load_xyzqr(fs), load_xyzqr(fv)
            rmsd[k] = math.sqrt(((X - Y) ** 2).sum(axis=0).mean())      # no superposition
            nat[k], qnet[k] = len(q), q.sum()
        print(f"Delta E_12 recomputed by correction.py: max |difference| from energies.csv {diff:.1e} kcal/mol")
        if len(have) > 2:
            ab = [abs(D[k]["dE12"]) for k in have]
            kr = max(rmsd, key=rmsd.get)
            print(f"structure-pair RMSD (no superposition): median {st.median(rmsd.values()):.2f} A, "
                  f"max {rmsd[kr]:.2f} A ({kr})")
            print(f"r(|Delta E_12|, RMSD)       = {pearson(ab, [rmsd[k] for k in have]):.2f}")
            rest = [k for k in have if k != kr]
            print(f"  without {kr}             = {pearson([abs(D[k]['dE12']) for k in rest], [rmsd[k] for k in rest]):.2f}")
            print(f"r(|Delta E_12|, atoms)      = {pearson(ab, [nat[k] for k in have]):.2f}")
            print(f"r(|Delta E_12|, |net charge|) = {pearson(ab, [abs(qnet[k]) for k in have]):.2f}")
    else:
        print("no coordinate files found; pass --coords-solvent and --coords-vacuum")

    # ------------------------------------------------------------ Table S4
    if args.structure_sets:
        header("Table S4: Delta E_12 for four choices of the two structures (kcal/mol)")
        print(f"{'solvent / vacuum structure':34s} {'mean':>7} {'median':>7} {'sd':>7}")
        for name, (sv, vc) in STRUCTURE_SETS.items():
            v = [dE12_from_files(os.path.join(args.structure_sets, sv, k + ".xyzqr"),
                                 os.path.join(args.structure_sets, vc, k + ".xyzqr")) for k in ids]
            s = (round(st.mean(v), 1), round(st.median(v), 1), round(st.stdev(v), 1))
            line = f"{name:34s} {s[0]:7.1f} {s[1]:7.1f} {s[2]:7.1f}"
            if args.check:
                ok = all(abs(x - p) < 1e-9 for x, p in zip(s, PUB_S4[name]))
                line += "  ok" if ok else "  <-- differs from published"
                if not ok:
                    bad.append(f"Table S4 {name}: published {PUB_S4[name]}, computed {s}")
            print(line)

    if args.check:
        print("\n" + "=" * 96)
        print("all checked values agree with the published ones" if not bad else f"{len(bad)} rows differ:")
        for b in bad:
            print("  - " + b)
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
