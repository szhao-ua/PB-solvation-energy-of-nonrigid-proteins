#!/usr/bin/env python3
"""
APBS runs for Example 1 (the perturbed two-atom system).

Two kinds of run are made:

  rigid  both elec blocks on the solvent structure, giving the rigid polar
         solvation energy Delta E_23; C-APBS = Delta E_23 + Delta E_12.
  dual   solvated block on the solvent structure data/solvent/<name>.pqr (mol 1),
         reference block on the perturbed vacuum structure data/vacuum/<name>.pqr
         (mol 2); this is the uncorrected APBS energy.

The solvent structure is the same for every pair, so the rigid runs, which do not
depend on the perturbation, use the pair of seed 12345, Cp = 0.1.

Both blocks are centred at the origin, so the two states share one grid.

Jobs (see the manuscript, Table 1, and the supplementary material, Tables S2-S3):
  Table 1, left   rigid + dual, h = 0.8 ... 0.1, 19.2 A box, Cp = 0.1, seed 12345
  Table 1, right  dual, h = 0.1, 19.2 A box, Cp = 0.1 ... 0.5, seed 12345
  Table S2        rigid + dual, h = 0.1, box 22.4 / 25.6 / 28.8 A
  Table S3        dual, h = 0.1, 19.2 A box, all Cp, five seeds

Usage:
    python run_apbs.py [--apbs apbs] [--only TAG]

Each job runs in runs/apbs/<name>/ and is skipped if its log already holds an
energy. All energies are collected into results/apbs.csv (kcal/mol).
The h = 0.1 jobs on the 28.8 A box need about 6 GB of memory.
"""
import argparse
import os
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RUNS = os.path.join(HERE, "runs", "apbs")
RESULTS = os.path.join(HERE, "results", "apbs.csv")

SEEDS = [12345, 2026, 99999, 123456789, 987654321]
CPS = [0.10, 0.20, 0.30, 0.40, 0.50]
BOX = 19.2                                          # A, held fixed under mesh refinement
DIME = {0.8: 25, 0.6: 33, 0.4: 49, 0.2: 97, 0.1: 193}
KJ_PER_KCAL = 4.184


def job_list():
    """(kind, h, dime, seed, cp) for every run used in the paper."""
    jobs = set()
    for h, d in DIME.items():
        jobs.add(("rigid", h, d, 0, 0.0))
        jobs.add(("dual", h, d, 12345, 0.10))
    for s in SEEDS:
        for cp in CPS:
            jobs.add(("dual", 0.1, 193, s, cp))
    for d in (225, 257, 289):
        jobs.add(("rigid", 0.1, d, 0, 0.0))
        jobs.add(("dual", 0.1, d, 12345, 0.10))
    # smallest grids first
    return sorted(jobs, key=lambda j: (j[2], j[0], j[3], j[4]))


def job_name(kind, h, dime, seed, cp):
    return f"{kind}_h{h}_d{dime}" + (f"_s{seed}_cp{cp:.2f}" if kind == "dual" else "")


def elec_block(name, mol, h, dime, sdie, ions):
    ion = ("  ion charge  1 conc 0.150 radius 2.0\n"
           "  ion charge -1 conc 0.150 radius 2.0\n") if ions else ""
    return (f"elec name {name}\n"
            f"  mg-manual\n"
            f"  dime {dime} {dime} {dime}\n"
            f"  grid {h} {h} {h}\n"
            f"  gcent 0.0 0.0 0.0\n"
            f"  mol {mol}\n"
            f"  npbe\n"
            f"  bcfl sdh\n"
            f"  pdie 1.0\n"
            f"  sdie {sdie}\n"
            f"{ion}"
            f"  chgm spl0\n"
            f"  srfm mol\n"
            f"  srad 1.4\n"
            f"  swin 0.3\n"
            f"  sdens 10.0\n"
            f"  temp 298.15\n"
            f"  calcenergy total\n"
            f"  calcforce no\n"
            f"end\n")


def apbs_input(kind, h, dime):
    mols = ["solvent.pqr"] + (["vacuum.pqr"] if kind == "dual" else [])
    read = "read\n" + "".join(f"    mol pqr {m}\n" for m in mols) + "end\n\n"
    return (read
            + elec_block("solvated", 1, h, dime, "80.0", True) + "\n"
            + elec_block("reference", 2 if kind == "dual" else 1, h, dime, "1.0", False)
            + "\nprint elecEnergy solvated - reference end\n\nquit\n")


def read_energy(log):
    """Global net ELEC energy in kcal/mol, or None if the run did not finish."""
    if not os.path.exists(log):
        return None
    E = None
    with open(log, errors="replace") as fh:
        for line in fh:
            if "Global net ELEC energy" in line:
                E = float(line.split("=")[1].split()[0]) / KJ_PER_KCAL
    return E


def run(job, apbs):
    kind, h, dime, seed, cp = job
    d = os.path.join(RUNS, job_name(*job))
    log = os.path.join(d, "run.log")
    if read_energy(log) is None:
        os.makedirs(d, exist_ok=True)
        # APBS cannot read paths containing spaces, so inputs are copied next to the .in file.
        name = f"2ATO_s{seed}_cp{cp:.2f}.pqr" if kind == "dual" else "2ATO_s12345_cp0.10.pqr"
        shutil.copy(os.path.join(DATA, "solvent", name), os.path.join(d, "solvent.pqr"))
        if kind == "dual":
            shutil.copy(os.path.join(DATA, "vacuum", name), os.path.join(d, "vacuum.pqr"))
        with open(os.path.join(d, "run.in"), "w") as fh:
            fh.write(apbs_input(kind, h, dime))
        with open(log, "w") as fh:
            subprocess.run([apbs, "run.in"], cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    return read_energy(log)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apbs", default="apbs", help="APBS executable (default: apbs on PATH)")
    ap.add_argument("--only", help="run only jobs whose name contains this string, e.g. h0.6")
    args = ap.parse_args()

    rows = []
    for job in job_list():
        if args.only and args.only not in job_name(*job):
            continue
        E = run(job, args.apbs)
        print(f"{job_name(*job):32s} {E if E is not None else 'FAILED'}", flush=True)
        rows.append((job, E))

    # merge with results already on file, so a partial run does not drop other rows
    table = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as fh:
            next(fh)
            for line in fh:
                k, h, d, _box, s, cp, E = line.strip().split(",")
                table[(k, float(h), int(d), int(s), round(float(cp), 2))] = E
    for (kind, h, dime, seed, cp), E in rows:
        if E is not None:
            table[(kind, h, dime, seed, round(cp, 2))] = repr(E)
    with open(RESULTS, "w") as fh:
        fh.write("kind,h,dime,box_A,seed,cp,energy_kcal\n")
        for k in sorted(table, key=lambda k: (k[0], -k[1], k[2], k[3], k[4])):
            kind, h, dime, seed, cp = k
            fh.write(f"{kind},{h},{dime},{(dime - 1) * h:.1f},{seed},{cp:.2f},{table[k]}\n")
    print(f"wrote {RESULTS}")


if __name__ == "__main__":
    main()
