#!/usr/bin/env python3
"""
APBS runs for Example 2 (the set of 70 proteins).

  rigid  apbs_rigid.in: both states on the solvent structure, giving the rigid
         polar solvation energy Delta E_23 (column dE23_APBS of results/energies.csv);
         C-APBS = Delta E_23 + Delta E_12.
  dual   apbs_dual.in: solvated state on the solvent structure (mol 1), reference
         state on the vacuum structure (mol 2), on one common grid; this is the
         uncorrected APBS energy (column APBS of results/energies.csv).

Both use a 225^3 grid at 0.5 A (a 112 A box), which holds every solute in the set
with at least 20 A of padding. Each run needs about 3 GB of memory.

Usage:
    python run_apbs.py [--solvent-dir ../../data-TIP3P] [--vacuum-dir ../../data-Xtal]
                       [--list ../../70set.txt] [--kind rigid|dual|both] [--apbs apbs]

Each run is made in runs/<PDB>_<kind>/ and skipped if its log already holds an
energy. Energies are collected into results/apbs_runs.csv, in kJ/mol as APBS
prints them and in kcal/mol (1 kcal/mol = 4.184 kJ/mol).
"""
import argparse
import os
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
RUNS = os.path.join(HERE, "runs")
RESULTS = os.path.join(HERE, "results", "apbs_runs.csv")
KJ_PER_KCAL = 4.184


def read_ids(path):
    """PDB IDs from a list with one '<index> <PDB ID>' pair per line."""
    with open(path) as fh:
        return [line.split()[-1] for line in fh if line.strip()]


def read_energy(log):
    """Global net ELEC energy in kJ/mol, or None if the run did not finish."""
    if not os.path.exists(log):
        return None
    E = None
    with open(log, errors="replace") as fh:
        for line in fh:
            if "Global net ELEC energy" in line:
                E = float(line.split("=")[1].split()[0])
    return E


def run(pdb, kind, args):
    d = os.path.join(RUNS, f"{pdb}_{kind}")
    log = os.path.join(d, "run.log")
    if read_energy(log) is None:
        os.makedirs(d, exist_ok=True)
        # APBS cannot read paths containing spaces, so inputs are copied next to the .in file.
        shutil.copy(os.path.join(args.solvent_dir, f"{pdb}.pqr"), os.path.join(d, "solvent.pqr"))
        if kind == "dual":
            shutil.copy(os.path.join(args.vacuum_dir, f"{pdb}.pqr"), os.path.join(d, "vacuum.pqr"))
        shutil.copy(os.path.join(HERE, f"apbs_{kind}.in"), os.path.join(d, "run.in"))
        with open(log, "w") as fh:
            subprocess.run([args.apbs, "run.in"], cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    return read_energy(log)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solvent-dir", default=os.path.join(REPO, "data-TIP3P"), help="solvent-state .pqr files")
    ap.add_argument("--vacuum-dir", default=os.path.join(REPO, "data-Xtal"), help="vacuum-state .pqr files")
    ap.add_argument("--list", default=os.path.join(REPO, "70set.txt"), help="list of PDB IDs")
    ap.add_argument("--kind", choices=["rigid", "dual", "both"], default="both")
    ap.add_argument("--apbs", default="apbs", help="APBS executable (default: apbs on PATH)")
    args = ap.parse_args()

    kinds = ["rigid", "dual"] if args.kind == "both" else [args.kind]
    table = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as fh:
            next(fh)
            for line in fh:
                pdb, kind, kj, _kcal = line.strip().split(",")
                table[(pdb, kind)] = float(kj)

    for pdb in read_ids(args.list):
        if not os.path.exists(os.path.join(args.solvent_dir, f"{pdb}.pqr")):
            print(f"{pdb}: no structure in {args.solvent_dir}, skipped")
            continue
        for kind in kinds:
            E = run(pdb, kind, args)
            print(f"{pdb} {kind:5s} " + (f"{E:.4f} kJ/mol = {E / KJ_PER_KCAL:.4f} kcal/mol" if E is not None else "FAILED"),
                  flush=True)
            if E is not None:
                table[(pdb, kind)] = E

    with open(RESULTS, "w") as fh:
        fh.write("pdb_id,kind,energy_kJ,energy_kcal\n")
        for (pdb, kind), E in sorted(table.items()):
            fh.write(f"{pdb},{kind},{E!r},{E / KJ_PER_KCAL!r}\n")
    print(f"wrote {RESULTS}")


if __name__ == "__main__":
    main()
