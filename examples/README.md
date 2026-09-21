# Numerical examples

Scripts, input structures and solver results behind the two numerical examples of the paper and its supplementary material. The committed results in each `results/` folder reproduce every table, and the `run_apbs.py` scripts regenerate the APBS results from scratch.

Install the Python dependencies from the repository root:

```bash
pip install -r requirements.txt
```

The corrected energies are formed as in the paper, $\Delta E = \Delta E_{23} + \Delta E_{12}$, where $\Delta E_{23}$ is the rigid polar solvation energy on the solvent structure and $\Delta E_{12}$ comes from `correction.py`. The uncorrected energies (DIPB, APBS) evaluate the reference state on the vacuum structure directly and are included to show that approach failing.

| Label | Meaning |
|---|---|
| DIPB, APBS | uncorrected dual-structure energy |
| C-DIPB, C-APBS, C-DelPhi | corrected energy, $\Delta E_{23} + \Delta E_{12}$ |
| REG | regularized diffuse-interface energy, $\Delta E_{23}^{\text{REG}} + \Delta E_{12}$ |

All energies are in kcal/mol.

---

## Example 1: perturbed two-atom system (`example1_two_atom/`)

Two unit charges of opposite sign at $x = \mp 2.2$ Å (radius 2 Å) form the solvent structure. The vacuum structure displaces each coordinate by $C_p(\gamma - 1/2)$ with $\gamma \sim \mathcal{U}(0,1)$, for $C_p = 0.1, \dots, 0.5$ and five random seeds.

| File | Contents |
|---|---|
| `data/solvent/2ATO_s<seed>_cp<Cp>.xyzqr` / `.pqr` | solvent structure of each pair |
| `data/vacuum/2ATO_s<seed>_cp<Cp>.xyzqr` / `.pqr` | perturbed vacuum structure of each pair |
| `make_positions.py` | regenerates `data/` with gfortran's random number generator, as the DIPB/REG code does; `--check` compares it with the committed files |
| `run_apbs.py` | all APBS runs → `results/apbs.csv` |
| `results/dipb_reg.csv` | DIPB and REG energies |
| `tables.py` | prints Table 1 and Tables S1–S3 |

A file in `data/solvent/` and the file of the same name in `data/vacuum/` form one solvent/vacuum pair, for one seed (12345, 2026, 99999, 123456789 or 987654321) and one value of $C_p$ (0.10 to 0.50). The solvent structure is the same in all 25 pairs; only the vacuum structure changes. For example, `data/solvent/2ATO_s12345_cp0.10.xyzqr` and `data/vacuum/2ATO_s12345_cp0.10.xyzqr` hold the pair used for $C_p = 0.1$ in Table 1.

```bash
cd examples/example1_two_atom
python tables.py --check        # tables from the committed results, compared with the published values
python run_apbs.py              # rerun APBS (h = 0.1 jobs on the largest box need ~6 GB of memory)
```

`results/apbs.csv` columns: `kind` (`rigid` = $\Delta E_{23}$, `dual` = uncorrected), mesh spacing `h` (Å), grid points per side `dime`, box side `box_A` (Å), `seed`, `cp`, `energy_kcal`.

`run_apbs.py` writes one APBS input per run. The input below is for the uncorrected (dual) run at $h = 0.1$ Å on the 19.2 Å box. Both blocks are centred at the origin, so the two states share one grid. The solvated state carries 0.15 M 1:1 salt. The other runs change only `dime` and `grid` (box side = (`dime` − 1) × `grid`). The rigid runs read `solvent.pqr` alone and use `mol 1` in the reference block.

```text
read
    mol pqr solvent.pqr
    mol pqr vacuum.pqr
end

elec name solvated
  mg-manual
  dime 193 193 193
  grid 0.1 0.1 0.1
  gcent 0.0 0.0 0.0
  mol 1
  npbe
  bcfl sdh
  pdie 1.0
  sdie 80.0
  ion charge  1 conc 0.150 radius 2.0
  ion charge -1 conc 0.150 radius 2.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

elec name reference
  mg-manual
  dime 193 193 193
  grid 0.1 0.1 0.1
  gcent 0.0 0.0 0.0
  mol 2
  npbe
  bcfl sdh
  pdie 1.0
  sdie 1.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

print elecEnergy solvated - reference end

quit
```

`results/dipb_reg.csv` columns: `config`, `h` (Å), `edge` (padding around the atoms in Å), `seed`, `cp`, `energy_kcal`, where `config` is

- `DIPB_dual`: uncorrected dual-structure DIPB energy
- `DIPB_rigid`: rigid DIPB energy $\Delta E_{23}$
- `REG_rigid`: rigid REG energy $\Delta E_{23}$

These energies come from the DIPB/REG code, which is not part of this repository.

---

## Example 2: 70 proteins (`example2_proteins/`)

The proteins are listed in `70set.txt` at the repository root. The solvent structure is the TIP3P-minimized structure (`data-TIP3P/`), and the vacuum structure is the restrained-minimized crystal structure (`data-Xtal/`); see the supplementary material for the preparation.

| File | Contents |
|---|---|
| `results/energies.csv` | energies of all 70 proteins, one row per protein (columns below) |
| `apbs_rigid.in`, `apbs_dual.in` | APBS inputs for the rigid ($\Delta E_{23}$) and the uncorrected dual-structure runs, shown below |
| `run_apbs.py` | runs APBS on every protein → `results/apbs_runs.csv` |
| `protein_stats.py` | Table 2, Table S5, sign counts, distribution of $\Delta E_{12}$; with coordinates also the RMSD correlation and Table S4 |
| `figures.py` | Fig. 2, Fig. S1 and Fig. S2 → `figures/` |

```bash
cd examples/example2_proteins
python protein_stats.py --check
python figures.py
python run_apbs.py --solvent-dir /path/to/data-TIP3P --vacuum-dir /path/to/data-Xtal
```

APBS inputs. The grid spacing of 0.5 Å on a 225³ mesh gives a cubic box of side 112 Å, which holds every solute in the set with at least 20 Å of padding. Both `elec` blocks are centred on the solvent structure (`gcent mol 1`), so the two states are discretized on one grid. The solvated state carries 0.15 M 1:1 salt, as in DelPhi. `run_apbs.py` copies each protein's structures to `solvent.pqr` and `vacuum.pqr` next to the input.

`apbs_rigid.in`, the rigid polar solvation energy $\Delta E_{23}$ (C-APBS = $\Delta E_{23} + \Delta E_{12}$):

```text
read
    mol pqr solvent.pqr
end

# solvated state
elec name solvated
  mg-manual
  dime 225 225 225
  grid 0.50 0.50 0.50
  gcent mol 1
  mol 1
  npbe
  bcfl sdh
  pdie 1.0
  sdie 80.0
  ion charge  1 conc 0.150 radius 2.0
  ion charge -1 conc 0.150 radius 2.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

# reference state in vacuum
elec name reference
  mg-manual
  dime 225 225 225
  grid 0.50 0.50 0.50
  gcent mol 1
  mol 1
  npbe
  bcfl sdh
  pdie 1.0
  sdie 1.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

print elecEnergy solvated - reference end

quit
```

`apbs_dual.in`, the uncorrected dual-structure energy (APBS). This file differs from `apbs_rigid.in` only in reading the vacuum structure as molecule 2 and evaluating the reference state on it with `mol 2`. `gcent` stays on molecule 1:

```text
read
    mol pqr solvent.pqr
    mol pqr vacuum.pqr
end

# solvated state
elec name solvated
  mg-manual
  dime 225 225 225
  grid 0.50 0.50 0.50
  gcent mol 1
  mol 1
  npbe
  bcfl sdh
  pdie 1.0
  sdie 80.0
  ion charge  1 conc 0.150 radius 2.0
  ion charge -1 conc 0.150 radius 2.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

# reference state in vacuum, on the vacuum structure (mol 2) and the same grid
elec name reference
  mg-manual
  dime 225 225 225
  grid 0.50 0.50 0.50
  gcent mol 1
  mol 2
  npbe
  bcfl sdh
  pdie 1.0
  sdie 1.0
  chgm spl0
  srfm mol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy total
  calcforce no
end

print elecEnergy solvated - reference end

quit
```

Columns of `results/energies.csv`:

| Column | Meaning |
|---|---|
| `dE12` | $\Delta E_{12}$ from `correction.py` |
| `DIPB`, `APBS` | uncorrected dual-structure energies |
| `REG` | regularized diffuse-interface energy |
| `C-DIPB`, `C-APBS`, `C-DelPhi` | corrected energies |
| `dE23_DIPB`, `dE23_APBS`, `dE23_DelPhi` | rigid energies $\Delta E_{23}$ on the solvent structure |

The APBS columns are the APBS output in kJ/mol divided by 4.184, and the DelPhi column is the corrected reaction-field energy in kT times 0.5922. DelPhi was run with `DelPhi.prm` at the repository root.

The atomic coordinates of the proteins are not part of this repository (apart from the 1TG0 example at the root). The RMSD correlation (Fig. S2) and Table S4 need the `.xyzqr` files of all 70 proteins, which you pass with `--coords-solvent`/`--coords-vacuum` or `--structure-sets` (a directory containing `data-TIP3P/`, `data-GBIS/`, `data-Xtal/` and `data-noW/`).
