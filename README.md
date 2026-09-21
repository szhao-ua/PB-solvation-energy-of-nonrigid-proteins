# Protein Solvation Energy Calculations under Nonrigid Conditions

This repository contains scripts, configuration files, and documentation for calculating polar solvation energy of nonrigid proteins in the Poisson-Boltzmann theory. The primary goal of this repository is to compute the correction term in the calculation of the polar solvation energy.

## Requirements

To run the tools and scripts in this repository, ensure you have the following installed:
- **Python 3.6+**
- **NumPy** (`pip install numpy`)
- **Matplotlib**, only for the figures in `examples/` (`pip install -r requirements.txt` installs both)
<!-- - **GROMACS v5.0.5** (or compatible version for structure preparation)
- **DelPhi v8.5.0** (for solving the Poisson-Boltzmann Equation via traditional method)
- **APBS** (Adaptive Poisson-Boltzmann Solver) -->

## Repository Structure

```text
.
├── APBS.in               # APBS template configuration file
├── DelPhi.prm            # DelPhi template configuration file
├── correction.py         # Main script to compute nonrigid energy contributions (ΔE_12)
├── readin.py             # Utility script to parse atomic coordinates (.xyzqr)
├── 70set.txt             # List of the 70 PDB IDs analyzed
├── requirements.txt      # Python dependencies
├── data-TIP3P/           # Directory for solvent-equilibrated atomic data (.xyzqr files)
├── data-Xtal/            # Directory for vacuum-equilibrated atomic data (.xyzqr files)
├── examples/             # Scripts and data for the numerical examples of the paper
│   ├── example1_two_atom/   # Example 1: perturbed two-atom system (Table 1, Tables S1-S3)
│   └── example2_proteins/   # Example 2: 70 proteins (Table 2, Figs. 2, S1, S2, Tables S4-S5)
└── README.md             # This documentation file
```

---

## Usage

### 1. Structure Preparation

All protein structures analyzed in this study were sourced from the Protein Data Bank (PDB) [1]. Structure preparation and molecular dynamics (MD) minimizations were conducted in **GROMACS v5.0.5** [2] using the **AMBER99SB** force field [3], with all titratable residues maintained in their standard charged states.

To evaluate the methods under nonrigid conditions, two distinct conformational states were utilized for each protein:

- **Solvent-equilibrated structure (Explicitly solvated state)**: The proteins were solvated using TIP3P water molecules, with ions added as necessary for neutralization. Further details on this structure preparation process can be found in Chakravorty et al. [4].
- **Vacuum-equilibrated structure (Reference state)**: Crystals were protonated and subjected to energy minimization with heavy atomic restraints (using a force constant of $1 \times 10^6\,\text{kJ}\,\text{mol}^{-1}\,\text{nm}^{-2}$) to preserve the native backbone topology. Further details on this structure preparation process can be found in Chakravorty et al. [4].

### 2. Calculating Solvation Energy

#### Using DelPhi

To run tests on the 70-protein set using DelPhi, we use the following input template. Detailed descriptions of these parameters can be found in the [DelPhi User Manual (v8.5.0)](https://compbio.clemson.edu/lab/media/download/DelPhiUserManual_v8_5_0.pdf).

```text
scale=2.0
perfil=70.0
in(pdb,file="1TG0.pdb")
in(siz,file="1TG0.siz")
in(crg,file="1TG0.crg")
indi=1.0
exdi=80.0
prbrad=1.4
salt=0.15
dencut=0.759
bndcon=2
maxc=0.0001
linit=1000
nonit=800
gaussian=0
energy(s)
out(energy,file="1TG0_TIP3P_TRAD_1.dat")
```

When calculating the energy, DelPhi provides the **Corrected reaction field energy**, which is the polar solvation energy in units of kT. This is converted to kcal/mol using the conversion factor:
$$1\,\text{kT} = 0.5922\,\text{kcal/mol}$$

#### Using APBS

Similarly, to run tests using APBS, we use the following input setup (also in `APBS.in`). Detailed descriptions of these parameters can be found in the [APBS User Manual](https://apbs.readthedocs.io/en/latest/using/examples/solvation-energies.html#polar-solvation). The grid spacing of 0.5 Å on a 225³ mesh gives a cubic box of side 112 Å, which holds every solute in the set with at least 20 Å of padding. Both `elec` blocks are centred on the same molecule, so the two states share one grid.

```text
read
    mol pqr data-TIP3P/1TG0.pqr
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

For the uncorrected dual-structure calculation, the vacuum structure is read as a second molecule, and the reference block is evaluated on it with `mol 2`. `gcent mol 1` is unchanged, so both states stay on the same grid:

```text
read
    mol pqr data-TIP3P/1TG0.pqr
    mol pqr data-Xtal/1TG0.pqr
end
...
# reference state in vacuum, on the vacuum structure (mol 2) and the same grid
elec name reference
  ...
  gcent mol 1
  mol 2
  ...
end
```

The complete rigid and dual-structure inputs used for the 70 proteins are in `examples/example2_proteins/` (see [`examples/README.md`](examples/README.md)).

Running APBS using this input file provides the **Global net ELEC energy** value in units of kJ/mol. This is converted to kcal/mol using the conversion factor: $$1\,\text{kcal/mol} = 4.184\,\text{kJ/mol}$$

### 3. Correcting for Nonrigid Contributions ($\Delta E_{12}$)

To evaluate conditions under nonrigid topologies, we must correct the rigid solvation energy by adding the $\Delta E_{12}$ nonrigid energy contribution. 

#### Data Preparation
The atomic coordinate data (`.xyzqr` files) must be placed in the respective directories:
- **Solvated data**: Place inside `data-TIP3P/`
- **Vacuum data**: Place inside `data-Xtal/`

*Note: The script `readin.py` is responsible for parsing this data. If different directory names are used, update the `base_path` variable inside `readin.py`.*

#### Running the Correction Script
The `correction.py` script computes the raw vacuum energy ($E_1$), the raw in-solvent energy ($E_2$), and the final nonrigid contribution ($\Delta E_{12}$). 

To run the calculation, use the following bash command:

```bash
python correction.py <PDB_ID>
```
*(e.g., `python correction.py 1TG0`)*

The output ($\Delta E_{12}$) is automatically saved to a text file in the working directory (e.g., `1TG0_delta_E12.txt`). 

#### Final Adjusted Solvation Energy
For both DelPhi and APBS, the final corrected polar solvation energy is simply the sum of the rigid polar solvation energy (converted to kcal/mol) and the nonrigid energy contribution ($\Delta E_{12}$):
$$\Delta E = \Delta E_{\text{rigid}} (\text{in kcal/mol}) + \Delta E_{12}$$

### 4. Reproducing the Numerical Examples

The `examples/` directory contains the input structures, solver results and scripts behind every table and figure of the two numerical examples. From the committed results:

```bash
cd examples/example1_two_atom && python tables.py --check     # Table 1, Tables S1-S3
cd examples/example2_proteins && python protein_stats.py --check   # Table 2, Table S5, statistics
cd examples/example2_proteins && python figures.py                # Fig. 2, Fig. S1
```

The `run_apbs.py` script in each example regenerates the APBS results. See [`examples/README.md`](examples/README.md) for details.

---
## References
[1] Berman HM, et al. *The Protein Data Bank*. Nucleic Acids Res. 2000.  
[2] Van Der Spoel D, et al. *GROMACS: fast, flexible, and free*. J Comput Chem. 2005.  
[3] Ponder JW, Case DA. *Force fields for protein simulations*. Adv Protein Chem. 2003.  
[4] Chakravorty A, et al. *Reproducing the ensemble average polar solvation energy of a protein from a single structure: Gaussian-based smooth dielectric function for macromolecular modeling*. J Chem Theory Comput. 2018.
