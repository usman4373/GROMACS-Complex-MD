# Protein–Ligand Molecular Dynamics Simulation with GROMACS

This guide provides a step‑by‑step protocol to set up and run a molecular dynamics (MD) simulation of a **protein–ligand complex** using GROMACS. It covers ligand topology generation, system building, solvation, ion addition, energy minimization, equilibration (NVT/NPT), production MD, and post‑processing.

## Prerequisites

- [GROMACS](https://www.gromacs.org/) (2020 or later) installed with GPU support (optional but recommended)
- Protein structure file (`protein.pdb`) and ligand structure file (`ligand.pdb` or `.mol2`)
- Parameter files (`ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`) – obtain from [standard GROMACS tutorials](http://www.mdtutorials.com/gmx/)
- For ligand topology: [CGenFF](https://cgenff.umaryland.edu/) (or other force field tools) and the `cgenff_charmm2gmx.py` script

> **Note:** The ligand topology generation is force‑field specific. This guide uses **CHARMM36** (via CGenFF) as an example. Adjust if using another force field.

---

## Step 1: Prepare topology and initial structure

### 1.1 Process the protein

```bash
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -ter -p protein_topol.top
```

During the run, you will be prompted to select:

- Force field (e.g., CHARMM36, Amber99SB, OPLS)
- Water model (e.g., TIP3P)
- Terminus type (usually none or both, depending on the protein)
- `-ter` ensures correct handling of charged termini.

### 1.2 Generate ligand topology

If your ligand is not a standard residue, pdb2gmx will fail. Follow this external workflow:

- Add hydrogens to the ligand using Avogadro and save as `.mol2`.
- Fix bond order in the `.mol2` file using a Perl script (provided in this repository).
- Generate ligand topology using CGenFF:
    - Upload your `.mol2` file → get `ligand.str` and `ligand_fix.mol2`.

- Convert CHARMM to GROMACS format:

```bash
python cgenff_charmm2gmx.py ligand ligand_fix.mol2 ligand.str charmm36-jul2022.ff
```

This produces a ligand topology file (e.g., ligand.itp) and a gro file.

### 1.3 Convert ligand PDB to GRO

```bash
gmx editconf -f ligand.pdb -o ligand.gro
```

### 1.4 Build the protein–ligand complex

- Combine the protein and ligand coordinates into a single `complex.gro`.
- You can manually edit a `.gro` file or use gmx editconf with `-cat`:

```bash
gmx editconf -f protein_processed.gro -o complex.gro -cat ligand.gro
```

Then, merge the topologies:

- Open `protein_topol.top` and `ligand.itp`.
- Include the ligand `.itp` file inside `topol.top` before the [ molecules ] directive.
- Add the ligand as a new molecule in the [ molecules ] section:

```
; Include ligand topology
#include "ligand.itp"

[ molecules ]
; name   number
Protein     1
Ligand      1
```

> Keep a backup of your topology file.

## Step 2: Solvation

### 2.1 Define the simulation box (dodecahedron recommended)

```bash
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -c -d 1.0
```

| Option     | Meaning                                                   |
|------------|-----------------------------------------------------------|
| `-bt`      | Rhombic dodecahedron box (better volume/surface ratio)    |
| `-c`       | Center the complex in the box                             |
| `-d 1.0`   | Minimum distance (nm) from solute to box edge             |

### 2.2 Fill the box with water

```bash
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

> `-cp` = solute configuration (complex in box) `-cs` = solvent configuration (SPC216 water model)

## Step 3: Add ions to neutralize the system
- First, assemble a `.tpr` file using `ions.mdp`:

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```

Then replace solvent molecules with ions:

```bash
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

> When prompted, select the SOL group (water) – this ensures ions are placed only in the solvent region.

## Step 4: Energy minimization

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -s em.tpr -deffnm em -nb gpu
```

## Step 5: Equilibration (NVT and NPT)

For a protein–ligand complex, we apply position restraints to the ligand (and optionally to the protein backbone) and use separate temperature coupling groups.

### 5.1 Create index file and ligand restraints

```bash
gmx make_ndx -f ligand.gro -o index_lig.ndx
```

Then generate restraints for the ligand:

```bash
gmx genrestr -f ligand.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
```

Now edit `topol.top` to include `#include "posre_lig.itp"` after the ligand topology include, and before [ molecules ].

### 5.2 Prepare a combined index for temperature coupling

Create an index file that merges protein and ligand into one group (for temperature coupling) or keeps them separate as needed.

```bash
gmx make_ndx -f em.gro -o index.ndx
```
Inside the interactive prompt:

- Type 1 | 13 (assuming 1 = protein, 13 = ligand) → creates group 14 (Protein_Ligand).
- Type q to quit.
- Note: The group numbers may differ. Use gmx make_ndx without arguments to list groups.

### 5.3 NVT equilibration (constant temperature)

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -v -s nvt.tpr -deffnm nvt -nb gpu -pme gpu -bonded gpu
```
> The `-r` flag uses the minimized coordinates as a reference for restraints.

### 5.4 NPT equilibration (constant pressure)

```bash
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -v -s npt.tpr -deffnm npt -nb gpu -pme gpu -bonded gpu
```

## Step 6: Production MD simulation

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_10ns.tpr
gmx mdrun -v -s md_10ns.tpr -deffnm md_10ns -nb gpu -pme gpu -bonded gpu
```

- To run `50 ns`, modify nsteps in `md.mdp` accordingly.
- Use `-deffnm` to set the base name of all output files (e.g., `md_50ns.gro`, `md_50ns.xtc`).

## Step 7: Post‑processing and analysis












