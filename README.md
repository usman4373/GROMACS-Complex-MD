# Protein-Ligand Molecular Dynamics Simulation with GROMACS

This guide provides a step‑by‑step protocol for setting up and running a molecular dynamics (MD) simulation of a **protein–ligand complex** using GROMACS. It covers ligand topology generation, system building, solvation, ion addition, energy minimization, equilibration (NVT/NPT), production MD, and post‑processing.

## Prerequisites

- [GROMACS](https://www.gromacs.org/) (2020 or later) installed with GPU support (optional but recommended)
- Protein structure file (`protein.pdb`) and ligand structure file (`ligand.pdb` or `.mol2`)
- Parameter files (`ions.mdp`, `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`) – obtain from [standard GROMACS tutorials](http://www.mdtutorials.com/gmx/)
- For ligand topology: [CGenFF](https://cgenff.com/) (or other force field tools) and the `cgenff_charmm2gmx.py` script

> **Note:** The ligand topology generation is force‑field specific. This guide uses **CHARMM36** (via CGenFF) as an example. Adjust if using another force field.

---

## Step 1: Prepare topology and initial structure

- If the ligand is not a standard residue, `pdb2gmx` will fail.
- Therefore, write topologies of both protein and ligand separately
- Create separate `protein.pdb` and `ligand.pdb` files from the `complex.pdb` file

### 1.1 Process the protein

### Download CHARMM36 Force field (Optional)
- Download CHARMM36 force field from [MacKerell lab website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs) (Extract the tarball in the working directrory)
- Also download `cgenff_charmm2gmx_py3_nx2.py` conversion script from the same website
- Now, generate the topology of the protein:

```bash
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p protein_topol.top
```

During the run, you will be prompted to select:
- Force field (e.g., CHARMM36, Amber99SB, OPLS)
- Water model (e.g., TIP3P)

> If writing the topology of a protein fails due to hydrogen atoms mismatch, then use `-ignh` in the above command. It will let GROMACS rebuild hydrogens correctly

### 1.2 Generate ligand topology
- Add hydrogens to the ligand using [Avogadro](https://avogadro.cc/) and save as `.mol2`.
- Fix bond order in the `.mol2` file using a [Perl script](http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.txt):

```bash
perl sort_mol2_bonds.pl ligand.mol2 ligand_fix.mol2
```

- Generate ligand topology using [CGenFF](https://cgenff.com/) (Create account to use CGenFF):
    - Upload your `ligand_fix.mol2` file → get `ligand_fix.str` and `ligand_fix.cgenff.mol2`.

- Convert CHARMM to GROMACS format:
- Use `cgenff_charmm2gmx_py3_nx2.py` script downloaded from [MacKerell lab website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs)

```bash
python cgenff_charmm2gmx_py3_nx2.py ligand ligand_fix.cgenff.mol2 ligand_fix.str charmm36-feb2026_cgenff-5.0.ff
```

This produces a ligand topology file (e.g., `ligand_fix.itp`) and a `.gro` file.

> **Tip:** If the above command does not run due to compatibility issues, then run this:
> ```bash
> conda create -n cgenff python=3.7 -y
> conda activate cgenff
> pip install networkx==2.3
> ```

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
<p align="justify">After the production MD run, the trajectory may contain artifacts due to periodic boundary conditions (PBC). Molecules can diffuse across box boundaries, making them appear “broken” or “jumping”. The first step is to correct this by centering the protein and removing PBC jumps.</p>

### 7.1 Remove periodic boundary effects

```bash
gmx trjconv -s md_10ns.tpr -f md_10ns.xtc -o md_50ns_center.xtc -center -pbc mol -ur compact
```

During execution, you will be prompted twice:

- Select group for centering – choose the protein (typically group 1).
- Select group for output – choose the system (group 0).
- Now use `md_50ns_center.xtc` for all downstream analysis.
 
### 7.2 Common structural analyses
All analyses below use the corrected trajectory (`md_50ns_center.xtc`) and the run input file (`md_10ns.tpr`).

### 7.2.1 RMSD – Root Mean Square Deviation

Measures how much the protein structure deviates from a reference (usually the starting structure) over time.

```bash
gmx rms -s md_10ns.tpr -f md_50ns_center.xtc -o RMSD.xvg -tu ns
```
At the prompts:
- For least squares fit – choose
- For RMSD calculation – choose

**Output:** `RMSD.xvg` (RMSD in `nm` vs. time in `ns`)

### 7.2.2 Radius of gyration (rGyr)
Indicates the compactness of the protein.

```bash
gmx gyrate -s md_10ns.tpr -f md_50ns_center.xtc -o gyrate.xvg
```

At prompt:
- Select group – choose

**Output:** `gyrate.xvg` (radius of gyration in `nm`)

### 7.2.3 RMSF – Root Mean Square Fluctuation
Per‑residue flexibility (requires an index file with residue numbers).

```bash
gmx rmsf -s md_10ns.tpr -f md_50ns_center.xtc -o rmsf.xvg -res
```

At prompt:
- Select group for fitting – usually 4 (Backbone)
- Select group for RMSF – 4 (Backbone) or 1 (Protein) – use `-res` to get per‑residue output.

**Output:** `rmsf.xvg` (RMSF in `nm` per residue)

### 7.2.4 SASA – Solvent Accessible Surface Area
Measures the surface area exposed to solvent.

```bash
gmx sasa -s md_10ns.tpr -f md_50ns_center.xtc -o sasa.xvg -or resarea.xvg
```
At prompt:
- Select group – 

**Outputs:**
- `sasa.xvg` – total SASA over time (`nm²`)
- `resarea.xvg` – residue‑wise SASA

### 7.3 Automated plotting and advanced analysis with Dynamics‑Visualizer

Instead of manually plotting each `.xvg` file, you can use the **[Dynamics‑Visualizer](https://github.com/usman4373/Dynamics-Visualizer)** – a Streamlit‑based app that automates post‑MD analysis and visualization for GROMACS trajectories.

#### What it does

- **Standard analyses** – Reads `.xvg` files (RMSD, RMSF, SASA, radius of gyration) and produces publication‑ready time‑series plots.
- **Principal Component Analysis (PCA)** – Computes and plots PCA for protein‑only or protein‑ligand systems with time‑based colour bars.
- **Dynamic Cross‑Correlation Matrix (DCCM)** – Generates heatmaps showing correlated motions between residues (and ligand, if present).
- **Trajectory visualizer** – Extracts frames from your `.xtc` trajectory, renders PNGs using PyMOL, and assembles MP4 videos of the simulation.

#### How to use it

1. **Install** the tool (see [its README](https://github.com/usman4373/Dynamics-Visualizer) for conda/pip instructions, including PyMOL).
2. **Prepare your data** – Use the **centered trajectory** (`md_10ns_noPBC.xtc`) and the run input file (`md_10ns.tpr`) from Step 7.1.
3. **Run the app**:
   ```bash
   streamlit run app.py
   ```
