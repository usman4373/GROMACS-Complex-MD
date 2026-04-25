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
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top
```

During the run, you will be prompted to select:
- Force field (e.g., CHARMM36, Amber99SB, OPLS)
- Water model (e.g., TIP3P)

> If writing the topology of a protein fails due to a mismatch of hydrogen atoms, then use `-ignh` in the above command. It will let GROMACS rebuild hydrogens correctly

### 1.2 Process the ligand
- Add hydrogens to the ligand using [Avogadro](https://avogadro.cc/) (Build > Hydrogens > Add hydrogens) and save as `ligand.mol2` (File > Export > Molecule).
- After adding hydrogens, open `ligand.mol2` file in a text editor and make 2 changes.
- File will look like this:
```bash
@<TRIPOS>MOLECULE
*****
 45 47 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C          17.1008   29.8156   -2.3788 C.3     1  UNL1        0.1918
      2 O          16.4529   30.9539   -1.8211 O.3     1  UNL1       -0.3404
      3 C          16.3137   31.9810   -2.7818 C.3     1  UNL1        0.1130
      4 C          15.1994   31.6049   -3.7917 C.3     1  UNL1        0.1119
      5 O          15.4453   32.2673   -5.0067 O.3     1  UNL1       -0.3865
      6 C          15.1480   30.0692   -4.0247 C.3     1  UNL1        0.1141
      7 O          14.1454   29.4672   -3.2398 O.3     1  UNL1       -0.3864
      8 C          16.5147   29.4035   -3.7657 C.3     1  UNL1        0.1342
      9 O          17.4002   29.7755   -4.7923 O.3     1  UNL1       -0.3841
     10 H          14.5897   32.2729   -5.5105 H       1  UNL1        0.2100
     11 H          14.3254   29.6794   -2.2879 H       1  UNL1        0.2100
     12 H          18.1482   29.1226   -4.7770 H       1  UNL1        0.2101
     13 C          15.9926   33.2931   -2.0512 C.3     1  UNL1        0.0730
     14 O          16.3484   34.4034   -2.8308 O.3     1  UNL1       -0.3924
     15 H          17.3381   34.4714   -2.7931 H       1  UNL1        0.2095
     16 O          16.9360   28.7251   -1.4919 O.3     1  UNL1       -0.3090
     17 C          17.9069   28.7532   -0.4612 C.3     1  UNL1        0.2137
     18 O          17.5769   29.7876    0.4670 O.2     1  UNL1       -0.4692
     19 C          18.5060   30.2895    1.3833 C.2     1  UNL1        0.1021
     20 C          19.8166   30.0694    1.2166 C.2     1  UNL1        0.0969
     21 C          20.8230   30.6384    2.1313 C.2     1  UNL1        0.3427
     22 O          20.4890   31.3458    3.1190 O.2     1  UNL1       -0.2448
     23 O          22.1740   30.4295    1.8566 O.3     1  UNL1       -0.4569
     24 C          22.6128   29.4199    0.9462 C.3     1  UNL1        0.1116
     25 C          21.6231   29.1794   -0.1658 C.2     1  UNL1       -0.0422
     26 C          20.2979   29.3512    0.0105 C.2     1  UNL1       -0.0185
     27 C          19.2938   29.0105   -1.0835 C.3     1  UNL1        0.0652
     28 C          19.7007   27.7890   -1.8825 C.2     1  UNL1       -0.0779
     29 C          19.7325   26.5700   -1.3361 C.2     1  UNL1       -0.1019
     30 H          18.1302   30.0710   -2.5201 H       1  UNL1        0.0939
     31 H          17.2217   32.1069   -3.3336 H       1  UNL1        0.0647
     32 H          14.2534   31.9063   -3.3928 H       1  UNL1        0.0647
     33 H          14.9037   29.9172   -5.0553 H       1  UNL1        0.0648
     34 H          16.3830   28.3417   -3.7549 H       1  UNL1        0.0671
     35 H          14.9431   33.3315   -1.8462 H       1  UNL1        0.0584
     36 H          16.5541   33.3222   -1.1408 H       1  UNL1        0.0584
     37 H          17.9252   27.8169    0.0564 H       1  UNL1        0.1066
     38 H          18.1667   30.8550    2.2260 H       1  UNL1        0.1033
     39 H          23.5442   29.7253    0.5172 H       1  UNL1        0.0743
     40 H          22.7147   28.5073    1.4954 H       1  UNL1        0.0743
     41 H          21.9833   28.8632   -1.1225 H       1  UNL1        0.0609
     42 H          19.2625   29.8484   -1.7482 H       1  UNL1        0.0453
     43 H          19.9697   27.9057   -2.9115 H       1  UNL1        0.0574
     44 H          19.4663   26.4373   -0.3083 H       1  UNL1        0.0532
     45 H          20.0242   25.7266   -1.9264 H       1  UNL1        0.0532
@<TRIPOS>BOND
     1     1    16    1
     2     1     2    1
     3     2     3    1
     4     3    13    1
     5     3     4    1
     6     4     5    1
     7     4     6    1
     8     5    10    1
     9     6     7    1
    10     6     8    1
    11     7    11    1
    12     8     9    1
    13     1     8    1
    14     9    12    1
    15    13    14    1
    16    14    15    1
    17    16    17    1
    18    17    18    1
    19    18    19    1
    20    19    20    2
    21    20    21    1
    22    21    22    2
    23    21    23    1
    24    23    24    1
    25    24    25    1
    26    25    26    2
    27    26    27    1
    28    20    26    1
    29    27    28    1
    30    17    27    1
    31    28    29    2
    32     1    30    1
    33     3    31    1
    34     4    32    1
    35     6    33    1
    36     8    34    1
    37    13    35    1
    38    13    36    1
    39    17    37    1
    40    19    38    1
    41    24    39    1
    42    24    40    1
    43    25    41    1
    44    27    42    1
    45    28    43    1
    46    29    44    1
    47    29    45    1
```
### First change
- In the "@<TRIPOS>MOLECULE" heading, replace "*****" with "LIG" or "UNL" etc

### Second change
- Fix the residue names and numbers such that they are all the same
> If the `ligand.mol2` file contains more than 1 molecule, then make changes accordingly

- After both changes, the file should look like this:

```bash
@<TRIPOS>MOLECULE
LIG
 45 47 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C          17.1008   29.8156   -2.3788 C.3     1  LIG         0.1918
      2 O          16.4529   30.9539   -1.8211 O.3     1  LIG        -0.3404
      3 C          16.3137   31.9810   -2.7818 C.3     1  LIG         0.1130
      4 C          15.1994   31.6049   -3.7917 C.3     1  LIG         0.1119
      5 O          15.4453   32.2673   -5.0067 O.3     1  LIG        -0.3865
      6 C          15.1480   30.0692   -4.0247 C.3     1  LIG         0.1141
      7 O          14.1454   29.4672   -3.2398 O.3     1  LIG        -0.3864
      8 C          16.5147   29.4035   -3.7657 C.3     1  LIG         0.1342
      9 O          17.4002   29.7755   -4.7923 O.3     1  LIG        -0.3841
     10 H          14.5897   32.2729   -5.5105 H       1  LIG         0.2100
     11 H          14.3254   29.6794   -2.2879 H       1  LIG         0.2100
     12 H          18.1482   29.1226   -4.7770 H       1  LIG         0.2101
     13 C          15.9926   33.2931   -2.0512 C.3     1  LIG         0.0730
     14 O          16.3484   34.4034   -2.8308 O.3     1  LIG        -0.3924
     15 H          17.3381   34.4714   -2.7931 H       1  LIG         0.2095
     16 O          16.9360   28.7251   -1.4919 O.3     1  LIG        -0.3090
     17 C          17.9069   28.7532   -0.4612 C.3     1  LIG         0.2137
     18 O          17.5769   29.7876    0.4670 O.2     1  LIG        -0.4692
     19 C          18.5060   30.2895    1.3833 C.2     1  LIG         0.1021
     20 C          19.8166   30.0694    1.2166 C.2     1  LIG         0.0969
     21 C          20.8230   30.6384    2.1313 C.2     1  LIG         0.3427
     22 O          20.4890   31.3458    3.1190 O.2     1  LIG        -0.2448
     23 O          22.1740   30.4295    1.8566 O.3     1  LIG        -0.4569
     24 C          22.6128   29.4199    0.9462 C.3     1  LIG         0.1116
     25 C          21.6231   29.1794   -0.1658 C.2     1  LIG        -0.0422
     26 C          20.2979   29.3512    0.0105 C.2     1  LIG        -0.0185
     27 C          19.2938   29.0105   -1.0835 C.3     1  LIG         0.0652
     28 C          19.7007   27.7890   -1.8825 C.2     1  LIG        -0.0779
     29 C          19.7325   26.5700   -1.3361 C.2     1  LIG        -0.1019
     30 H          18.1302   30.0710   -2.5201 H       1  LIG         0.0939
     31 H          17.2217   32.1069   -3.3336 H       1  LIG         0.0647
     32 H          14.2534   31.9063   -3.3928 H       1  LIG         0.0647
     33 H          14.9037   29.9172   -5.0553 H       1  LIG         0.0648
     34 H          16.3830   28.3417   -3.7549 H       1  LIG         0.0671
     35 H          14.9431   33.3315   -1.8462 H       1  LIG         0.0584
     36 H          16.5541   33.3222   -1.1408 H       1  LIG         0.0584
     37 H          17.9252   27.8169    0.0564 H       1  LIG         0.1066
     38 H          18.1667   30.8550    2.2260 H       1  LIG         0.1033
     39 H          23.5442   29.7253    0.5172 H       1  LIG         0.0743
     40 H          22.7147   28.5073    1.4954 H       1  LIG         0.0743
     41 H          21.9833   28.8632   -1.1225 H       1  LIG         0.0609
     42 H          19.2625   29.8484   -1.7482 H       1  LIG         0.0453
     43 H          19.9697   27.9057   -2.9115 H       1  LIG         0.0574
     44 H          19.4663   26.4373   -0.3083 H       1  LIG         0.0532
     45 H          20.0242   25.7266   -1.9264 H       1  LIG         0.0532

The rest of the file will remain identical
```
  
- Now, fix bond order in the `ligand.mol2` file using a [Perl script](http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.txt):

```bash
perl sort_mol2_bonds.pl ligand.mol2 ligand_fix.mol2
```

## 1.3 Generate ligand topology

- Generate ligand topology using [CGenFF](https://cgenff.com/) (Create account to use CGenFF):
    - Upload your `ligand_fix.mol2` file → select parameters → Run Cgenff engine → Convert results to GROMACS format and download
    - Output will contain: `ligand_fix_gmx.pdb`, `ligand_fix_gmx.top`, and charmm36.ff directory
> CGenFF gives a score to each parameter to show how reliable it is.
> - Scores below 10 are usually safe to use directly.
> - Scores between 10 and 50 mean you should check or validate them.
> - Scores above 50 are unreliable and often need to be re-parameterized manually.
>
> These scores are very important because they tell you how much you can trust the generated parameters.

- The files from the CGenFF server need a few changes before use
- The topology file `ligand_fix_gmx.top` is made for full system use in GROMACS
- It needs to be modified further
- First, make a copy of this file and rename it to `ligand_gmx.itp`:

```bash
cp ligand_fix_gmx.top ligand_gmx.itp
```

- Now, make further changes to the `ligand_gmx.itp` file (open in text editor):
1. Remove `#include "./charmm36.ff/forcefield.itp"` line
2. Replace `#include "./charmm36.ff/lig_ffbonded.itp"` with the actual contents of the `charmm36.ff/lig_ffbonded.itp` file
3. Delete everything from `; Include water topology` to the end of the file
4. Rename the `[moleculetype]` from `Other` to `LIG`
5. Replace `#include "posre.itp"` to `#include "posre_ligand.itp"`

## 1.4 Build the protein-ligand complex
- The `protein_processed.gro` file, which was generated during the protein processing step, will be used to create the complex
- It contains the processed, force field-compliant structure of the protein
- For ligand, `ligand_fix_gmx.pdb` file will be used, which was obtained from the CGenFF server that has all of the necessary H atoms and matches the atom names in the LIG topology
- First convert `ligand_fix_gmx.pdb` file into `.gro` format:

```bash
gmx editconf -f ligand_fix_gmx.pdb -o ligand.gro
```
- Now, make a copy of the `protein_processed.gro` file:

```bash
cp protein_processed.gro complex.gro
```

- Open the `ligand.gro` file in a text editor, copy its coordinate part, and paste it into the `complex.gro` file
- Place it right after the last protein atom line and before the box vector lines:

```cpp
  374PHE    OT2 5811   1.473   1.005  -1.391     # Last protein atom line
    1lig     C1    1   1.710   2.982  -0.238
    1lig     O2    2   1.645   3.095  -0.182
    1lig     C3    3   1.631   3.198  -0.278
    1lig     C4    4   1.520   3.161  -0.379
    1lig     O5    5   1.544   3.227  -0.501
    1lig     C6    6   1.515   3.007  -0.403
    1lig     O7    7   1.414   2.947  -0.324
    1lig     C8    8   1.651   2.940  -0.377
    1lig     O9    9   1.740   2.978  -0.479
    1lig    H10   10   1.459   3.227  -0.551
    1lig    H11   11   1.433   2.968  -0.229
    1lig    H12   12   1.815   2.912  -0.478
    1lig    C13   13   1.599   3.329  -0.205
    1lig    O14   14   1.635   3.440  -0.283
    1lig    H15   15   1.734   3.447  -0.279
    1lig    O16   16   1.694   2.872  -0.149
    1lig    C17   17   1.791   2.875  -0.046
    1lig    O18   18   1.758   2.979   0.047
    1lig    C19   19   1.851   3.029   0.138
    1lig    C20   20   1.982   3.007   0.122
    1lig    C21   21   2.082   3.064   0.213
    1lig    O22   22   2.049   3.135   0.312
    1lig    O23   23   2.217   3.043   0.186
    1lig    C24   24   2.261   2.942   0.095
    1lig    C25   25   2.162   2.918  -0.017
    1lig    C26   26   2.030   2.935   0.001
    1lig    C27   27   1.929   2.901  -0.108
    1lig    C28   28   1.970   2.779  -0.188
    1lig    C29   29   1.973   2.657  -0.134
    1lig    H30   30   1.813   3.007  -0.252
    1lig    H31   31   1.722   3.211  -0.333
    1lig    H32   32   1.425   3.191  -0.339
    1lig    H33   33   1.490   2.992  -0.506
    1lig    H34   34   1.638   2.834  -0.375
    1lig    H35   35   1.494   3.333  -0.185
    1lig    H36   36   1.655   3.332  -0.114
    1lig    H37   37   1.793   2.782   0.006
    1lig    H38   38   1.817   3.086   0.223
    1lig    H39   39   2.354   2.973   0.052
    1lig    H40   40   2.272   2.851   0.149
    1lig    H41   41   2.198   2.886  -0.112
    1lig    H42   42   1.926   2.985  -0.175
    1lig    H43   43   1.997   2.791  -0.291
    1lig    H44   44   1.947   2.644  -0.031
    1lig    H45   45   2.002   2.573  -0.193
   5.77620   6.44341   6.98742   # box vector line
```
- Because 45 atoms were added to the `complex.gro` file, update the second line of `complex.gro` by increasing the atom count by 45 (In general, add the number of atoms based on how many atoms are in your ligand)

## 1.5 Build the topology of protein-ligand complex
- The ligand introduces new dihedral parameters, so the ligand topology must be included at the top of `topol.top`
- Add the following after `#include "./charmm36-feb2026_cgenff-5.0.ff/forcefield.itp"` line, as shown below:

```cpp
; Include forcefield parameters
#include "./charmm36-feb2026_cgenff-5.0.ff/forcefield.itp"

; Include ligand parameters    | Add these |
#include "ligand_gmx.itp"      | lines     |
```
> <p style="justify;"> The "Include ligand parameters" statement must be positioned carefully. It should be placed above any [ moleculetype ] sections to ensure all parameters are available before any molecules are defined. At the same time, it needs to come after the main force field "Include forcefield parameters", since atom types must already be declared before any bonded parameters that rely on them are introduced </p>

- The final change is required in the [ molecules ] section
- Since a new molecule has been added to `complex.gro`, it must also be included in this directive:

```cpp
[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG                 1
```

The topology and coordinate files now match in terms of the system composition, and the system is ready for solvation

## Step 2: Solvation
### 2.1 Define the simulation box

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

`-cp` = solute configuration (complex in box) `-cs` = solvent configuration (SPC216 water model)

> `Back Off! I just backed up topol.top to ./#topol.top.1#` message while running the above command is normal and indicates that a backup of the original file (e.g., #topol.top.1#) has been created. The original file is then updated with new changes, while the previous version remains safely stored as a backup.

## Step 3: Add ions to neutralize the system
- First, assemble a `.tpr` file using [ions.mdp](http://www.mdtutorials.com/gmx/complex/Files/ions.mdp):

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```

Then replace solvent molecules with ions:

```bash
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

> When prompted, select the `SOL` group (water) – this ensures ions are placed only in the solvent region

## Step 4: Energy minimization
- After assembling the system, generate the binary input file using `grompp` with the provided input parameter file ([em.mdp](http://www.mdtutorials.com/gmx/complex/Files/em.mdp))
```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
```
```bash
gmx mdrun -v -s em.tpr -deffnm em -nb gpu
```

## Step 5: Equilibration (NVT and NPT)

For a protein–ligand complex, position restraints are applied to the ligand (and optionally to the protein backbone) and use separate temperature coupling groups.

### 5.1 Create index file and ligand restraints
```bash
gmx make_ndx -f ligand.gro -o index_lig.ndx
```
- You will be prompted while running the above command. Type the following options:
```bash
> 0 & ! a H*  # Type these and press enter
> q           # Then type this and press enter
```
> This creates a new group containing all non-hydrogen atoms of the entire system

Then generate restraints for the ligand:

```bash
gmx genrestr -f ligand.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
```
- The above command will prompt to select a group:
```cpp
Reading structure file
Select group to position restrain
Group     0 (         System) has    45 elements
Group     1 (          Other) has    45 elements
Group     2 (            lig) has    45 elements
Group     3 (   System_&_!H*) has    25 elements
Select a group: 2
Selected 2: 'lig'
```
> Select ligand group (2)

- Now edit `topol.top` to include `#include "posre_lig.itp"` before [ molecules ]:

```cpp
; Include topology for ions
#include "./charmm36-feb2026_cgenff-5.0.ff/ions.itp"

; Ligand position restraints     | Add   |
#ifdef POSRES_LIG                | these |
#include "posre_lig.itp"         | lines |
#endif                           | here  |

[ system ]
; Name
Protein in water

[ molecules ]
; Compound        #mols
Protein_chain_A     1
LIG                 1
SOL         20461
NA               12
```

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
- Run the NVT equilibration step using the [nvt.mdp](http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp) file
```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
```
```bash
gmx mdrun -v -s nvt.tpr -deffnm nvt -nb gpu -pme gpu -bonded gpu
```
> The `-r` flag uses the minimized coordinates as a reference for restraints.

### 5.4 NPT equilibration (constant pressure)
- Run the NPT equilibration step using the [npt.mdp](http://www.mdtutorials.com/gmx/complex/Files/npt.mdp) file
```bash
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
```
```bash
gmx mdrun -v -s npt.tpr -deffnm npt -nb gpu -pme gpu -bonded gpu
```

## Step 6: Production MD simulation
- Run the NPT equilibration step using the [md.mdp](http://www.mdtutorials.com/gmx/complex/Files/md.mdp) file
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_10ns.tpr
```
```bash
gmx mdrun -v -s md_10ns.tpr -deffnm md_10ns -nb gpu -pme gpu -bonded gpu
```

> Use `-deffnm` to set the base name of all output files (e.g., `md_50ns.gro`, `md_50ns.xtc`)
>
> To run `50 ns` or a different simulation time, modify nsteps in `md.mdp` accordingly
>
> 10,000 ps = 10 ns
>
> To run 50 ns:
>
> nsteps = 50000/0.002 ​= 25000000 nsteps


## Resume crashed simulation
- To resume a crashed simulation run, add `-cpi` argument and its checkpoint file `md_10ns.cpt`
```bash
gmx mdrun -v -s md_10ns.tpr -deffnm md_10ns -cpi md_10ns.cpt -nb gpu -pme gpu -bonded gpu
```

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
