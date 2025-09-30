# Example01 - Bulk MoTe2

This example covers the full calculations of the bulk MoTe2 system, including geometry optimisation, convergence tests, band structures and density of states (DOS).
This guide covers how to set up the VASP settings, prepare files and run the calculations on the cluster (with slurm system).

---

## Input files

```
.
├── POSCAR   # VASP input file containing the crystal structure
├── KPOINTS  # VASP input file containing the k-mesh for SCF calculations
├── input.py # toolkit input file containing the setting specified by the user
```

## Running example01

### 1. Configure defaults
Set paths and global defaults in **`defaults.py`**, such as the path to your folder containing all VASP POTCAR files, the settings of your cluster:

```python
PSEUDOPOTENTIALS_PATH = "/path/to/vasp/pseudopotentials"
CALCULATION_PATH = "your/default/Calculation/path"
CPU_PER_NODE = 192
CPU_PER_SOCKET = 24
```

### 2. Define the input.py
Specify system-specific overrides in **`input.py`**. In this example the full workflow will be covered, with the steps in the order as defined below. Each steps in the `STEPS` variable is defined in the way that the value for each key corresponds to the path from which the needed input files will be copied. Additionally some system specific VASP variables need to be set up in INCAR_SETTINGS, as well as the k-meshes and ENCUT values for convergence tests. The full setup file should look like this:

```python
from defaults import *

CALCULATION_PATH = "example01_calc"

STEPS = {
    't_scf': 'your POSCAR & KPOINTS path',    
    't_geo': 't_scf',
    't_conv_test': 't_geo',
    't_bs': 'converged',
    't_dos': 'converged',
    't_scf_so': 'converged',
    't_geo_so': 't_scf_so',
    't_bs_so': 't_geo_so',
    't_dos_so': 't_geo_so'
}

INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
INCAR_SETTINGS['IVDW'] = 11
INCAR_SETTINGS['LMAXMIX'] = 4
INCAR_SETTINGS['LASPH'] = '.TRUE.'

CONV_KMESH = [(6, 6, 2), (7, 7, 2), (8, 8, 3)] 
CONV_ENCUT = [250, 300, 350, 400]
```
After setting up the input files, now you can run the steps: `t_scf`, `t_geo`, `t_conv_test`; as the output of convergence tests will be needed to determine the k-mesh and ENCUT parameters. Running the first 3 steps is described in Section 3.

### 3. Run the first steps
TO DO

---

### 4. Verify the convergence
After the `t_conv_test` step is finished, check the file `t_conv_test/convergence_report.csv` file, and choose the parameters for the rest of the calculations based on the lattice constants and if the particular calculation converged (the external pressure lower than 1 kilobar for all directions and the atomic forces below 1e-3). Here the k-mesh of (8, 8, 3) and ENCUT=350 should be enough. Now that you need to specify the path to SCF calculations with this parameters for the next calculations (band structure, DOS and spin-orbit calculations), by changing the STEPS with 'converged' in the **`input.py`** to the correct path. Additionally you need to specify the `K_PATH` variable for the band structure calculations, here it is defined as Γ-M-K-Γ-A-L-H-A, you can also specify the number of points between each high symmetry point using the `POINTS_PER_SEGMENT` variable. The final **`input.py`** file should look like something like this.

```python
from defaults import *

CALCULATION_PATH = "example01_calc"
STEPS = {
    't_scf': 'path containing POSCAR & KPOINTS',    
    't_geo': 't_scf',
    't_conv_test': 't_geo',
    't_bs': 'example01/t_conv_test/k8x8x3/e350/03_t_scf',
    't_dos': 'example01/t_conv_test/k8x8x3/e350/03_t_scf',
    't_scf_so': 'example01/t_conv_test/k8x8x3/e350/03_t_scf',
    't_geo_so': 't_scf_so',
    't_bs_so': 't_geo_so',
    't_dos_so': 't_geo_so'
}

INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
INCAR_SETTINGS['IVDW'] = 11
INCAR_SETTINGS['LMAXMIX'] = 4
INCAR_SETTINGS['LASPH'] = '.TRUE.'

CONV_KMESH = [(6, 6, 2), (7, 7, 2), (8, 8, 3)] 
CONV_ENCUT = [250, 300, 350, 400]

K_PATH = [
    (0.00000,  0.00000,  0.00000), # Γ
    (0.50000,  0.00000,  0.00000), # M 
    (0.33333,  0.33333,  0.00000), # K
    (0.00000,  0.00000,  0.00000), # Γ
    (0.00000,  0.00000,  0.50000), # A
    (0.50000,  0.00000,  0.50000), # L
    (0.33333,  0.33333,  0.50000), # H
    (0.00000,  0.00000,  0.50000)  # A
]
```
Now that the file for further calculations is prepared you can run 6 remaining steps, which will be described in Section 5.

### 5. Run the rest of the calculations
TO DO

### 6. Analyze the results
Refer to **Postprocessing.ipynb**

