# VASP Workflow Toolkit

This repository provides a **Python-based toolkit** for setting up, managing, and analyzing **VASP density functional theory (DFT) calculations**.  
The central script, `toolkit.py`, automates common VASP workflows such as SCF runs, geometry optimization, band structure, density of states, and convergence tests (with and without spin–orbit coupling).

---

## 📂 Repository Structure

```
.
├── toolkit.py   # Main workflow manager
├── input.py     # User-specific calculation settings
├── defaults.py  # Default parameters (used if not overridden in input.py)
```

---

## ⚙️ Features

- **Step-based workflow control**: each stage (SCF, geometry optimization, DOS, etc.) is modular.  
- **Supports SOC & non-SOC workflows**.  
- **Convergence tests** over ENCUT and k-meshes.  
- **Automatic INCAR handling**:
  - Adds/removes tags depending on the step.  
  - Generates missing input files when possible.  
- **POTCAR management** from recommended pseudopotentials.  
- **Automatic parallelization settings** (`NCORE`, `KPAR`).  
- **Report generation**:
  - PDF summaries for geometry optimization and SCF runs.  
  - CSV convergence reports.  

---

### Dependencies
- Python ≥ 3.9  
- `numpy`  
- `matplotlib`  
- `h5py`  

---

## 📝 Usage

### 1. Configure defaults
Set paths and global defaults in **`defaults.py`**, e.g.:

```python
PSEUDOPOTENTIALS_PATH = "/path/to/vasp/pseudopotentials"
CALCULATION_PATH = "your/Calculation/path"
CPU_PER_NODE = 192
CPU_PER_SOCKET = 24
```

### 2. Define your calculation
Specify system-specific overrides in **`input.py`**.  
At minimum, the `STEPS` dictionary must be defined. For each step either the previous one, or a path to the folder with necessary files needs to be provided:

```python
STEPS = {
    't_scf': 'path/to/POSCAR',    
    't_geo': 't_scf',
    't_bs': 't_geo'
}

INCAR_SETTINGS['SYSTEM'] = 'bulk MoTe2'
INCAR_SETTINGS['ISMEAR'] = 0
INCAR_SETTINGS['SIGMA'] = 0.05
```

### 3. Run a workflow step
Use `toolkit.py` with `--step` and `--part`:

```bash
python toolkit.py --step t_scf --part dry
python toolkit.py --step t_scf --part scf
python toolkit.py --step t_geo --part cg_opt
```

Steps and parts correspond to different stages of a workflow.

Alternatively you can you the functions corresponding to specific steps in explicitly in python:

```python
path_to_calc = t_bs(folder="your/calc/dir", part="scf", previous="your/scf/folder") 
```
This way of using the toolkit can be beneficial for preprocessing files for multiple similar calculations at the same time.
However, this implementation still needs the settings to be defined as previously described.

---

## 📖 Supported Steps

### Steps (non SOC & SOC variants)
- `t_scf`: non SOC SCF calculation, includes parts: 'dry' & 'scf'
- `t_scf_so`:  SOC SCF calculation, includes parts: 'dry' & 'scf'  
- `t_geo`: non SOC geometry optimisation (Conjugate Gradient → quasi-Newton → SCF), includes parts: 'dry', 'cg_opt', 'nw_opt', 'scf' & 'report'
- `t_geo_so`: SOC geometry optimisation (Conjugate Gradient → quasi-Newton → SCF), includes parts: 'dry', 'cg_opt', 'nw_opt', 'scf' & 'report'
- `t_conv_test`: non SOC convergence studies (with geometry optimisation) over ENCUT and k-mesh, includes parts: 'dry', 'cg_opt', 'nw_opt', 'scf' & 'report'  
- `t_bs`: non SOC band structure, includes parts: 'dry', 'bs'
- `t_bs_so`: SOC band structure, includes parts: 'dry', 'bs'  
- `t_dos`: non SOC density of states (DOS), includes parts: 'dry', 'dos'
- `t_dos_so`: SOC density of states (DOS), includes parts: 'dry', 'dos'   

---

## 📊 Reports

- **PDF reports** include stress, forces, and SCF convergence plots.  
- **Convergence reports** (CSV) summarize lattice constants, forces, energy, and convergence flags for different ENCUT and k-mesh values.  

---

## ⚠️ Notes

- Ensure **POSCAR** is always present.  
- If **POTCAR** is missing, the toolkit will attempt to generate it from recommended pseudopotentials.  
- **CHGCAR** is required only for NSCF runs (band structure/DOS).  
- Missing files are handled with fallbacks where possible.
- For band structure calculations a **KPOINTS** file containing a path should be present in the previous directory or K_PATH variable in **input.py** need to be specified
  
