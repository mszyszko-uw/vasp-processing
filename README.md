# vasp-processing — VASP preprocessing & postprocessing toolkit

**Purpose:** this repository bundles tools for preparing, submitting and postprocessing VASP calculations on HPC (SLURM) clusters and on small local machines.
It combines two major components:

* **DFT-toolkit** — pipeline manager for creating SLURM jobs, arrays, and multi-step workflows.
* **vaspout_h5.py** — focused postprocessing utilities that read `vaspout.h5`, extract MAGMOM, band structure (BS), density-of-states (DOS) and produce publication-ready plots (PNG) and optional gnuplot `.txt` files.
* **toolkit.py** — preprocessing module for creating VASP input files, handling missing files and some reports of the calculations (geometry optimisation and convergence tests).

---

## Table of contents

* [Assumptions (global)](#assumptions-global)
* [Installation & environment](#installation--environment)
* [Repository contents & important files](#repository-contents--important-files)
* [Per-task resource guidelines & IO rules](#per-task-resource-guidelines--io-rules)
* [Using DFT-toolkit (pipeline manager) — quickstart](#using-dft-toolkit-pipeline-manager---quickstart)
* [SLURM job templates & examples (single job, array job)](#slurm-job-templates--examples-single-job-array-job)
* [vaspout_h5.py — postprocessing & plotting (usage and options)](#vaspout_h5py---postprocessing--plotting-usage-and-options)
* [Batch processing helper: `mass_process.sh` (verbatim)](#batch-processing-helper-mass_processsh-verbatim)
* [Outputs & reports](#outputs--reports)
* [Troubleshooting & tips](#troubleshooting--tips)

---

## Assumptions (global)

* VASP **6.5.1 or newer** is used.
* Calculations, pre- and postprocessing are normally performed on **HPC clusters with SLURM**. Local, few-core machines can be used for lightweight tasks / debugging.
* Postprocessing is done primarily from **`vaspout.h5`** whenever possible.
* Support for `vaspwave.h5` (replacement for WAVECAR/CHGCAR) is planned but **not yet available** — treat WAVECAR/CHGCAR as necessary fallbacks until `vaspwave.h5` becomes supported.
* Use **SLURM array jobs** for ensembles of independent calculations.
* The user provides starting **POSCAR** files.
* Use recommended **PBE POTCARs** kept in a known directory on the cluster (configured in `defaults.py`).
* After each task finishes (successfully or crashed) a short **report** is produced (PDF/CSV/LOG depending on step).

---

## Installation & environment

1. **Create a Python environment** (recommended: `venv` or `conda`).

   ```bash
   python -m venv ./venv
   source ./venv/bin/activate
   ```

2. **Install requirements** from the repository:

   ```bash
   pip install -r requirements.txt
   ```

   Requirements: Python ≥ 3.9, `simple_slurm`, `PyYAML`, `h5py`, `numpy`, `matplotlib`.

3. **Cluster modules**: your SLURM script should load the appropriate compiler/MPI/VASP modules (example later).

4. **Pseudopotentials**: make sure the PBE POTCAR directory is available on the cluster and update `defaults.py` to point to it.
   
5. **Configure defaults.py**: go to `defaults.py` and configure your cluster settings, paths for log file, pseudopotentials, as well as your default VASP settings for various calculations. 

---

## Repository contents & important files

* `toolkit.py` — preprocessing module (preparing VASP input files fo calculations).
* `dft_toolkit.py` — pipeline manager for building and submitting SLURM scripts (see `DFT-toolkit.pdf`).
* `defaults.py` — global defaults (POTCAR path, CPU per node/socket, common INCAR defaults).
* `input.py` — user-specific overrides; define `STEPS` and per-step settings here.
* `vaspout_h5.py` — postprocessing module (band/DOS plotting, MAGMOM extraction).
* `mass_process.sh` — parallel wrapper to run `vaspout_h5.py` in many directories.
* `README.md` and `DFT-toolkit.pdf` — project documentation (this README builds on both).
* `Plotting.py`, `InformationCollector.py` — private modules used by `vaspout_h5.py` (must be available on `PYTHONPATH`).

---

## Per-task resource guidelines & IO rules

These are *recommended* resources and IO behaviors. Adjust for your cluster.

### Single-task summaries (recommended)

* **t_dry / t_dry_so**

  * cores: `1`
  * mem/core: `< 1 GB`
  * time: seconds
  * vasp: any
  * exec: `std` or `ncl`
  * outputs: `vaspout.h5` (also `OUTCAR`, `IBZKPT` optional)
  * purpose: quick check; run *without SLURM* for local checks.

* **t_scf / t_scf_so**

  * cores: `10–100`
  * mem/core: `< 2 GB`
  * time: `10 s – 10 h` (longer for `so`)
  * exec: `std` (or `ncl` for SOC)
  * outputs: `vaspout.h5`, `vaspwave.h5` (planned), `CHGCAR`, `WAVECAR`, `EIGENVAL`, logs
  * purpose: obtain TOTEN, charge density, wavefunction.

* **t_geo / t_geo_so**

  * cores: `10–100`
  * mem/core: `< 2 GB`
  * time: `1 min – 100 h` (optimizations are expensive)
  * exec: `std` (or `ncl` for SOC)
  * outputs: `CONTCAR`, `OUTCAR`, `vaspout.h5`, `vaspwave.h5`

* **t_bs / t_bs_so**, **t_dos / t_dos_so**

  * cores: `10–100`
  * mem/core `< 2 GB`
  * time `1 min – 10 h`
  * inputs require `CHGCAR` from previous SCF
  * outputs: `EIGENVAL`, `vaspout.h5`, `PROCAR` (optional), `DOSCAR` for DOS

### File retention rules (policy: what we KEEP / REMOVE)

**We ALWAYS KEEP:**

* `INCAR`, `POSCAR`, `POTCAR`, `KPOINTS`, `EIGENVAL`, `log`, `OSZICAR`, `OUTCAR`, `*.h5` (including `vaspout.h5`), `IBZKPT`, `err`, `job.sh`

**We ALWAYS REMOVE:**

* `REPORT`, `XDATACAR`, `PCDAT`

**Depending on task (parsed and then removed):**

* `WAVECAR`, `CHG*`, `WAVEDER*`, `PROCAR`, `vasprun.xml`, `BSEFATBAND`, `CONTCAR`
  (these are removed when they are not needed to keep repository size reasonable)

> Note: keep `vaspout.h5` unless you have explicit reasons to prune; the toolkit & postprocessing work from `vaspout.h5`.

---

## Using DFT-toolkit (pipeline manager) — quickstart

DFT-toolkit supports creating SLURM scripts from configurable machine and steps YAML files and can submit/monitor/cancel jobs.

1. **Create machine file** (YAML) describing cluster: partitions, modules, environment path. Example `config_scp.yaml`:

   ```yaml
   slurm:
     nodes: 1
     partition: short
     time: "1:00:00"
   script:
     module: vasp/6.5.1
   env:
     type: venv
     path: /path/to/venv
   ```

2. **Create steps file** (`steps.yaml`) with per-step `slurm` and `cmd` sections. Example minimal `scf` step:

   ```yaml
   scf:
     slurm:
       time: "0:05:00"
       nodes: 1
       ntasks: 28
       mem_per_cpu: "2GB"
     cmd: |
       ulimit -s unlimited
       START_DIR="$(pwd)"
       output=$(python3 $START_DIR/Preprocessing/toolkit.py --step scf --part dry)
       cd "$output"
       mpiexec vasp_std > log
   ```

3. **Run toolkit actions**:

   * Create scripts: `python dft_toolkit.py --path /path/to/work --config config_scp.yaml create`
   * Submit: `python dft_toolkit.py --steps scf submit`
   * Check queue: `python dft_toolkit.py checkqueue`
   * Run arrays: `python dft_toolkit.py --array --steps a_test submit`

4. **Local usage** (no SLURM): call the functions from `toolkit.py` or run the CLI to produce input files and run VASP locally (manually or via `mpiexec`), described in `toolkit/Preprocessing/README.md`.

---

## SLURM job templates & examples

### Single-job (t_scf) job script template (example)

```bash
#!/bin/bash -l
#SBATCH --job-name=t_scf
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=out_t_scf
#SBATCH --error=err_t_scf

module load intel/2023 impi/2021
module load VASP/6.5.1
source /path/to/venv/bin/activate

ulimit -s unlimited
START_DIR="$(pwd)"

# Preprocessing (create inputs)
output=$(python3 $START_DIR/Preprocessing/toolkit.pytoolkit.py --step t_scf --part dry --path "$START_DIR")
cd "$output"

# Run VASP
mpiexec vasp_std > log

# Postprocess
python3 postprocessing_plot.py
```

### Array-job example (ENCUT / k-point convergence)

Use `a_test` with a step file defining the array of input folders or with `--array` mode:

```bash
#SBATCH --array=1-20%8   # run up to 8 parallel array tasks
```

The toolkit provides `--array` actions to generate per-task folders and submission commands.

---

## vaspout_h5.py — postprocessing & plotting (usage and options)

`vaspout_h5.py` works from `vaspout.h5` and provides:

* `vaspout_h5(vaspout: str='vaspout.h5')` class with:

  * `extract_magmom()` — writes `MAGMOM` (INCAR-compatible)
  * `read_BS()` / `read_DOS()` — populate internal `InformationCollector`
  * `plot_BS(description='', kpaths='', bands='', **kwargs)` — wrapper to produce BS plots
  * `plot_DOS(description='', ax=None, E0=None, gnuplot=False, **kwargs)` — DOS plots
  * `plot_BS_default()` / `plot_DOS_default()` — quick saves

* `plot_BS_DOS(BS, DOS, ...)` — side-by-side BS + DOS figure

### CLI usage

```bash
python vaspout_h5.py band "1,2 Mo dxy" --kpaths "1 2 3 -4" --bands "3-10" --E0 fermi --save
python vaspout_h5.py dos "total" --E0 0 --gnuplot --save
```

### `description` / `kpaths` / `bands` grammar (summary)

* `description` — e.g. `"1,2 Mo x p dxy"`: atom indices (1-based), species, direction (`x`,`y`,`z`) optional, orbitals (`s,p,px,py,pz,dxy,...`). **Do not mix species** in one description string (e.g., `"Mo Te"` is invalid). `"total"` plots total DOS.
* `--kpaths` — space-separated segment indices; negative reverses path (e.g. `"1 2 3 -4"`).
* `--bands` — `"5"`, `"3-10"`, `"3,5,8"`, `"3-6,9"`.

### Band-plot kwargs (detailed)

`plot_BS` delegates to `BandStructurePlot(ic, **kwargs)`. Key kwargs:

* `min_diff` (float or None): threshold for grouping bands for numbering (default: `(Emax-Emin)/100`)
* `ax` (matplotlib axis): draw onto this axis
* `E0` (`None` | `'fermi'` | float): energy reference (`None` subtracts valence-band maximum; `'fermi'` subtracts `Efermi`)
* `gnuplot` (bool): write gnuplot-ready `.txt` files via `BandToGnuplot(...)`
* `folder` (str): output folder (defaults to `vaspout.h5` folder)
* `save` (bool): save figure (filename derived from `description`)
* `bandnums` (bool): annotate band numbers / ranges at the right edge
* `color` (matplotlib color)
* `mult` (float): multiplier for scatter marker sizes (dot area ∝ `abs(Pmat)*100*mult`)
* `alpha` (float): scatter opacity
* `linestyle` (str or dash tuple): used for `"clean"` (non-projected) plots

**Behavior notes:**

* If `description` contains `'clean'`, plain band lines are plotted (no projection dots).
* If not `'clean'`, band lines are black dotted and projections are shown as sized scatter dots.
* If `gnuplot=True`, `BandToGnuplot` receives processed matrices and writes files to `folder`.
* `save=True` -> saved file base: `folder + description.replace(' ', '_')`.

---

## Batch processing helper: `mass_process.sh` (verbatim)

This script finds subdirectories by name and runs `vaspout_h5.py` inside them, throttling the number of concurrent jobs.

```bash
#!/bin/bash
# Usage: ./mass_process.sh max_jobs parent_folder calc_folder [args...]
# Runs ./vaspout_h5.py inside every matching subdir, keeping at most max_jobs in parallel.

MAX_JOBS=$1
parent_folder=$2
calc_folder=$3
shift 3
script_args=("$@")

if [[ -z "$MAX_JOBS" || -z "$parent_folder" || -z "$calc_folder" ]]; then
    echo "Usage: $0 max_jobs parent_folder calc_folder [args...]"
    exit 1
fi

# Fixed script to execute inside each subdir
START_DIR="$(pwd)/${0%/*}"
SCRIPT="$START_DIR/vaspout_h5.py"

# Function to throttle jobs to MAX_JOBS
function throttle() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
}

# Find all subdirectories named $calc_folder and run in parallel
find "$parent_folder" -type d -name "$calc_folder" | while read -r dir; do
    throttle
    (
        echo ">>> Entering $dir"
        cd "$dir" || exit 1
        "python3" "$SCRIPT" "${script_args[@]}"
    ) &
done

wait
echo "All jobs finished."
```

**Important:** this script **matches directory names**, not directories that contain `vaspout.h5`. If a matched directory lacks `vaspout.h5` the invoked `vaspout_h5.py` will likely fail. To process only directories containing `vaspout.h5`, use:

```bash
find /path/to/parent -type f -name vaspout.h5 -printf '%h\n' | sort -u | while read -r dir; do
  # run per-dir actions
done
```

---

## Outputs & reports

**Reports:** after most job finishes (success or crash), the toolkit produces a short report summarizing:

* convergence (SCF iterations, forces),
* final energies (TOTEN, bandgap if relevant),
* whether the job crashed / reason if parseable,
* links to plots (BS/DOS) and raw `vaspout.h5`.

Reports may be a combination of:

* PDF summary (geometry optimization/SCF),
* CSV convergence tables (for array tests),
* plain-text job logs.

---

## Troubleshooting & tips

* Ensure `Plotting.py` and `InformationCollector.py` are present and importable by `vaspout_h5.py`.
* For BS/DOS postprocessing always try to use `vaspout.h5`. If `vaspwave.h5` exists, note: *full vaspwave support is not yet implemented*.
* For array jobs on clusters, use `a_*` tasks and the DFT-toolkit `--array` feature; monitor SLURM queue with `squeue -u $USER` or DFT-toolkit `checkqueue`.
* If band numbering looks odd, adjust `min_diff` in `plot_BS(...)`.
* If plotting fails due to missing keys in `vaspout.h5`, run `t_dry` locally to verify inputs before submitting heavy jobs.
* Make sure your `config_scp.yaml` points to correct VASP module and Python environment.
