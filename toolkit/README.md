# DFT-toolkit

**DFT-toolkit** is a Python package for automating Density Functional Theory (DFT) calculations across multiple DFT software packages.
It allows users to:

* Create, submit, and monitor jobs in the **SLURM** queue.
* Build multi-step calculation pipelines with defined dependencies.
* Work seamlessly with any DFT software on any HPC cluster using SLURM.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/mszyszko-uw/vasp-processing.git
```

Install requirements:

* Python ≥ 3.6
* pip ≥ 19.0
* simple_slurm==0.3.6
* PyYAML==6.0.2

```bash
pip install -r requirements.txt
```

> It is strongly recommended to use a separate Python environment (either `venv` or `conda`).

---

## Quickstart

Run with Python:

```bash
python dft_toolkit.py [options] action
```

List available options:

```bash
python dft_toolkit.py --help
```

---

## Available Actions

* `freenodes` — list all available nodes with free CPUs and memory
* `waiting` — report waiting time in the SLURM queue
* `print` — print SLURM scripts
* `create` — create SLURM scripts
* `submit` — create and submit SLURM scripts
* `array` — run job arrays from all subdirectories
* `checkqueue` — list all user jobs
* `canceljob` — cancel a job
* `jobinfo` — print job details from SLURM
* `print_avalable_steps` — list defined jobs

---

## Available Options

* `-h, --help` — show help message and exit
* `--path PATH` — set a path to the working directory or resource
* `--config CONFIGURE_FILE` — provide a configuration file
* `--id JOB_ID` — specify a SLURM Job ID (for details or cancellation)
* `--steps STEPS` — define list of jobs to run
* `--array` — run job array based on step and path
* `--dependency_step DEPENDENCY_STEP` — specify a blocking step

---

## Configuration Files

DFT-toolkit uses two configuration files:

1. **Machine file** — defines cluster details (YAML or JSON format)
2. **Step file (`steps.yaml`)** — defines pipeline steps

### Machine File

Defines cluster details such as SLURM partitions, module names, and Python environment.

Sections:

* **slurm** — cluster details (nodes, partitions, time, etc.)
* **script** — modules to load in job submission
* **env** — Python environment type (`venv` or `conda`) and path

Example `config.yaml`:

```yaml
slurm:
  nodes: 1
  partition: short
  time: "1:00:00"

script:
  module: vasp/22

env:
  type: venv   # venv or conda
  path: /temp/dft-toolkit/environments/venv
```

### Steps Definition

The `steps.yaml` defines pipeline steps.
Each step includes:

* **slurm** — SLURM job parameters
* **cmd** — SLURM script body (calculations, pre/postprocessing commands)

Example `steps.yaml`:

```yaml
scf:    # step 1
  slurm:
    time: "0:05:00"
    nodes: 1
    ntasks: 28
    mem_per_cpu: "2GB"
  cmd: |
    ulimit -s unlimited
    START_DIR="$(pwd)"
    cd "$output"
    mpiexec vasp_std > log

postprocessing:  # step 2
  slurm:
    time: "0:20:00"
    nodes: 1
    ntasks: 1
    mem_per_cpu: "1GB"
  cmd: |
    cd OUTPUTS
    python ../postprocessing_plot.py
```

---

## Example Workflow

1. **Create virtual environment**

   ```bash
   python -m venv ./toolbox_env
   ```

2. **Clone repository**

   ```bash
   git clone https://github.com/mszyszko-uw/vasp-processing.git
   ```

3. **Install dependencies**

   ```bash
   source toolbox_env/bin/activate
   cd vasp-processing
   pip install -r requirements.txt
   ```

4. **Create machine file**
   Example `config.yaml`:

   ```yaml
   slurm:
     nodes: 1
     partition: __PARTITION__
     time: "1:00:00"
   script:
     module: vasp/22
   env:
     type: venv
     path: __PATH__/toolbox_env
   ```

5. **Create steps file**
   Example `steps.yaml`:

   ```yaml
   scf:
     slurm:
       time: "0:05:00"
       nodes: 1
       ntasks: 28
       mem_per_cpu: "1GB"
     cmd: |
       ulimit -s unlimited
       module load intel-compilers/2023.2.1 impi/2021.10.0 VASP/6.5.1-Dsingle_prec_bse
       module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

       START_DIR="$(pwd)"
       output=$(python3 toolkit.py --step step01 --part dry)
       cd "$output"
       mpiexec vasp_std > log

       cd "$START_DIR"
       output=$(python3 toolkit.py --step step01 --part scf)
       cd "$output"
       mpiexec vasp_std > log
   ```

6. **Examine SLURM script**

   ```bash
   python dft_toolkit.py --steps scf print
   ```

7. **Run calculations**

   ```bash
   python dft_toolkit.py --steps scf submit
   ```

---
