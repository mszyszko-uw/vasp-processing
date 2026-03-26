# DFTbox

Toolkit for managing and postprocessing VASP DFT calculations on HPC/SLURM clusters.

## Installation

Using conda:

```
bashconda create --name dftbox_env python==3.12
conda activate dftbox_env
conda install conda-build
CONDA_BLD_PATH=~/conda-bld conda-build recipe/ -c conda-forge -c defaults
conda install --use-local dftbox

