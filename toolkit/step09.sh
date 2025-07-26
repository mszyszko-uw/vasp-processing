#!/bin/bash -l
#SBATCH --job-name=s01e01scf
#SBATCH --time=0:02:00
#SBATCH --account=plg2dmagsem-cpu 
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=out_step09
#SBATCH --error=err_step09

ulimit -s unlimited
module load  intel-compilers/2023.2.1  impi/2021.10.0   VASP/6.5.1-Dsingle_prec_bse
module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

START_DIR="$(pwd)"
output=$(python3 toolkit.py --step step09 --part dry)
cd "$output"
mpiexec vasp_ncl > log

cd "$START_DIR"
output=$(python3 toolkit.py --step step09 --part dos)
cd "$output"
mpiexec vasp_ncl > log

