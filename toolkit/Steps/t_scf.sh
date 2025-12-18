#!/bin/bash -l
#SBATCH --job-name=t_scf
#SBATCH --time=0:05:00
#SBATCH --account=plg2d4cat-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=out
#SBATCH --error=err

ulimit -s unlimited
module load  intel-compilers/2023.2.1  impi/2021.10.0   VASP/6.5.1-Dsingle_prec_bse
module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

START_DIR="$(pwd)"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --input "$INPUT" --step t_scf --part dry)
cd "$output"
mpiexec vasp_std > log

cd "$START_DIR"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --input "$INPUT" --step t_scf --part scf)
cd "$output"
mpiexec vasp_std > log
