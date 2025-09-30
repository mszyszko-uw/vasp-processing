#!/bin/bash -l
#SBATCH --job-name=t_geo
#SBATCH --time=0:20:00
#SBATCH --account=plg2dmagsem-cpu 
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=out_step01
#SBATCH --error=err_step01

ulimit -s unlimited
module load  intel-compilers/2023.2.1  impi/2021.10.0   VASP/6.5.1-Dsingle_prec_bse
module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

START_DIR="$(pwd)"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --step t_geo --part dry)
cd "$output"
mpiexec vasp_std > log

cd "$START_DIR"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --step t_geo --part cg_opt)
cd "$output"
mpiexec vasp_std > log

cd "$START_DIR"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --step t_geo --part nw_opt)
cd "$output"
mpiexec vasp_std > log

cd "$START_DIR"
output=$(python3 $START_DIR/Preprocessing/toolkit.py --step t_geo --part scf)
cd "$output"
mpiexec vasp_std > log

cd "$START_DIR"
python3 $START_DIR/Preprocessing/toolkit.py --step t_geo --part report

