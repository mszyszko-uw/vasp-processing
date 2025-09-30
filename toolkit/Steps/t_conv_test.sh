#!/bin/bash -l
#SBATCH --job-name=t_conv_test
#SBATCH --time=1:30:00
#SBATCH --account=plg2dmagsem-cpu 
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=out_step03
#SBATCH --error=err_step03

ulimit -s unlimited
module load  intel-compilers/2023.2.1  impi/2021.10.0   VASP/6.5.1-Dsingle_prec_bse
module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

START_DIR="$(pwd)"


mapfile -t my_paths < <(python3 $START_DIR/Preprocessing/toolkit.py --step t_conv_test --part dry)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  mpiexec vasp_std > log
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 $START_DIR/Preprocessing/toolkit.py --step t_conv_test--part cg_opt)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  mpiexec vasp_std > log
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 $START_DIR/Preprocessing/toolkit.py --step t_conv_test --part nw_opt)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  mpiexec vasp_std > log
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 $START_DIR/Preprocessing/toolkit.py --step t_conv_test --part scf)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  mpiexec vasp_std > log
  cd "$START_DIR"
done

cd "$START_DIR"
python3 $START_DIR/Preprocessing/toolkit.py --step t_conv_test --part report

