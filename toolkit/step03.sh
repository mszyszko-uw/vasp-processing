#!/bin/bash -l
#SBATCH --job-name=s03
#SBATCH --time=1:30:00
#SBATCH --account=plg2dmagsem-cpu 
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=out_step03_test
#SBATCH --error=err_step03_test

ulimit -s unlimited
module load  intel-compilers/2023.2.1  impi/2021.10.0   VASP/6.5.1-Dsingle_prec_bse
module load Python/3.11.5 matplotlib/3.8.2 h5py/3.11.0

START_DIR="$(pwd)"


mapfile -t my_paths < <(python3 toolkit.py --step step03 --part dry)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  timeout 2m mpiexec vasp_std > log
  if [[ $? -eq 124 ]]; then
    echo "VASP run in $p timed out after 5 minutes." >> "$START_DIR/timed_out.log"
  fi
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 toolkit.py --step step03 --part cg_opt)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  timeout 10m mpiexec vasp_std > log
  if [[ $? -eq 124 ]]; then
    echo "VASP run in $p timed out after 5 minutes." >> "$START_DIR/timed_out.log"
  fi
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 toolkit.py --step step03 --part nw_opt)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  timeout 10m mpiexec vasp_std > log
  if [[ $? -eq 124 ]]; then
    echo "VASP run in $p timed out after 5 minutes." >> "$START_DIR/timed_out.log"
  fi
  cd "$START_DIR"
done

cd "$START_DIR"
mapfile -t my_paths < <(python3 toolkit.py --step step03 --part scf)
for p in "${my_paths[@]}"; do
  cd "$p"
  echo "$p"
  timeout 5m mpiexec vasp_std > log
  if [[ $? -eq 124 ]]; then
    echo "VASP run in $p timed out after 5 minutes." >> "$START_DIR/timed_out.log"
  fi
  cd "$START_DIR"
done

cd "$START_DIR"
python3 toolkit.py --step step03 --part report

