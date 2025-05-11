#!/bin/bash
#SBATCH --job-name=name
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40GB
#SBATCH --output=log_test
#SBATCH --error=err_test
#SBATCH --export=none

# Autor Maciej Szyszko
# prosty job do slurm przydatny do puszczania obliczeń dla każdego z podfolderów (np. różne koncentracje domieszek)
# można w miare kontrolować w którym momencie obliczeń jesteśmy na podstawie głównego pliku wyjściowego (--output)
# raczej można to zrobić optymalniej przez job array

ulimit -l unlimited

module load vasp/6.3.2-intel-2021b
#module load vasp/5.4.4-intel-2021b-vtst
cd $SLURM_SUBMIT_DIR

for subdirectory in */; do

    dir_name="${subdirectory%/}"
    echo "Entering directory: $dir_name"
    cd "$dir_name" || { echo "Faile to enter $dir_name"; continue; }


    mpirun vasp_std > "${dir_name}_vasp.out" 2> "${dir_name}_vasp.err"

    echo "Finished VASP in $dir_name"

    cd ..
done
