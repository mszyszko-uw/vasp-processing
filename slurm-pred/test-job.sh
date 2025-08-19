#!/bin/sh

#SBATCH --cpus-per-task       1
#SBATCH --error               outpur-%j.err
##SBATCH --hold
#SBATCH --job-name            my_job
#SBATCH --mem                 1G
#SBATCH --nodes               4
#SBATCH --ntasks              16
#SBATCH --output              output-%j.out
#SBATCH --partition           long-gpu
#SBATCH --gres		      gpu:4
#SBATCH --time                48:30:00

module load VASP
vasp_parallel
