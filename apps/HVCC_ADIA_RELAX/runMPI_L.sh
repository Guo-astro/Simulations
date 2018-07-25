#!/bin/bash


#SBATCH -p M 
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 64
#SBATCH -J HRS_HVCC_IMBH_CHEM
#SBATCH -o stdout.log
#SBATCH -e stderr.log

#module load gnu/openmpi165
export OMP_NUM_THREADS=32
srun ./hvcc_dens_relax.out




