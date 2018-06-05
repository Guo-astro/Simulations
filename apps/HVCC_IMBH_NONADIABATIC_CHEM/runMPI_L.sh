#!/bin/bash


#SBATCH -p M 
#SBATCH -n 8
#SBATCH -N 4
#SBATCH -c 32
#SBATCH -J HYBRID_HRS_HVCC_IMBH_CHEM
#SBATCH -o stdout.log
#SBATCH -e stderr.log

#module load gnu/openmpi165
export OMP_NUM_THREADS=16
srun ./hydrostatic_gsph.out


