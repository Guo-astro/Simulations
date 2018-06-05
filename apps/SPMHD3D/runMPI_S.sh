#!/bin/bash

#SBATCH -p S
#SBATCH -n 1
#SBATCH -c 32 
#SBATCH -J Density_relax_single_timestep_1e6
#SBATCH -o stdout.log
#SBATCH -e stderr.log

#module load gnu/openmpi165
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
srun ./gspmhd.out
#mpirun -mca btl ^openib -np `expr ${SLURM_NTASKS_PER_NODE} \* ${SLURM_JOB_NUM_NODES}` ./sph.out
RETCODE=$?
exit ${RETCODE}
