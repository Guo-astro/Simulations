#!/bin/bash

#SBATCH -p S
#SBATCH -n 1
#SBATCH -c 32 
#SBATCH -J hydrostatic_gsph
#SBATCH -o stdout.log
#SBATCH -e stderr.log

#module load gnu/openmpi165
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
srun ./hydrostatic_gsph.out
#mpirun -mca btl ^openib -np `expr ${SLURM_NTASKS_PER_NODE} \* ${SLURM_JOB_NUM_NODES}` ./sph.out
RETCODE=$?
exit ${RETCODE}
