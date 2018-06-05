#!/bin/bash

#SBATCH -p L 
#SBATCH -n 256
#SBATCH -c 16
#SBATCH -J MBH5_engmod
#SBATCH -o stdout.log
#SBATCH -e stderr.log

#module load gnu/openmpi165
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
srun ./sph.out
#mpirun -mca btl ^openib -np `expr ${SLURM_NTASKS_PER_NODE} \* ${SLURM_JOB_NUM_NODES}` ./sph.out
RETCODE=$?
exit ${RETCODE}


