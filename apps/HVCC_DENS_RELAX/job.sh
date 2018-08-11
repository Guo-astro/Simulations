#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l ncpus=64
#PBS -q madduxA
#PBS -N job-isorelax
./hvcc_dens_relax.out > out.log
