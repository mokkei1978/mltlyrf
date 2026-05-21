#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=4

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4

./flat2case2.exe
