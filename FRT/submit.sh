#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1

set -e

cd ${PBS_O_WORKDIR}

time ./main.exe

exit 0
