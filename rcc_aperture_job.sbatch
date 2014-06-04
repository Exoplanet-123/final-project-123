#!/bin/bash

#SBATCH --nodes=16
#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --job-name=aperture_calc
#SBATCH --output=aperture_calc.out

module load mpi4py

time mpirun -np 16 python parallelizer.py /tmp/hannah\&zakir/dataset01 A output_A.txt
time mpirun -np 16 python parallelizer.py /tmp/hannah\&zakir/dataset01 B output_B.txt
time mpirun -np 16 python parallelizer.py /tmp/hannah\&zakir/dataset01 C output_C.txt
time mpirun -np 16 python parallelizer.py /tmp/hannah\&zakir/dataset01 D output_D.txt