#!/bin/bash

#SBATCH --nodes=16
#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --job-name=aperture_calc
#SBATCH --output=aperture_calc.out

module load mpi4py

time mpirun -np 16 python parallelizer.py bigdata A output_A.txt
time mpirun -np 16 python parallelizer.py bigdata B output_B.txt
time mpirun -np 16 python parallelizer.py bigdata C output_C.txt
time mpirun -np 16 python parallelizer.py bigdata D output_D.txt