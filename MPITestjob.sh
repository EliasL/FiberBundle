#!/bin/bash

#SBATCH -J CFL_S00
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --exclusive=user
#SBATCH --time=1-00:00:00

ml eb
ml Julia/1.7.2-linux-x86_64
mpiexecjl -n 4 julia testing/MPI_testing/hello.jl 

#git add --all
#git commit -am "job done"
#git push

wait