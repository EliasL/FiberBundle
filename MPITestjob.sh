#!/bin/bash

#SBATCH -J Testing
#SBATCH -p porelab
#SBATCH -N 4
#SBATCH -n 3
#SBATCH --exclusive=user
#SBATCH --time=0-00:01:00

ml eb
ml Julia/1.7.2-linux-x86_64

# Update path to include mpiexecjl
export PATH="/home/vemundlu/.julia/bin/:$PATH"

mpiexecjl -n 4 julia testing/MPI_testing/hello.jl 

#git add --all
#git commit -am "job done"
#git push
wait