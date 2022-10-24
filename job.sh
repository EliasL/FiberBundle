#!/bin/bash

#SBATCH -J estimateTime
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive=user
#SBATCH --time=0-05:00:00

ml eb
#ml av
ml Julia/1.7.2-linux-x86_64
julia --threads 1 support/timeEstimator.jl

#mpiexecjl -n 4 julia --threads 30 generateData.jl

#git add --all
#git commit -am "job done"
#git push

wait