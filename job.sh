#!/bin/bash

#SBATCH -J 256-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --nice
#SBATCH --time=3-7:69:21

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 50 generateData.jl L 256 t 0.9 a 2.0 NR LLS CLS s 0 1000 dist ConstantAverageUniform 

wait
