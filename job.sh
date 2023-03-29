#!/bin/bash

#SBATCH -J 512-200
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --nice=50000
#SBATCH --time=9-24:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 50 generateData.jl L 512 t 1.0 a 2.0 NR LLS CLS s 1 200 dist ConstantAverageUniform 

wait
