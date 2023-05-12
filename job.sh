#!/bin/bash

#SBATCH -J 512-200
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --nice=50000
#SBATCH --time=9-24:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 2 generateData.jl L 512 t 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 a 2.0 NR LLS CLS s 0 200 dist ConstantAverageUniform 

wait
