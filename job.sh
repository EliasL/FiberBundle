#!/bin/bash

#SBATCH -J 128-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --nice=50000
#SBATCH --time=9-24:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 64 generateData.jl L 128 t 1.4 a 2.0 NR LLS CLS s 0 1000 dist Weibull 

wait
