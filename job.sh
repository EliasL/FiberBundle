#!/bin/bash

#SBATCH -J 128-4000
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --exclusive=user
#SBATCH --time=0-0:27:60

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 128 t 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 a 2.0 NR CLS LLS s 0 4000 

wait
