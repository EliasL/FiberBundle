#!/bin/bash

#SBATCH -J 512-100
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --exclusive=user
#SBATCH --time=8-28:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 8 16 32 64 128 256 512 t 0.0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.5 0.6 0.7 0.8 0.9 a 2.0 NR CLS LLS s 0 100 

wait
