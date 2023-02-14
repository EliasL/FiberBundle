#!/bin/bash

#SBATCH -J 8-40000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=4-2:68:10

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 8 t 0.0 0.1 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.5 0.6 0.7 0.8 0.9 a 2.0 NR LLS CLS s 10000 40000 

wait
