#!/bin/bash

#SBATCH -J 256-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=3-0:15:37

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 256 t 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 a 2.0 NR CLS LLS s 0 1000 

wait
