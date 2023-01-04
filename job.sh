#!/bin/bash

#SBATCH -J 128-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=3-0:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 128 t 0.5 0.6 0.7 0.8 0.9 a 2.0 NR CLS LLS s 0 1000 

wait
