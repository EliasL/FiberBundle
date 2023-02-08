#!/bin/bash

#SBATCH -J 32-10000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=3-0:15:1

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 16 32 t 0.0 a 2.0 NR LLS CLS s 0 10000 

wait
