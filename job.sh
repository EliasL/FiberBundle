#!/bin/bash

#SBATCH -J 256-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=3-0:15:27

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 256 t 0.8 0.9 a 2.0 NR LLS CLS s 0 1000 

wait
