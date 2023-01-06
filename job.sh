#!/bin/bash

#SBATCH -J 1024-30
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=9-24:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 1024 t 0.0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 a 2.0 NR CLS LLS s 0 30 

wait
