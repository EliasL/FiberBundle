#!/bin/bash

#SBATCH -J 128-1000
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --time=3-0:43:7

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 128 t 0.0 a 2.0 NR ELS s 0 1000 

wait
