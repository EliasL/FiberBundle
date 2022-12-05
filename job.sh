#!/bin/bash

#SBATCH -J 256-100
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --exclusive=user
#SBATCH --time=0-5:53:1

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 8 16 32 64 128 256 t 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 NR SNR UNR s 0 100 

wait
