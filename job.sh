#!/bin/bash

#SBATCH -J 512-500
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --exclusive=user
#SBATCH --time=0-12:87:60

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 40 generateData.jl L 512 t 0.0 a 2.0 NR CLS LLS s 0 500 

wait
