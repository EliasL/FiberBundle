#!/bin/bash

#SBATCH -J 512-200
#SBATCH -p porelab
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --nice=50000
#SBATCH --time=9-24:15:0

ml eb
ml Julia/1.7.2-linux-x86_64
julia --threads 50 generateData.jl L 512 t 0.05 0.1 0.15 0.2 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 a 2.0 NR LLS CLS s 1 200 dist ConstantAverageUniform 

wait
