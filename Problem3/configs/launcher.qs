#!/bin/bash
#SBATCH --job-name=ICE_DYNAMICS_TEST
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --partition=x20core
#SBATCH --reservation=sirius

mpirun ./DYNAMICS_TESTS ./config.json > ./output.txt
