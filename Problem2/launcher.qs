#!/bin/bash
#SBATCH --job-name=Transport_test
#SBATCH --ntasks=32
#SBATCH --time=12:00:00
#SBATCH --partition=x20core
#SBATCH --reservation=sirius

mpirun ./run ./config.json > ./output.txt
