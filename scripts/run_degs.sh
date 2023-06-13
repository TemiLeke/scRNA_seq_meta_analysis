#!/bin/bash
#SBATCH --job-name=running_degs
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=96:00:00
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tadeoye@usf.edu

module purge
module add apps/R/4.0.2-el7-gcc-openblas
module add mpi/openmpi/3.1.6

mpirun Rscript applyMAST.R