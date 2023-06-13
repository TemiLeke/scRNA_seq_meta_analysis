#!/bin/bash
#SBATCH --job-name=sequence_quality
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tadeoye@usf.edu

bash sequence_quality.sh
