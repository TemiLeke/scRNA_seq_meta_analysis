#!/bin/bash
#SBATCH --job-name=mathys_fastq        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tadeoye@usf.edu

module purge
module load apps/python/3.8.5

python get_synapse_data.py
