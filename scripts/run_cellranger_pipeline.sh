#!/bin/bash
#SBATCH --job-name=get_mathys_pfc_counts
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=output.%j.get_mathys_pfc_counts
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tadeoye@usf.edu

bash cellranger_count_aggr_all.sh
