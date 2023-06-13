#!/bin/bash

cd /work/t/tadeoye/scRNAseq_AD_meta_analysis/data/mathys_pfc/raw/fastq/

echo starting to update all samples

i=1
for file in /work/t/tadeoye/scRNAseq_AD_meta_analysis/data/mathys_pfc/raw/fastq/*
do 
    echo y | /work/t/tadeoye/scRNAseq_AD_meta_analysis/src/txg-linux-v1.3.0/txg fastqs upload --project-id 0PNKznfzQE2GaFMp03fNkQ $file

	((i=i+1))
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
