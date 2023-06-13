#!/bin/bash 


#!/bin/bash

export PATH="/work/t/tadeoye/scRNAseq_AD_meta_analysis/src/FastQC:$PATH"

SECONDS=0

# change working directory

cd /work/t/tadeoye/scRNAseq_AD_meta_analysis/data/mathys_pfc/raw/bam/

samples=$(ls */ | cut -d'/' -f 1)

mkdir /work/t/tadeoye/scRNAseq_AD_meta_analysis/results/mathys_pfc/fastqc/

echo starting to process all samples

i=1

for sample in $samples
do
	echo processing sample num $sample
	
    fastqc "$sample/possorted_genome_bam.bam" --outdir=/work/t/tadeoye/scRNAseq_AD_meta_analysis/results/mathys_pfc/fastqc/

	mv /work/t/tadeoye/scRNAseq_AD_meta_analysis/results/mathys_pfc/fastqc/possorted_genome_bam_fastqc.html "/work/t/tadeoye/scRNAseq_AD_meta_analysis/results/mathys_pfc/fastqc/sample_$sample_qc.html"
				 
	((i=i+1))
done

echo DONE QC step for all samples

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

