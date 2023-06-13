#!/bin/bash

export PATH="/work/t/tadeoye/scRNAseq_AD_meta_analysis/src/cell_ranger/cellranger-7.0.1:$PATH"

echo PLEASE MAKE SURE THAT YOU ARE USING APPROPRIATE SETTINGS FOR THIS STEP BASED ON THE REFERENCE PAPER FOR WHICH DATA IS BEING PROCESSED
SECONDS=0

# change working directory

cd /work/t/tadeoye/scRNAseq_AD_meta_analysis/data/mathys_pfc/raw/

samples=$(ls fastq/*.gz|cut -d'_' -f1|cut -d'/' -f2|uniq)

#also create csv file for aggr run
#First check if file exits

[-f mathys_pfc_raw.csv] && rm mathys_pfc_info.csv
echo sample_id,molecule_h5 > mathys_pfc_info.csv

echo starting to process all samples

i=1
for sample in $samples
do
	echo processing sample num $i and sample $sample

	cellranger count --id=mathys_pfc_$sample \
		 --transcriptome=/work/t/tadeoye/scRNAseq_AD_meta_analysis/data/reference_genome/refdata-gex-GRCh38-2020-A/ \
        	 --fastqs=/work/t/tadeoye/scRNAseq_AD_meta_analysis/data/mathys_pfc/raw/fastq/ \
                 --sample=$sample \ 
                 --localcores=16 \
                 --localmem=16 \
				 --chemistry=SC3Pv2 
				 
	echo $sample, mathys_pfc_$sample/outs/molecule_info.h5 >> mathys_pfc_info.csv
	((i=i+1))
done
echo DONE count step for all samples, NOW STARTING AGGREGATE step

echo PLEASE MAKE SURE THAT YOU ARE USING APPROPRIATE SETTINGS FOR THIS STEP BASED ON THE REFERENCE PAPER FOR WHICH DATA IS BEING PROCESSED

cellranger aggr --id=mathys_pfc_aggr \
                --csv=mathys_pfc_info.csv

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

