I processed few fastq samples (script is still running but will kill it soon) on circe from the dataset you share here using the attached scripts (scr_count_aggr_all.sh and run.scr_sc_all.sh). 

Synapse | Sage Bionetworks
Sage Bionetworks, info@sagebase.org

Synapse is a platform for supporting scientific collaborations centered around shared biomedical data sets. Our...


Few things to note:
1. each sample was processed via 4 lanes (L1, L2, L3 and L3 in the sample names) and has three fastq.gz files (I1, R1 and R2 in the sample names) and there are total 48 samples (DI17-XXXX_SX in the sample names).
  Exp. For sample D17-8753_S1, the 7 files are: fastqs/D17-8753_S1_L001_I1_001.fastq.gz  fastqs/D17-8753_S1_L002_R1_001.fastq.gz  fastqs/D17-8753_S1_L003_R2_001.fastq.gz
fastqs/D17-8753_S1_L001_R1_001.fastq.gz  fastqs/D17-8753_S1_L002_R2_001.fastq.gz  fastqs/D17-8753_S1_L004_I1_001.fastq.gz
fastqs/D17-8753_S1_L001_R2_001.fastq.gz  fastqs/D17-8753_S1_L003_I1_001.fastq.gz  fastqs/D17-8753_S1_L004_R1_001.fastq.gz
fastqs/D17-8753_S1_L002_I1_001.fastq.gz  fastqs/D17-8753_S1_L003_R1_001.fastq.gz  fastqs/D17-8753_S1_L004_R2_001.fastq.gz

2. script gets all sample names D17-XXXX_SX and feeds it to the cellranger count command one by one in a loop.
3. I have used statndard commad option, but please check the paper and use exactly same options (if they have used different that standard options) in the cellranger count command in the script.
4. After completing the cellranger count, code call aggr to aggregate results for all samples processed in the step above. Please check the paper and use exactly same options (if they have used different that standard options) in the cellranger count command in the script.

You might also want to change ntasks and mem-per-cpu in the run (scr script) script to speed up calculations.

I have followed this blog and this website to make the script, please refer to these resources (any other that you might find on the web) if you hit any issues and let me know if you have any questions 

Getting started with Cell Ranger - Dave Tang's blog
Davo

Cell Ranger is a set of analysis pipelines that process Chromium single cell 3' RNA-seq data. The pipelines proc...




Regards
Dr. Syed Islamuddin Shah