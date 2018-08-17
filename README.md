# Pipeline-practice #

## Overview ##
This is a practice pipeline for sample WES data that was already analyzed by CRI. The purpose is to develop usable and portable code for preprocessing and variant calling and to test the stability of the code before using it in primary analysis. The ereprocessing pipeline follows very closely to GATK best practices, mostly using legacy tutorials from GATK 3.x and modifying where necessary for GATK 4.x

## Inputs ##
This pipeline assumes paired-end WES input in FASTQ format (file_1 and file_2). 

## Usage ##
Each step in the pipeline has a PBS script for batch submission to Gardner HPC at UChicago which points to a single BASH shell. To use the scripts, ensure the scripts are in the correct directories and that all necessary scripts have had permissions changed to make executable:
`chmod +x <SCRIPTNAME>`

Go to the directory holding the executable PBS script and submit the job: 
`qsub <SCRIPTNAME.pbs>`

## Tool versions ##
Java v1.8
Picard tools version v2.8.1
GATK v4.0.6.0
BWA 0.7.17

## Reference sequences ##

### BWA ###
Reference alignment genome is version prepared by Heng Li (BWA author): GRCh38 with viral decoys and no alt contigs
Original filename: GRCh38_full_plus_hs38d1_analysis_set
See http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use for more info about the reference genome

### BaseRecalibrater ###
dbSNP release 151: https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi
Mills_and_1000G_gold_standard.indels.hg38.vcf and 000G_phase1.snps.high_confidence.hg38.vcf from GATK Resource Bundle (hg38)



