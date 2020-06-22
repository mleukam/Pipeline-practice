## This script runs a two-stage process called Base Quality Score Recalibration (BQSR).
## Specifically, it produces a recalibrated table with the BaseRecalibrator tool
## It then outputs a recalibrated BAM or CRAM file.
## Input is aligned BAM files with duplicates marked
## Outputs are a recalibration table and a recalibrated BAM
## Usage = ./mlrecalbases.sh inside PBS script
## Built specifically for WES data that was transferred 8/2/18
## Designed for batch submission to Gardner HPC at UChicago
## First line = shebang to specify interpretor (bash)
## Before using, use chmod to make executable

#!/bin/bash

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t'

## Load compiler
module load java-jdk/1.8.0_92

## Load necessary modules
## NB: the path to GATK is automatically generated as an environment variable
## Location of GenomeAnalysis.jar = ${GATK}
module load gatk/4.0.6.0

## Navigate to directory with input aligned and merged BAM files
cd /scratch/mleukam/james_wes/

## Define the input files as an array
RGBAMLIST=($(ls *_addRG.bam))

## Pull the sample name from the input file names and make new array
SMLIST=(${RGBAMLIST[*]%_*})

## Analyze patterns of covariation in the sequence dataset
for SAMPLE in ${SMLIST[*]};
do 
	java -Xmx16G -jar ${GATK} BaseRecalibrator \
    	-R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    	-I ${SAMPLE}_addRG.bam \
   		--known-sites /group/kline-lab/ref/dbsnp_151.vcf \
    	--known-sites /group/kline-lab/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    	--known-sites /group/kline-lab/ref/1000G_phase1.snps.high_confidence.hg38.vcf \
    	-O ${SAMPLE}_recal_data.table;
## Apply the recalibration to the sequence data
	java -Xmx16G -jar ${GATK} ApplyBQSR \
		-R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    	-I ${SAMPLE}_addRG.bam \
    	--bqsr-recal-file ${SAMPLE}_recal_data.table \
    	-O ${SAMPLE}_recal_reads.bam;
 done
