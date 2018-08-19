## Script to mark duplicates in merged "clean" BAM files
## Input is aligned BAM files after merger with uBAM (output of mlalign.sh)
## Output is marked BAM file
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
module load GATK/4.0.6.0

## Navigate to directory with input aligned and merged BAM files
cd /scratch/mleukam/james_wes/

## Define the input files as an array
RGBAMLIST=($(ls *_addRG.bam))

## Pull the sample name from the input file names and make new array
SMLIST=(${RGBAMLIST[*]%_*})

## Analyze patterns of covariation in the sequence dataset
for SAMPLE in ${SMLIST[*]} 
do 
	java -Xmx16G -jar ${GATK} BaseRecalibrator \
    -R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    -I ${SAMPLE}_addRG.bam \
    --known-sites /group/kline-lab/ref/dbsnp_151.vcf \
    --known-sites /group/kline-lab/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    --known-sites /group/kline-lab/ref/1000G_phase1.snps.high_confidence.hg38.vcf \
    -O ${SAMPLE}_recal_data.table 
done

## Do a second pass to analyze covariation remaining after recalibration
for SAMPLE in ${SMLIST[*]}
do
	java -Xmx16G -jar ${GATK} BaseRecalibrator \
    -R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    -I ${SAMPLE}_addRG.bam \
    --known-sites /group/kline-lab/ref/dbsnp_151.vcf \
    --known-sites /group/kline-lab/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    --known-sites /group/kline-lab/ref/1000G_phase1.snps.high_confidence.hg38.vcf \
    -BQSR ${SAMPLE}_recal_data.table \
    -o ${SAMPLE}_post_recal_data.table 
done

## Generate before/after plots
for SAMPLE in ${SMLIST[*]}
do
	java -Xmx16G -jar ${GATK} AnalyzeCovariates \
	-R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    -before ${SAMPLE}_recal_data.table \
    -after ${SAMPLE}_post_recal_data.table \
    -plots ${SAMPLE}_recalibration_plots.pdf
done

## Apply the recalibration to the sequence data
for SAMPLE in ${SMLIST[*]}
do
	java -Xmx16G -jar ${GATK} PrintReads \
	-R /group/kline-lab/ref/GRCh38_full_plus_decoy.fa \
    -I ${SAMPLE}_addRG.bam \
    -BQSR ${SAMPLE}_recal_data.table \
    -O ${SAMPLE}_recal_reads.bam 
done
