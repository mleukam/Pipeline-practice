## Script to make GATK-usable index files from reference lists
## Inputs are reference VCF files (see repo README for details)
## Outputs are .tbi files
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

## Navigate to reference folder
cd /group/kline-lab/ref

## This script will not use a for loop due to the one-time nature of its use
## For future use this will have to be customized again
java -Xmx8G -jar ${GATK} IndexFeatureFile -F dbsnp_151.vcf
java -Xmx8G -jar ${GATK} IndexFeatureFile -F 1000G_phase1.snps.high_confidence.hg38.vcf
java -Xmx8G -jar ${GATK} IndexFeatureFile -F Mills_and_1000G_gold_standard.indels.hg38.vcf
