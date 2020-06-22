#!/bin/bash

#####################
## VARIANT CALLING ##
#####################

## README

# This script is part of a custom pipeline for:
# calling variants from mouse WGS from Balbc/J background
# intended to run on gardner HPC with PBS wrapper

# The following inputs are required:
# 1. Aligned tumor BAM file with BQSR done for specific mouse strain
# 2. Aligned normal BAM file with BQSR done for specific mouse strain
# 3. Known indels from mouse genome project for filtering -- indexed with GATK IndexFeatureFile
# 4. Known snps from mouse genome project for filtering -- indexed with GATK IndexFeatureFile

# The following outputs are obtained:
# 1. VCF of variants, ready for annotation and further filtering
# 2. BAM of assembled haplotypes

## PBS script must call for 16G of mem and 1 thread

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Assign variables given as arguments
progname=$(basename $0)
sample=$1
norm_name=$2

## Assign additional variables including working directory containing files
wkdir=/scratch/mleukam/mouse/
ref=/scratch/mleukam/mouse/genome.fa
knownsites1=/scratch/mleukam/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf
knownsites2=/scratch/mleukam/mouse/mgp.v5.merged.snps_all.dbSNP142.vcf
tmpdir=/scratch/mleukam/temp
tumor=A20_bqsr.bam
normal=BALB_cJ_bqsr.bam

## Error handling
error_exit()
{
#	----------------------------------------------------------------
#	Function for exit due to fatal program error
#		Accepts 1 argument:
#			string containing descriptive error message
#	----------------------------------------------------------------
	echo "${progname}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

## Load compiler
module load java-jdk/1.8.0_92

## Load modules
module load gatk/4.0.6.0

## navigate to file containing bam files
cd $wkdir || error_exit "unable to navigate to working directory"

## generate VCF and bamout with mutect2
## can use four threads
java -Xmx16G -jar ${GATK} Mutect2 \
-R ${ref} \
-I ${tumor} \
-I ${normal} \
--tumor-sample ${sample} \
--normal-sample ${norm_name} \
--germline-resource ${knownsites2} \
--output ${sample}.vcf.gz \
-bamout tumor_normal_${sample}.bam \
--TMP_DIR ${tmpdir}

# Exit
echo "Variant calling completed. VCF is ready for further analysis"
exit 0
