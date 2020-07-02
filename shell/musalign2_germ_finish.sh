#!/bin/bash

###########################################
## PRE-PROCESSING AND ALIGNING MOUSE WGS ##
###########################################

## README

# purpose: realign 100bp *paired-end* whole genome reads from balb/c sanger mouse genome project
# in addition to aligning, script will
# -- mark illumina adapters
# -- filter unmapped reads and multimapped reads
# -- index results of alignment and sort by coordinate
# -- merge results with unaligned bam to restore read group info and correct hard clipping
# -- mark duplicates after alignment
# intended to run on gardner HPC with PBS wrapper
# before using, make executable with chmod
# before starting, update the readgroup information (not currently encoded as a variable)

# The following inputs are required:
# 1. BAM file containing all the reads
# 2. Reference sequence for alignment -- indexed for bwa and with picard dictionary (see musprep script)
# 3. Known indels from mouse genome project for BQSR
# 4. Known snps from mouse genome project for BQSR

# The following outputs are obtained:
# 1. ${sample}_aligned.bam --> ready for variant calling

## PBS script must call for 16G of mem and 1 thread

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Assign variables given as arguments
progname=$(basename $0)
sample=$1
inbam=$2

## Assign additional variables including working directory containing files
wkdir=/scratch/mleukam/mouse/
ref=/scratch/mleukam/mouse/genome.fa
# Note: mm10, patch 6 from Ensembl via igenomes; Ensembl contigs are named: 1, 2, etc
# This script expects the reference genome, picard dict, and bwa index files in the working directory
knownsites1=/scratch/mleukam/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf
knownsites2=/scratch/mleukam/mouse/mgp.v5.merged.snps_all.dbSNP142.vcf
fqdir=/scratch/mleukam/mouse/fastqc
tmpdir=/scratch/mleukam/temp

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

## Load compilers
module load java-jdk/1.8.0_92
module load gcc/6.2.0

## Load modules
module load picard/2.8.1
module load gatk/4.0.6.0
module load bwa/0.7.17
module load fastqc/0.11.5

## navigate to file containing fastq files
cd ${wkdir} || error_exit "unable to navigate to working directory"


## Base quality score re-alignment
## Input: marked duplicates BAM
## Output: BQSR table and recalibrated BAM (bqsr.bam)
java -Xmx16G -jar ${GATK} BaseRecalibrator \
-R ${ref} \
-I ${sample}_markduplicates.bam \
--known-sites ${knownsites1} \
--known-sites ${knownsites2} \
-O ${sample}_bqsr.table || error_exit "failed to create BQSR table"

java -Xmx8G -jar ${GATK} ApplyBQSR \
-R ${ref} \
-I ${sample}_markduplicates.bam \
--bqsr-recal-file ${sample}_bqsr.table \
-O ${sample}_bqsr.bam || error_exit "failed to apply BQSR table"

# Clean up: deleme marked duplicates BAM
rm ${sample}_markduplicates.bam

# Exit
echo "Alignment and BQSR completed. ${sample}_bqsr.bam is ready for further analysis"
exit 0
