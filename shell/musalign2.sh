#!/bin/bash

###########################################
## PRE-PROCESSING AND ALIGNING MOUSE WGS ##
###########################################

## README

# This script is a custom pipeline for:
# align short (50bp) SINGLE END whole genome reads from mouse tumor
# in addition to aligning, script will
# -- mark illumina adapters
# -- filter unmapped reads and multimapped reads
# -- index results of alignment and sort by coordinate
# -- merge results with unaligned bam to restore read group info and correct hard clipping
# -- remove duplicates after alignment
# intended to run on gardner HPC with PBS wrapper
# before using, make executable with chmod
# before starting, update the readgroup information (not currently encoded as a variable)

# The following inputs are required:
# 1. Single read fastq file containing all the reads
# 2. Reference sequence for alignment -- indexed for bwa and with picard dictionary (see musprep script)
# NB: this script expects the reference genome, index files and known variants in working directory

# The following outputs are obtained:
# 1. ${sample}_aligned.bam --> ready for variant calling

## PBS script must call for 8G of mem and 12 threads

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Assign variables given as arguments
progname=$(basename $0)
sample=$1
fq=$2

## Assign additional variables including working directory containing files
wkdir=/scratch/mleukam/mouse/
ref=/scratch/mleukam/mouse/genome.fa
knownsites1=/scratch/mleukam/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf
knownsites2=/scratch/mleukam/mouse/mgp.v5.merged.snps_all.dbSNP142.vcf
fqdir=/scratch/mleukam/mouse/fastqc
tmpdir=/scratch/mleukam/temp
mem=Xmx16G

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

## get fastq quality report
fastqc -o ${fqdir} ${sample}.fq || error_exit "fastQC failed"
echo fastqc completed for ${sample}

## Convert fastq to uBAM
## Inputs: fq file from arguments
## Output: unaligned BAM
## NB: include necessary read information - this will need to be added by hand!
java -${mem} -jar ${PICARD} FastqToSam \
FASTQ=${sample}.fq \
O=${sample}_unaligned.bam \
READ_GROUP_NAME=HW5FWBBXX.6 \
SAMPLE_NAME=${sample} \
LIBRARY_NAME=coreWGS \
PLATFORM_UNIT=K00242 \
PLATFORM=illumina \
SEQUENCING_CENTER=UCHICAGOCORE \
RUN_DATE=2018-09-30T00:00:00-0400 || error_exit "fastq conversion to BAM failed" 

## Mark Illumina adapters
## Input: unaligned BAM
## Output: marked Illumina adapters
java -${mem} -jar ${PICARD} MarkIlluminaAdapters \
I=${sample}_unaligned.bam \
O=${sample}_markilluminaadapters.bam \
M=${sample}_markilluminaadapters_metrics.txt \
TMP_DIR=${tmpdir} || error_exit "failed to mark Illumina adapters" 

## Save uBAM for later merger

## Initial alignment, merge with uBAM
## Input: marked Illumina adapters BAM
## Output: Piped (aligned and merged) BAM
java -${mem} -jar ${PICARD} SamToFastq \
I=${sample}_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=${tmpdir} | \
bwa mem -M -t 12 -p ${ref} /dev/stdin | \
java -${mem} -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${sample}_unaligned.bam \
OUTPUT=${sample}_piped.bam \
R=${ref} \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=${tmpdir} || error_exit "alignment and/or merger failed"

## Clean up: delete marked Illumina adapter BAM
rm ${sample}_markilluminaadapters.bam

## Mark duplicates
## Input: piped BAM from alignment
## Output: marked duplicates BAM
java -${mem} -jar ${PICARD} MarkDuplicates \
INPUT=${sample}_piped.bam \
OUTPUT=${sample}_markduplicates.bam \
METRICS_FILE=${sample}_markduplicates_metrics.txt \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
CREATE_INDEX=true \
TMP_DIR=${tmpdir} || error_exit "mark duplicates failed" 

## Clean up: delete unmerged and merged BAMs
rm ${sample}_unaligned.bam
rm ${sample}_piped.bam

## Base quality score re-alignment
## Input: marked duplicates BAM
## Output: BQSR table and recalibrated BAM (bqsr.bam)
java -${mem} -jar ${GATK} BaseRecalibrator \
-R ${ref} \
-I ${sample}_markduplicates.bam \
--known-sites ${knownsites1} \
--known-sites ${knownsites2} \
-O ${sample}_bqsr.table || error_exit "failed to create BQSR table"

java -${mem} -jar ${GATK} ApplyBQSR \
-R ${ref} \
-I ${sample}_markduplicates.bam \
--bqsr-recal-file ${sample}_bqsr.table \
-O ${sample}_bqsr.bam || error_exit "failed to apply BQSR table"

# Clean up: deleme marked duplicates BAM
rm ${sample}_markduplicates.bam

# Exit
echo "Alignment and BQSR completed. ${sample}_bqsr.bam is ready for further analysis"
exit 0
