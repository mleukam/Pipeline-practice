#!/bin/bash

###########################################
## PRE-PROCESSING AND ALIGNING HUMAN WES ##
###########################################

## PBS script must call for 8G of mem and 12 threads


## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Assign variables given as arguments
progname=$(basename $0)
sample=$1
fq1=$2
fq2=$3

## Assign memory
mem=Xmx16G

## Assign working and temp directories 
wkdir=/scratch/mleukam/wes_controls
tmpdir=/scratch/mleukam/temp

## Assign references
ref=/gpfs/data/kline-lab/ref/hg19/hg19/ucsc.hg19.fasta
knownsites1=/gpfs/data/kline-lab/ref/hg19/hg19/dbsnp_138.hg19.vcf
knownsites2=/gpfs/data/kline-lab/ref/hg19/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf
knownsites3=/gpfs/data/kline-lab/ref/hg19/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

## Assign read group information given as arguments
instrument=$4
readgroup=$5

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

## navigate to file containing WES files
cd ${wkdir}

## Convert fastq to uBAM
## Inputs: fq1 and fq2 from arguments
## Output: unaligned BAM
java -${mem} -jar ${PICARD} FastqToSam \
F1=${fq1} \
F2=${fq2} \
O=${sample}_unaligned.bam \
READ_GROUP_NAME=${readgroup} \
SAMPLE_NAME=${sample} \
LIBRARY_NAME=library1 \
PLATFORM_UNIT=${instrument} \
PLATFORM=illumina \
SEQUENCING_CENTER=Theragenetex \
RUN_DATE=2018-09-30 || error_exit "failed to convert fastq to uBAM with readgroup info" 

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
--known-sites ${knownsites3} \
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





