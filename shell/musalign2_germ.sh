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

## PBS script must call for 8G of mem and 12 threads

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

# convert bam to unaligned bam file
# include necessary read information - this will need to be added by hand!
java -Xmx8G -jar ${PICARD} RevertSam \
    I=${inbam} \
    O=${sample}_unaligned.bam \
    SANITIZE=true \
    MAX_DISCARD_FRACTION=0.010 \
    ATTRIBUTE_TO_CLEAR=XT \
    ATTRIBUTE_TO_CLEAR=XN \
    ATTRIBUTE_TO_CLEAR=AS \
    ATTRIBUTE_TO_CLEAR=OC \
    ATTRIBUTE_TO_CLEAR=OP \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=true \
    REMOVE_DUPLICATE_INFORMATION=true \
    REMOVE_ALIGNMENT_INFORMATION=true || error_exit "conversion to uBAM failed"

## Mark Illumina adapters
## Input: unaligned BAM
## Output: marked Illumina adapters
java -Xmx8G -jar ${PICARD} MarkIlluminaAdapters \
I=${sample}_unaligned.bam \
O=${sample}_markilluminaadapters.bam \
M=${sample}_markilluminaadapters_metrics.txt \
TMP_DIR=${tmpdir} || error_exit "failed to mark Illumina adapters" 

## Save uBAM for later merger

## Initial alignment, merge with uBAM
## Input: marked Illumina adapters BAM
## Output: Piped (aligned and merged) BAM
java -Xmx8G -jar ${PICARD} SamToFastq \
I=${sample}_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=${tmpdir} | \
bwa mem -M -t 12 -p ${ref} /dev/stdin | \
java -Xmx8G -jar ${PICARD} MergeBamAlignment \
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
java -Xmx8G -jar ${PICARD} MarkDuplicates \
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
java -Xmx8G -jar ${GATK} BaseRecalibrator \
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