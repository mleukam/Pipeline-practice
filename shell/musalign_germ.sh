#!/bin/bash

##########
# README #
##########

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
# 1. Single read fastq file containing all the reads
# 2. Reference sequence for alignment -- indexed for bwa and with picard dictionary (see musprep script)
# NB: this script expects the reference genome, index files and known variants in working directory

# The following outputs are obtained:
# 1. ${SAMPLE}_aligned.bam --> ready for variant calling
# 2. ${SAMPLE}_unaligned.bam --> bam file containing raw reads for archival purposes

#########
# SETUP #
#########

# set script to fail if any command, variable, or output fails
set -euo pipefail

# set IFS to split only on newline and tab
IFS=$'\n\t' 

# load compilers
module load java-jdk/1.8.0_92 
module load gcc/6.2.0
 
# load modules
# NB: the path to picard.jar and gatk.jar is loaded with the module
# location of picard.jar = ${PICARD}
# location of gatk.jar = ${GATK}
module load fastqc/0.11.5
module load picard/2.8.1
module load bwa/0.7.17
module load samtools/1.6.0

# navigate to working directory
cd /scratch/mleukam/mouse

####################
# DEFINE VARIABLES #
####################

# sample name (do not include .fq suffix)
SAMPLE=BALB_cJ
# Note: Fastq file containing all reads for sample: A20.fq in working directory

# name of reference genome file (do not include .fa suffix)
REF=genome
# Note: mm10, patch 6 from Ensembl via igenomes; Ensembl contigs are named: 1, 2, etc
# This script expects the reference genome, picard dict, and bwa index files in the working directory

# path to fastqc output directory
FQDIR=/scratch/mleukam/mouse/fastqc

# path to temporary working directory
TMPDIR=/scratch/mleukam/temp

##########
# FASTQC #
##########

echo starting run for ${SAMPLE}

# get fastqc report on raw sequences before proceeding
# note: fastqc won't create the output directory; has to be done beforehand
# creates zip file and html file in the output directory
fastqc -o ${FQDIR} ${SAMPLE}.fq

echo fastqc completed for ${SAMPLE}

#############
# MAKE UBAM #
#############

# convert bam to unaligned bam file
# include necessary read information - this will need to be added by hand!
java -Xmx32G -jar ${PICARD} RevertSam \
    I=${SAMPLE}.bam \
    O=${SAMPLE}_unaligned.bam \
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
    REMOVE_ALIGNMENT_INFORMATION=true

echo ${SAMPLE} converted to unaligned BAM

# mark Illumina adapters
java -Xmx32G -jar ${PICARD} MarkIlluminaAdapters \
I=${SAMPLE}_unaligned.bam \
O=${SAMPLE}_markilluminaadapters.bam \
M=${SAMPLE}_markilluminaadapters_metrics.txt \
TMP_DIR=${TMPDIR}

echo illumina adapters marked

###################
# ALIGN SEQUENCES #
###################

# revert BAM file temporarily back to fastq
java -Xmx32G -jar ${PICARD} SamToFastq \
I=${SAMPLE}_markilluminaadapters.bam \
FASTQ=${SAMPLE}_temp.fq \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=${TMPDIR}

# align with bwa mem in paired/interleaved mode
bwa mem -M -t 12 -p ${REF}.fa balbc.fq > balbc_aligned.sam

echo ${SAMPLE} is aligned

# filter out unmapped and multimapped reads
# these reads are not useful for variant calling and create downstream errors
samtools view -F 4 -q 1 ${SAMPLE}_temp.sam > ${SAMPLE}_temp.aligned.filtered.sam

# merge in picard sequence dictionary created in setup
# create a new file (unsortedtemp.sam) that has both the dictionary and the aligned reads.
# using this template: cat dictionary.sam > unsorted_file.sam && cat file.sam >> unsorted_file.sam
cat {REF}.dict > ${SAMPLE}_dict.sam && cat ${SAMPLE}_temp.aligned.filtered.sam >> ${SAMPLE}_dict.sam

# samtools sort creates downstream errors with picard tools
# only recommend sorting using picard tools
# sort aligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=${SAMPLE}_dict.sam \
O=${SAMPLE}_temp.aligned.query.bam \
SORT_ORDER=queryname

# sort unaligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=${SAMPLE}_unaligned.bam \
O=${SAMPLE}_unaligned.query.bam \
SORT_ORDER=queryname

# merge aligned bam with ubam to restore headers, quality and read group information
# merging also removes hardclips from BWA for discordance between best matching kmer and read
java -Xmx32G -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=${SAMPLE}_temp.aligned.query.bam \
UNMAPPED_BAM=${SAMPLE}_unaligned.query.bam \
OUTPUT=${SAMPLE}_merged.bam \
R=${REF}.fa \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=${TMPDIR}

# sort merged BAM by coordinate for later analysis
java -Xmx32G -jar ${PICARD} SortSam \
I=${SAMPLE}_merged.bam \
O=${SAMPLE}_merged.sorted.bam \
SORT_ORDER=coordinate

echo unmatched reads filtered, ${SAMPLE} is merged

# mark duplicates
java -Xmx32G -jar ${PICARD} MarkDuplicates \
INPUT=${SAMPLE}_merged.sorted.bam \
OUTPUT=${SAMPLE}_aligned.bam \
METRICS_FILE=${SAMPLE}_markduplicates_metrics.txt \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
CREATE_INDEX=true \
TMP_DIR=${TMPDIR}

echo duplicates marked

# clean up
rm ${SAMPLE}_markduplicates.bam
rm ${SAMPLE}_temp.fq
rm ${SAMPLE}_temp.sai
rm ${SAMPLE}_temp.sam
rm ${SAMPLE}_temp.aligned.bam
rm ${SAMPLE}_filteredtemp.sam
rm ${SAMPLE}_unsortedtemp.sam
rm ${SAMPLE}_aligned.filtered.sam
rm ${SAMPLE}_unaligned.query.bam
rm ${SAMPLE}_aligned.query.bam
rm ${SAMPLE}_dict.sam
rm ${SAMPLE}_markilluminaadapters.bam

# end the run
echo alignment script completed
echo ${SAMPLE} is ready for analysis
