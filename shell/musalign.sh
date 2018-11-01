#!/bin/bash

##########
# README #
##########

# align short (50bp) *single end* whole genome reads from mouse tumor
# in addition to aligning, script will
# -- mark illumina adapters
# -- filter unmapped reads and multimapped reads
# -- index results of alignment and sort by coordinate
# -- merge results with unaligned bam to restore read group info and correct hard clipping
# -- base quality score recalibration after alignment
# intended to run on gardner HPC with PBS wrapper
# before using, make executable with chmod
# before starting, update the readgroup information (not currently encoded as a variable)

# The following inputs are required:
# 1. Single read fastq file containing all the reads
# 2. Reference sequence for alignment -- indexed for bwa and with picard dictionary (see musprep script)
# 3. Known germline SNPs and Indels for BQSR
# NB: this script expects the reference genome, index files and known variants in working directory

# The following outputs are obtained:
# 1. ${SAMPLE}_bqsr.bam --> ready for variant calling
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
module load gatk/4.0.6.0

# navigate to working directory
cd /scratch/mleukam/mouse

####################
# DEFINE VARIABLES #
####################

# sample name (do not include .fq suffix)
SAMPLE=A20
# Note: Fastq file containing all reads for sample: A20.fq in working directory

# name of reference genome file (do not include .fa suffix)
REF=genome
# Note: mm10, patch 6 from Ensembl via igenomes; Ensembl contigs are named: 1, 2, etc
# This script expects the reference genome, picard dict, and bwa index files in the working directory

# name of known indels (include suffix)
INDELS=mgp.v5.merged.indels.dbSNP142.normed.vcf
# Note: known germline SNPs in all 18 sequenced mouse strains. Uses Ensembl contigs.

# name of known SNPs (include suffix)
SNP=mgp.v5.merged.snps_all.dbSNP142.vcf 
# Note: known germline indels in all 18 sequenced mouse strains. Uses Ensembl contigs.

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

# convert fastq to unaligned bam file
# include necessary read information - this will need to be added by hand!
java -Xmx32G -jar ${PICARD} FastqToSam \
FASTQ=${SAMPLE}.fq \
O=${SAMPLE}_unaligned.bam \
READ_GROUP_NAME=HW5FWBBXX.6 \
SAMPLE_NAME=${SAMPLE} \
LIBRARY_NAME=coreWGS \
PLATFORM_UNIT=K00242 \
PLATFORM=illumina \
SEQUENCING_CENTER=UCHICAGOCORE \
RUN_DATE=2018-09-30T00:00:00-0400

echo ${SAMPLE}.fq converted to unaligned BAM

# mark Illumina adapters
java -Xmx32G -jar ${PICARD} MarkIlluminaAdapters \
I=${SAMPLE}_unaligned.bam \
O=${SAMPLE}_markilluminaadapters.bam \
M=${SAMPLE}_markilluminaadapters_metrics.txt \
TMP_DIR=/scratch/mleukam/${TMPDIR}

echo illumina adapters marked

###################
# ALIGN SEQUENCES #
###################

# revert BAM file temporarily back to fastq
java -Xmx32G -jar ${PICARD} SamToFastq \
I=${SAMPLE}_markilluminaadapters.bam \
FASTQ=${SAMPLE}_temp.fq \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/scratch/mleukam/${TMPDIR}

# align sequences with BWA aln+samse (best for single end reads <70bp)
# note that t flag in bwa aln is set to 12, max is 28 threads/cores on single node in Gardner
# will only work if PBS script asks for 12 ppn
bwa aln -t 12 ${REF}.fa \
${SAMPLE}_temp.fq > ${SAMPLE}_temp.sai

bwa samse ${REF}.fa \
${SAMPLE}_temp.sai \
${SAMPLE}_temp.fq > ${SAMPLE}_temp.sam

echo ${SAMPLE} is aligned

# filter out unmapped and multimapped reads
# these reads are not useful for variant calling and create downstream errors
samtools view -F 4 -q 1 ${SAMPLE}_temp.sam > ${SAMPLE}_aligned.filtered.sam

# merge in picard sequence dictionary created in setup
# create a new file (unsortedtemp.sam) that has both the dictionary and the aligned reads.
# using this template: cat dictionary.sam > unsorted_file.sam && cat file.sam >> unsorted_file.sam
cat {REF}.dict > ${SAMPLE}_dict.sam && cat ${SAMPLE}_aligned.filtered.sam >> ${SAMPLE}_dict.sam

# samtools sort creates downstream errors with picard tools
# only recommend sorting using picard tools
# sort aligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=${SAMPLE}_dict.sam \
O=${SAMPLE}_aligned.query.bam \
SORT_ORDER=queryname

# sort unaligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=${SAMPLE}_unaligned.bam \
O=${SAMPLE}_unaligned.query.bam \
SORT_ORDER=queryname

# merge aligned bam with ubam to restore headers, quality and read group information
# merging also removes hardclips from BWA for discordance between best matching kmer and read
java -Xmx32G -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=${SAMPLE}_aligned.query.bam \
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
I=${SAMPLE}_merged.sam \
O=${SAMPLE}_merged.sorted.bam \
SORT_ORDER=coordinate

echo unmatched reads filtered, ${SAMPLE} is merged

####################################
# BASE QUALITY SCORE RECALIBRATION #
####################################

# mark duplicates
java -Xmx32G -jar ${PICARD} MarkDuplicates \
INPUT=${SAMPLE}_merged.sorted.bam \
OUTPUT=${SAMPLE}_markduplicates.bam \
METRICS_FILE=${SAMPLE}_markduplicates_metrics.txt \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
CREATE_INDEX=true \
TMP_DIR=${TMPDIR}

# generate bqsr table
java -Xmx32G -jar ${GATK} BaseRecalibrator \
-R ${REF}.fa \
-I ${SAMPLE}_markduplicates.bam \
--known-sites ${SNP} \
--known-sites ${INDELS} \
-O ${SAMPLE}_bqsr.table

# apply bqsr table
java -Xmx32G -jar ${GATK} ApplyBQSR \
-R ${REF}.fa \
-I ${SAMPLE}_markduplicates.bam \
--bqsr-recal-file ${SAMPLE}_bqsr.table \
-O ${SAMPLE}_bqsr.bam

echo BQSR completed for ${SAMPLE}

# clean up
rm ${SAMPLE}_markduplicates.bam
rm ${SAMPLE}_temp.fq
rm ${SAMPLE}_temp.sai
rm ${SAMPLE}_temp.sam
rm ${SAMPLE}_aligned.bam
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
