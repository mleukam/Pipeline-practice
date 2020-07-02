#!/bin/bash

##########
# README #
##########

# align short (50bp) *single end* whole genome reads from mouse tumor
# mark illumina adapters, filter unmapped reads and multimapped reads, index and sort
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

####################################
# BASE QUALITY SCORE RECALIBRATION #
####################################

# generate bqsr table
java -Xmx32G -jar ${GATK} BaseRecalibrator \
-R ${REF}.fa \
--known-sites ${SNP} \
--known-sites ${INDELS} \
-I ${SAMPLE}_markduplicates.bam \
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
