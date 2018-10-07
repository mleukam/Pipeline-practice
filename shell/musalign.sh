#!/bin/bash
#
# alignment script optimized for mouse whole genome
# this script is not yet optimized for general use 
# this script is customized for Sravya's WGS experiment 10/2018
# intended to run on gardner HPC with PBS wrapper
# before using, make executable with chmod
# 
#### INPUTS ####
# 1. Fastq file containing all tumor reads: A20.fq
# 2. Reference sequence for alingment: GRCm38p6_ref.fa
#
# set script to fail if any command, variable, or output fails
set -euo pipefail
#
# set IFS to split only on newline and tab
IFS=$'\n\t' 
#
# load compilers
module load java-jdk/1.8.0_92
module load gcc/6.2.0
# 
#### FASTQC #####
#
# load module
module load fastqc/0.11.5
#
# navigate to directory containing fastq files
cd /scratch/mleukam/mouse
#
# run fastqc on input sample
## note: fastqc won't create the output directory; has to be done beforehand
## creates zip file and html file in the output directory
fastqc -o /scratch/mleukam/mouse/fastqc A20.fq

#### CONVERT FASTQ TO UBAM ####
#
# Load necessary modules
# NB: the path to picard.jar is automatically generated as an environment variable when picard tool module is loaded
# location of picard.jar = ${PICARD}
module load picard/2.8.1
#
# no need to change directory
#
# run FastqToSam on suppled sample's fastq file
# returns unaligned BAM to the same directory
# -Xmx2G asks for 8GB RAM
# include necessary read information
java -Xmx16G -jar ${PICARD} FastqToSam \
FASTQ=A20.fq \
O=A20_unaligned.bam \
READ_GROUP_NAME=HW5FWBBXX.6 \
SAMPLE_NAME=A20 \
LIBRARY_NAME=coreWGS \
PLATFORM_UNIT=K00242 \
PLATFORM=illumina \
SEQUENCING_CENTER=UCHICAGOCORE \
RUN_DATE=2018-09-30T00:00:00-0400 \
#
#### MARK ILLUMINA ADAPTERS ####
#
# picard tools loaded above
# sample list created above
# run MarkIlluminaAdapters on uBAM created from provided sample
java -Xmx16G -jar ${PICARD} MarkIlluminaAdapters \
I=A20_unaligned.bam \
O=A20_markilluminaadapters.bam \
M=A20_markilluminaadapters_metrics.txt \
TMP_DIR=/scratch/mleukam/temp/
#
#### ALIGN SEQUENCES ####
#
# load necessary modules
module load bwa/0.7.17
module load samtools/1.6.0
# 
# generate index files from reference
bwa index -a bwtsw GRCm38p6_ref.fa 
samtools faidx GRCm38p6_ref.fa
java -jar ${PICARD} CreateSequenceDictionary \
    REFERENCE=GRCm38p6_ref.fa \
    OUTPUT=GRCm38p6_ref.dict
#
# loop to run pipeline on all of the called files in the directory
# note that t flag in bwa is set to 28 for number of cores in each Gardner node
# the reference genome used in this case is GRCm38 (mm10) patch 6 (most recent)
# index file for bwa downloaded with reference genome from NCBI
java -Xmx16G -jar ${PICARD} SamToFastq \
I=A20_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/scratch/mleukam/temp | \
bwa mem -M -t 28 -p GRCm38p6_ref.fa /dev/stdin | \
java -Xmx16G -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=A20_unaligned.bam \
OUTPUT=A20_piped.bam \
R=/group/kline-lab/ref/GRCm38p6_ref.fa \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/scratch/mleukam/temp








