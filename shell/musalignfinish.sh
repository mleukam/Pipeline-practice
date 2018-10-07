#!/bin/bash
#
# screwed up a path to a jarfile in an earlier fersion of musalign.sh and needed to finish
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
# navigate to directory
cd /scratch/mleukam/mouse
#
#### ALIGN SEQUENCES ####
#
# load necessary modules
module load bwa/0.7.17
module load picard/2.8.1
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
R=/scratch/mleukam/mouse/GRCm38p6_ref.fa \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/scratch/mleukam/temp