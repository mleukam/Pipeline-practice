#!/bin/bash
# 
# script to align short (50bp) single end whole genome reads
# custom project for mouse tumor sequencing project, no generalizable variables in current version
#
#### INPUTS ####
# 1. Unaligned BAM file containing read information with illumina adapters marked: A20_markilluminaadapters.bam
# 2. Reference sequence for alingment: GRCm38p6_ref.fa
# 3. BWA index files for reference sequence stored in working directory
#
#### OUTPUTS ####
# 1. A20_merged.bam
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
# navigate to working directory
cd /scratch/mleukam/mouse
#
# load necessary modules
module load bwa/0.7.17
module load picard/2.8.1
module load samtools/1.6.0
# 
# revert BAM file temporarily back to fastq
java -Xmx16G -jar ${PICARD} SamToFastq \
I=A20_markilluminaadapters.bam \
FASTQ=A20_temp.fq \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/scratch/mleukam/temp
#
# align sequences with BWA
# note that t flag in bwa is set to 28 for number of cores in each Gardner node
# index file for bwa downloaded with reference genome from NCBI
# the reference genome used in this case is GRCm38 (mm10) patch 6 (most recent)
bwa aln -t 12 GRCm38p6_ref.fa A20_temp.fq > A20_temp.sai
bwa samse -t 12 GRCm38p6_ref.fa A20_temp.sai A20_temp.fq > A20_temp.sam
#
#### FILTER, SORT AND MERGE ####
#
# filter out unmapped and multimapped reads
# these reads are not useful for variant calling and create downstream errors
samtools view -F 4 -q 1 A20_temp.sam > A20_aligned.filtered.sam
#
# merge in picard sequence dictionary created in setup
# create a new file (unsortedtemp.sam) that has both the dictionary and the aligned reads.
# using this template: cat dictionary.sam > unsorted_file.sam && cat file.sam >> unsorted_file.sam
cat GRCm38p6_ref.dict > A20_dict.sam && cat A20_aligned.filtered.sam >> A20_dict.sam
#
# samtools sort creates downstream errors with picard tools
# only recommend sorting using picard tools
java -Xmx16G -jar ${PICARD} SortSam \
I=A20_dict.sam \
O=A20_aligned.query.bam \
SORT_ORDER=queryname
#
java -Xmx16G -jar ${PICARD} SortSam \
I=A20_unaligned.bam \
O=A20_unaligned.query.bam \
SORT_ORDER=queryname
#
# merge aligned bam with ubam to restore headers, quality and read group information
java -Xmx16G -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=A20_aligned.query.bam \
UNMAPPED_BAM=A20_unaligned.query.bam \
OUTPUT=A20_merged.bam \
R=/scratch/mleukam/mouse/GRCm38p6_ref.fa \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/scratch/mleukam/temp
#
# sort merged BAM by coordinate for later analysis
java -Xmx16G -jar ${PICARD} SortSam \
I=A20_merged.sam \
O=A20_merged.sorted.bam \
SORT_ORDER=coordinate
#
#### CLEAN UP ####
#
# remove temporary files
rm A20_temp.fq
rm A20_temp.sai
rm A20_temp.sam
rm A20_aligned.bam
rm A20_filteredtemp.sam
rm A20_unsortedtemp.sam
rm A20_faligned.filtered.sam
rm A20_unaligned.query.bam
rm A20_aligned.query.bam
rm A20_dict.sam