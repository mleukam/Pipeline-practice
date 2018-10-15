#!/bin/bash

##########
# README #
##########

# align short (50bp) single end whole genome reads from mouse tumor
# mark illumina adapters, filter unmapped reads and multimapped reads, index and sort
# this script is customized for A20 mouse WGS experiment 10/2018 not for general use
# intended to run on gardner HPC with PBS wrapper
# before using, make executable with chmod

#### INPUTS
# 1. Fastq file containing all reads for sample: A20.fq
# 2. Reference sequence for alingment (mm10): GRCm38p6_ref.fa
# 3. Known germline SNPs in all 18 sequenced mouse strains: mgp.v5.merged.snps_all.dbSNP142.vcf
# 4. Known germline indels in all 18 sequenced mouse strains: mgp.v5.merged.indels.dbSNP142.normed.vcf

#### OUTPUTS
# 1. A20_bqsr.bam --> ready for variant calling
# 2. A20_unaligned.bam --> intermediate file

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
# NB: the path to picard.jar is automatically generated as an environment variable when picard tool module is loaded
# location of picard.jar = ${PICARD}
module load fastqc/0.11.5
module load picard/2.8.1
module load bwa/0.7.17
module load samtools/1.6.0

# navigate to working directory
cd /scratch/mleukam/mouse

# the reference genome used in this case is GRCm38 (mm10) patch 6 (most recent)
# generate bwa index files from reference
bwa index -a bwtsw GRCm38p6_ref.fa 

# generate samtools index from reference
samtools faidx GRCm38p6_ref.fa

# generate picard dictionary from reference
java -jar ${PICARD} CreateSequenceDictionary \
REFERENCE=GRCm38p6_ref.fa \
OUTPUT=GRCm38p6_ref.dict

# get fastqc report on raw sequences before proceeding
# note: fastqc won't create the output directory; has to be done beforehand
# creates zip file and html file in the output directory
fastqc -o /scratch/mleukam/mouse/fastqc A20.fq

#############
# MAKE UBAM #
#############

# convert fastq to unaligned bam file
# include necessary read information
java -Xmx32G -jar ${PICARD} FastqToSam \
FASTQ=A20.fq \
O=A20_unaligned.bam \
READ_GROUP_NAME=HW5FWBBXX.6 \
SAMPLE_NAME=A20 \
LIBRARY_NAME=coreWGS \
PLATFORM_UNIT=K00242 \
PLATFORM=illumina \
SEQUENCING_CENTER=UCHICAGOCORE \
RUN_DATE=2018-09-30T00:00:00-0400 \

# mark Illumina adapters
java -Xmx32G -jar ${PICARD} MarkIlluminaAdapters \
I=A20_unaligned.bam \
O=A20_markilluminaadapters.bam \
M=A20_markilluminaadapters_metrics.txt \
TMP_DIR=/scratch/mleukam/temp

###################
# ALIGN SEQUENCES #
###################

# revert BAM file temporarily back to fastq
java -Xmx32G -jar ${PICARD} SamToFastq \
I=A20_markilluminaadapters.bam \
FASTQ=A20_temp.fq \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/scratch/mleukam/temp

# align sequences with BWA aln+samse (best for single end reads <70bp)
# note that t flag in bwa aln is set to 12, max is 28 threads/cores on single node in Gardner
# will only work if PBS script asks for 12 ppn
bwa aln -t 12 GRCm38p6_ref.fa A20_temp.fq > A20_temp.sai
bwa samse GRCm38p6_ref.fa A20_temp.sai A20_temp.fq > A20_temp.sam

# filter out unmapped and multimapped reads
# these reads are not useful for variant calling and create downstream errors
samtools view -F 4 -q 1 A20_temp.sam > A20_aligned.filtered.sam

# merge in picard sequence dictionary created in setup
# create a new file (unsortedtemp.sam) that has both the dictionary and the aligned reads.
# using this template: cat dictionary.sam > unsorted_file.sam && cat file.sam >> unsorted_file.sam
cat GRCm38p6_ref.dict > A20_dict.sam && cat A20_aligned.filtered.sam >> A20_dict.sam

# samtools sort creates downstream errors with picard tools
# only recommend sorting using picard tools
# sort aligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=A20_dict.sam \
O=A20_aligned.query.bam \
SORT_ORDER=queryname

# sort unaligned sample by queryname
java -Xmx32G -jar ${PICARD} SortSam \
I=A20_unaligned.bam \
O=A20_unaligned.query.bam \
SORT_ORDER=queryname

# merge aligned bam with ubam to restore headers, quality and read group information
# merging also removes hardclips from BWA for discordance between best matching kmer and read
java -Xmx32G -jar ${PICARD} MergeBamAlignment \
ALIGNED_BAM=A20_aligned.query.bam \
UNMAPPED_BAM=A20_unaligned.query.bam \
OUTPUT=A20_merged.bam \
R=/scratch/mleukam/mouse/GRCm38p6_ref.fa \
CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/scratch/mleukam/temp

# sort merged BAM by coordinate for later analysis
java -Xmx32G -jar ${PICARD} SortSam \
I=A20_merged.sam \
O=A20_merged.sorted.bam \
SORT_ORDER=coordinate

####################################
# BASE QUALITY SCORE RECALIBRATION #
####################################

# mark duplicates
java -Xmx32G -jar ${PICARD} MarkDuplicates \
INPUT=A20_merged.sorted.bam \
OUTPUT=A20_markduplicates.bam \
METRICS_FILE=A20_markduplicates_metrics.txt \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
CREATE_INDEX=true \
TMP_DIR=/scratch/mleukam/temp

# generate bqsr table
java -Xmx32G -jar ${GATK} BaseRecalibrator \
-R /scratch/mleukam/mouse/GRCm38p6_ref.fa \
-I A20_markduplicates.bam \
--known-sites /scratch/mleukam/mouse/mgp.v5.merged.snps_all.dbSNP142.vcf \
--known-sites /scratch/mleukam/mouse/mgp.v5.merged.indels.dbSNP142.normed.vcf \
-O A20_bqsr.table

# apply bqsr table
java -Xmx32G -jar ${GATK} ApplyBQSR \
-R /scratch/mleukam/mouse/GRCm38p6_ref.fa \
-I A20_markduplicates.bam \
--bqsr-recal-file A20_bqsr.table \
-O A20_bqsr.bam

# clean up
rm A20_markduplicates.bam
rm A20_temp.fq
rm A20_temp.sai
rm A20_temp.sam
rm A20_aligned.bam
rm A20_filteredtemp.sam
rm A20_unsortedtemp.sam
rm A20_aligned.filtered.sam
rm A20_unaligned.query.bam
rm A20_aligned.query.bam
rm A20_dict.sam
rm A20_markilluminaadapters.bam
