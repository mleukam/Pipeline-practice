# #!/bin/bash
# 
# script for BQSR of short mouse (50bp) single end whole genome reads on Balb/c background
# custom project for mouse tumor sequencing project, no generalizable variables in current version
# specifically, it produces a recalibrated table with the BaseRecalibrator tool
# it then outputs a recalibrated BAM or CRAM file.
# input is aligned BAM files with duplicates marked
# outputs are a recalibration table and a recalibrated BAM
#
# BQSR walker recalibrates base quality scores based on alignment
# BQSR masks sites of known variation and degrades quality of mismatched bases outside regions of known variation
# Best to use a broad database of known variants at this step to prevent unintended degredation
# Source of known variants = Sanger mouse genome project merged variant list (all strains)
#
# Input sequence data for this script is output of mark duplicates script: A20_merged.sorted.bam
#
# set script to fail if any command, variable, or output fails
set -euo pipefail
#
# set IFS to split only on newline and tab
IFS=$'\n\t'
#
# load compiler
module load java-jdk/1.8.0_92
#
# load necessary modules
# NB: the path to GATK and Picard is automatically generated as an environment variable
# location of GenomeAnalysis.jar = ${GATK}
# location of Picard.jar = ${PICARD}
module load gatk/4.0.6.0
module load picard/2.8.1
#
# navigate to working directory
cd /scratch/mleukam/mouse
#
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
--known-sites /scratch/mleukam/mouse/mgp.v5.merged.snps_all.dbSNP142.renamed.vcf \
--known-sites /scratch/mleukam/mouse/mgp.v5.merged.indels.dbSNP142.normed.renamed.vcf \
-O A20_bqsr.table

# apply bqsr table
java -Xmx32G -jar ${GATK} ApplyBQSR \
-R /scratch/mleukam/mouse/GRCm38p6_ref.fa \
-I A20_markduplicates.bam \
--bqsr-recal-file A20_bqsr.table \
-O A20_bqsr.bam

# clean up
rm A20_markduplicates.bam

