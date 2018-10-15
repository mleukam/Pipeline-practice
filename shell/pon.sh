##!/bin/bash
#
# script for creation of a panel of normals
# 
#### INPUTS ####
# 1. 50 germline CRAM files selected with random generator from 1000 genome project
# 2. Specific reference genome used to generate CRAM files from EBI
# 3. Reference genome used to align tumor samples
#
# set script to fail if any command, variable, or output fails
set -euo pipefail
#
# set IFS to split only on newline and tab
IFS=$'\n\t' 
#
# load compilers
module load java-jdk/1.8.0_92
#
# load modules
module load gatk/4.0.6.0
# 
# navigate to the directory containing panel of normal input files
cd /gpfs/data/kline-lab/ref/1000G_PoN/
#
# gather the desired input files in the directory as an array
BAMLIST=($(ls *.bam))

# pull the sample name from the input file names and make new array
SMLIST=(${BAMLIST[*]%.*})
#
# run Mutect2 in tumor-only mode to generate VCF of germline variants compared to reference
# use same reference that was used for alignment of tumor sequences
for SAMPLE in ${SMLIST[*]};
do
	java -Xmx16G -jar ${GATK} Mutect2 \
	-R /gpfs/data/kline-lab/ref/GRCh38_full_plus_decoy.fa \
	-I ${SAMPLE}.bam \
	-tumor ${SAMPLE%%.*} \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-O ${SAMPLE}.vcf.gz;
done
#
# collect all the output VCFs into a single panel of normals
java -Xmx16G -jar ${GATK} CreateSomaticPanelOfNormals \
-vcfs *.vcf.gz \
-O panelofnormals.vcf.gz

