#!/bin/bash

##########
# README #
##########

# script to set up and format all the necessary files for alignment and variant calling for mouse whole genome tumor sequencing project
# each of these subscripts should only need to be run one time
# custom built for A20 project, not generalizable in current format

#########
# SETUP #
#########

# set script to fail if any command, variable, or output fails
set -euo pipefail

# set IFS to split only on newline and tab
IFS=$'\n\t' 

# load compilers
module load gcc/6.2.0
 
# load modules
# NB: the path to picard.jar is automatically generated as an environment variable when picard tool module is loaded
# location of picard.jar = ${PICARD}
module load bcftools/1.6.0     

# navigate to working directory
cd /scratch/mleukam/mouse

################################################
# RENAME CONTIGS IN GERMLINE VARIANT REFERENCE #
################################################

# VCF files from Sanger Mouse Genome Project contain known SNPs from 18 strains of inbred lab mice including balb/c
# will use all known SNPs for masking in base quality score recalibration (all 18 strains)
# reference genome (GRCm38/mm10) uses NCBI names for chromosomes
# VCF from sanger uses a bare number for chromosomes
# need to rename VCF contigs to match reference sequence
# replace_chr.file is text file containing space-separated columns with each line following format "old_name new_name/n"
# NCBI/refseq chromosome contig names taken from website: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26
# 1 NC_000067.6
# 2 NC_000068.7
# 3 NC_000069.6
# 4 NC_000070.6
# 5 NC_000071.6
# 6 NC_000072.6
# 7 NC_000073.6
# 8 NC_000074.6
# 9 NC_000075.6
# 10 NC_000076.6
# 11 NC_000077.6
# 12 NC_000078.6
# 13 NC_000079.6
# 14 NC_000080.6
# 15 NC_000081.6
# 16 NC_000082.6
# 17 NC_000083.6
# 18 NC_000084.6
# 19 NC_000085.6
# X NC_000086.7
# Y NC_000087.7
# MT NC_005089.1
bcftools annotate --rename-chrs replace_chr.file mgp.v5.merged.snps_all.dbSNP142.vcf > mgp.v5.merged.snps_all.dbSNP142.renamed.vcf
bcftools annotate --rename-chrs replace_chr.file mgp.v5.merged.indels.dbSNP142.normed.vcf > mgp.v5.merged.indels.dbSNP142.normed.renamed.vcf

########################
# GENERATE INDEX FILES #
########################

# the reference genome used in this case is GRCm38 (mm10) patch 6 (most recent)
# generate bwa index files from reference
bwa index -a bwtsw GRCm38p6_ref.fa 

# generate samtools index from reference
samtools faidx GRCm38p6_ref.fa

# generate picard dictionary from reference
java -jar ${PICARD} CreateSequenceDictionary \
REFERENCE=GRCm38p6_ref.fa \
OUTPUT=GRCm38p6_ref.dict