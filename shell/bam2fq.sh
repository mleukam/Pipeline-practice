#!/bin/bash

###################
## BAM2FQ SCRIPT ##
###################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## load compiler
module load gcc/6.2.0

## load module
module load samtools/1.6.0

## Navigate to folder containing sequence data
cd /scratch/mleukam/duke_rna/

## define variables
in=$1
prefix=${in%.bam}
out1=${in/%.bam/_1.fq}
out2=${in/%.bam/_2.fq}

echo "modules loaded, starting analysis"
 
## shuffle and group reads together by their names 
## convert to two fastq files with samtools
## discard reads that are not paired
samtools collate -O -@ 8 ${in} ${prefix} | \
samtools fastq -1 ${out1} -2 ${out2} -0 /dev/null -s /dev/null -n -F 0x900 -@ 8 /dev/stdin

echo "${in} converted to ${out1} and ${out2}"

## exit
echo "script complete for ${in}"
exit 0