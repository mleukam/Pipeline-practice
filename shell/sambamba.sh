#!/bin/bash

######################
## SAMBAMBMA SCRIPT ##
######################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## load compiler
module load dmd/2.072.1

## load module
module load sambamba/0.6.5

## Navigate to folder containing sequence data
cd /scratch/mleukam/duke_wes/

## define variables
in=$1
out=${in/.bam/.rmdup.mapq20.bam}
 
## index bams
sambamba index ${in}

echo "index for ${in} complete"
 
## the first part filters out reads of mapping quality < 20
## then pipe the output on data stream as input for the second part
## second part removes duplicates
sambamba markdup -r -t 8 ${in} /dev/stdout | sambamba view -t 8 -F "mapping_quality >= 20" -f bam /dev/stdin -o ${out}

echo "low quality reads removed from ${in}"
echo "duplicates removed from ${in}"

## exit
echo "script complete for ${in}"
exit 0