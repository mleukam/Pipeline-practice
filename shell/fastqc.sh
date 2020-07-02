#!/bin/bash

###################
## FASTQC SCRIPT ##
###################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load java-jdk/1.8.0_92

## Load necessary modules
module load fastqc/0.11.5

## Navigate to folder containing sequence data
cd /scratch/mleukam/duke_rna/

## Define variables
fq=${1}

## Note: fastqc won't create the output directory; has to be done beforehand
## Creates zip file and html file in the output directory
## fastqc documentation suggests no more than 6 threads on 32bit OS
## need at least 0.25G memory per thread
fastqc -t 6 -o /scratch/mleukam/fastqc ${fq}

## Exit
echo "script complete for ${fq}"
exit 0
