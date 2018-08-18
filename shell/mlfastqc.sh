## Script to run fastq on a group of files with a specific output
## Built specifically for WES data that was transferred 8/2/18
## First line = shebang to specify interpretor (bash)
## Before using, use chmod to make executable

#!/bin/bash

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t'

## Load compiler
module load java-jdk/1.8.0_92

## Load necessary modules
module load fastqc/0.11.5

## Navigate to directory containing fastq files
cd /scratch/mleukam/james_wes/

## loop to run fastq on each fq file in directory
## note: fastqc won't create the output directory; has to be done beforehand
## creates zip file and html file in the output directory
for a in *.fq; 
do 
	fastqc -o /scratch/mleukam/james_wes/fastqc $a;
done
