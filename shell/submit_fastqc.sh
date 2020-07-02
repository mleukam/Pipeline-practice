#!/bin/bash

##################################
## SUBMISSION SCRIPT FOR FASTQC ##
##################################

## To be run on a head node
## Will spawn a PBS script for each individual file

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Navigate to folder containing script
cd /home/mleukam/pbs/

## loop to input files and create child processes
while read -r line || [[ -n "$line" ]]; 
do
	qsub -S /bin/bash \
	-v file=${line} \
	-N fastqc_${line} \
	-l walltime=1:00:00 \
	-l nodes=1:ppn=6 \
	-l mem=4gb \
	-j oe \
	-o /home/mleukam/logs/fastqc/fastqc_${line}.log \
	run_fastqc.pbs
	echo "submitting bam2fq script for $line";
done < /gpfs/data/kline-lab/fq_files.txt

## exit code
exit 0