#!/bin/bash

############################################
## SUBMISSION SCRIPT FOR SAMBAMBMA SCRIPT ##
############################################

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
	-N sambamba_${line} \
	-l walltime=3:00:00 \
	-l nodes=1:ppn=1 \
	-l mem=4gb \
	-j oe \
	-o /home/mleukam/logs/sambamba/sambamba_${line}.log \
	run_sambamba.pbs
	echo "submitting sambamba script for $line";
done < /gpfs/data/kline-lab/wes_files.txt

## exit code
exit 0