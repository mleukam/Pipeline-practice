#!/bin/bash 

##################################################
## SUBMISSION SCRIPT FOR EXCAVATOR DATA PREPARE ##
##################################################

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
	-v infile=${line} \
	-N excavator_dataprep_${line} \
	-l walltime=8:00:00 \
	-l nodes=1:ppn=12 \
	-l mem=8gb \
	-j oe \
	-o /home/mleukam/logs/excavator${line}.log \
	run_exp_dataprepare.pbs
	echo "submitting dataprepare script for $line";
done < /gpfs/data/kline-lab/inputs/exp_prep_list.txt

## exit code if error
exit 0
