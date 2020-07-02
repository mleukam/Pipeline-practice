#!/bin/bash

############################################
## SUBMISSION SCRIPT FOR CBS SEGMENTATION ##
############################################

## To be run on a head node
## Will spawn a PBS script for each individual file

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Navigate to folder containing script
cd /home/mleukam/shell/

## assign input file to a variable name
file=/gpfs/data/kline-lab/inputs/snpfiletable.tsv

## loop whose input is awk-parsed sample name from column 2 of tsv assigned to ${file}
for sample in `awk -F"\t" '{print $2}' $file | sort | uniq`; do 
	## exported variables will be passed to PBS script by using qsub -V
	## assign columns of tsv to variable (sample and fastq file names)
	export sample=`awk -F"\t" -v s=$sample '$2==s {print $2}' $file`
	export infile=`awk -F"\t" -v s=$sample '$2==s {print $1}' $file`
	echo "submitting CBS segmentation script for "${sample}
	echo "filename is "${infile}
	## invoke run script
	qsub -S /bin/bash -V -N CBS_${sample} \
	-l walltime=1:00:00 -l nodes=1:ppn=1 -l mem=8gb \
	-j oe -o /home/mleukam/logs/CBS/CBS_${sample}.log \
	run_CBSseq.sh;
done && exit 0

## exit code if error
exit 1
