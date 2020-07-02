#!/bin/bash

##############################################
## SUBMISSION SCRIPT FOR GETTING READCOUNTS ##
##############################################

## To be run on a head node
## Will spawn a PBS script for each individual file

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Navigate to folder containing script
cd /home/mleukam/pbs/

## assign input file to a variable name
file=/gpfs/data/kline-lab/sample_list.tsv

## loop whose input is awk-parsed sample name from column 1 of tsv assigned to ${file}
for sample in `awk -F"\t" '{print $1}' $file | sort | uniq`; do 
	## exported variables will be passed to PBS script by using qsub -V
	## assign columns of tsv to variable (sample and fastq file names)
	export sam=`awk -F"\t" -v s=$sample '$1==s {print $1}' $file`
	export fq1=`awk -F"\t" -v s=$sample '$1==s {print $2}' $file`
	export fq2=`awk -F"\t" -v s=$sample '$1==s {print $3}' $file`
	echo "submitting kallisto script for "${sam}
	echo "file 1 is "${fq1}
	echo "file 2 is "${fq2}
	## invoke run script
	qsub -S /bin/bash -V -N kallisto_${sam} \
	-l walltime=6:00:00 -l nodes=1:ppn=8 -l mem=8gb \
	-j oe -o /home/mleukam/logs/kallisto/kallisto_${sam}.log \
	run_kallisto.pbs;
done && exit 0

## exit code
exit 0