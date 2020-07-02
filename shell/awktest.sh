#!/bin/bash

## assign input file to a variable name
file=/gpfs/data/kline-lab/sample_list.tsv

## loop whose input is awk-parsed sample name from column 1 of tsv
for sample in `awk -F"\t" '{print $1}' $file | sort | uniq`; do 
	## exported variables will be passed to PBS script by using qsub -V
	## assign columns of tsv to variable (sample and fastq file names)
	export sample=`awk -F"\t" -v s=$sample '$1==s {print $1}' $file`
	export fq1=`awk -F"\t" -v s=$sample '$1==s {print $2}' $file`
	export fq2=`awk -F"\t" -v s=$sample '$1==s {print $3}' $file`
	echo "submitting kallisto script for "${sample}
	echo ${fq1}
	echo ${fq2}
	echo ${sample};
done


