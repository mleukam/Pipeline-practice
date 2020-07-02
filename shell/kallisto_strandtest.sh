#!/bin/bash

#################################
## KALLISTO COUNT READS SCRIPT ##
#################################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0

## Load necessary modules
module load kallisto/0.44.0 

## Assign variables
pathto=/scratch/mleukam/duke_rna
sample=${1}
fq1=${2}
fq2=${3}
fidx=/gpfs/data/kline-lab/ref/transcriptome/Homo_sapiens.GRCh37.67.cdna.all.kidx

## navigate to directory containing reads
cd ${pathto}

## Make output directories
mkdir -p ${sample}/readcounts

## execute quant program
kallisto quant -i ${fidx} -o ${pathto}/${sample}/readcounts/frstrand_test -t 8 --bias -b 100 --seed=10 --fr-stranded ${fq1} ${fq2}
kallisto quant -i ${fidx} -o ${pathto}/${sample}/readcounts/rfstrand_test -t 8 --bias -b 100 --seed=10 --rf-stranded ${fq1} ${fq2}
kallisto quant -i ${fidx} -o ${pathto}/${sample}/readcounts/nostrand_test -t 8 --bias -b 100 --seed=10 ${fq1} ${fq2}

## exit code
exit 0