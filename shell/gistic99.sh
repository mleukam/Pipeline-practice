#!/bin/bash

########################
## GISTIC2.0 ANALYSIS ##
########################
## for copy number variation
## requires output of segmentation algorithm such as CBS or GLAD

## SETUP -------------------------------------------------------

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0

## load necessary modules
module load gistic/2.0.23

## Navigate to working directory
cd /gpfs/data/kline-lab/NIHdata/gladfiles/

## VARIABLES -----------------------------------------------------

## output folder created beforehand: gladfiles/gistic_output
hg19_ref=/apps/software/gcc-6.2.0/gistic/2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat
seg=/gpfs/data/kline-lab/NIHdata/gladfiles/seg.tsv
markers=/gpfs/data/kline-lab/NIHdata/gladfiles/markerfile.tsv
gisticout=/gpfs/data/kline-lab/NIHdata/gladfiles/gistic_output

## RUN PROGRAM ---------------------------------------------------

# run gistic
gistic2 -b ${gisticout} \
-seg ${seg} \
-refgene ${hg19_ref} \
-mk ${markers} \
-genegistic 1 \
-ta 0.1 \
-td 0.1 \
-js 4 \
-qvt 0.25 \
-rx 1 \
-conf 0.99 \
-broad 1 \
-brlen 0.98 \
-armpeel 1 \
-savegene 1 \
-fname nih_gistic_results_99 && exit 0

## if error
exit 1
