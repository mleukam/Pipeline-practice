#!/bin/bash

#######################################
## Shell script for CBS segmentation ##
#######################################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Load compiler
module load gcc/6.2.0

## Load module
module load R/3.5.0

## Navigate to directory containing R script
cd /home/mleukam/R/Rscripts/

## Call R script
Rscript CBS_segmentation.R ${infile} ${sample} && exit 0

## Exit if error
exit 1
