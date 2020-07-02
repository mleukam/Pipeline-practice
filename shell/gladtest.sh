#!/bin/bash

###############################
## Shell script for tximport ##
###############################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0
module load gsl/2.3

## Load module
module load R/3.5.0

## Navigate to directory containing R script
cd /home/mleukam/R/Rscripts/

## Call R script
Rscript gladtest.R && exit 0

## Exit if error
exit 1
