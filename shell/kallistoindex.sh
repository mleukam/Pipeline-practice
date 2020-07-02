#!/bin/bash

####################
## KALLISTO INDEX ##
####################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0

## Load necessary modules
module load kallisto/0.44.0 

## Navigate to folder containing reference sequence
cd /gpfs/data/kline-lab/ref/Gencode/

## Make index
kallisto index -i gencode.v29.transcripts.kidx gencode.v29.transcripts.fa