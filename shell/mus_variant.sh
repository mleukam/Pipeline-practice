#!/bin/bash

##########
# README #
##########

# the purpose of this script is to call somatic variants in a non-standard mouse tumor
# this version of the script is customized for a single sample and is not for general use

# inputs:
# 1. tumor sequences: A20_bqsr.bam  (output from alignment and bqsr scripts)
# 2. "normal pain" sequences:
# 3. reference sequence:
# 4. panel of normal inputs:
# 5.    


#########
# SETUP #
#########

# Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t'

## Load compiler
module load java-jdk/1.8.0_92

## Load necessary modules
module load gatk/4.0.6.0