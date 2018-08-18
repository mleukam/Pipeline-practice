## Make BWA index files
## Input is reference genome in FAST-A format
## Output is 5 BWA index files
## Designed for batch submission to Gardner HPC at UChicago
## First line = shebang to specify interpretor (bash)
## Before using, use chmod to make executable

#!/bin/bash

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t'

## Load compilers
module load gcc/6.2.0

## Load necessary modules
module load bwa/0.7.17

## Navigate to directory containing reference genome
cd /group/kline-lab/ref/

## Flag -a bwtsw indicates to optimize for whole genome
bwa index -a bwtsw GRCh38_full_plus_decoy.fa
