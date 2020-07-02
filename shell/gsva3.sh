###############################
## Shell script for gsva    ##
###############################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0

## Load module
module load R/3.4.1

## Navigate to directory containing R script
cd /home/mleukam/R/Rscripts/

## Call R script
Rscript gsva3.r && exit 0

## Exit if error
exit 1
