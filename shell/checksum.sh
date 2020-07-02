#!/bin/bash
#
# One-off script for generating hash list
#
#
# Navigate to file containing bam files that need hashing
#
cd /scratch/mleukam/wes_data
#
# Generate hash list and pipe to text file
md5sum *final.bam > checksum.txt
