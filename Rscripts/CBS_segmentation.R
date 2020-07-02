#!/usr/bin/env Rscript

# load packages
library(DNAcopy)
library(tidyverse)

# assign path variable
pathto <- "/gpfs/data/kline-lab/NIHdata/gladfiles/"

# assign variables passed from function call arguments in command line
infile <- commandArgs(trailingOnly = TRUE)[1]
samplename <- commandArgs(trailingOnly = TRUE)[2]

# read in csv generated from python script extractcols
snps <- read_csv(
	paste0(pathto, infile)
		)

# remove any values with unknown chromosome (NA)
snps <- drop_na(snps, Chromosome)

# convert to CNA object
snps.CNA <- CNA(genomdat = snps$LogRatio,
                      chrom = snps$Chromosome,
                      maploc = snps$Position,
                      data.type = "logratio",
                      sampleid = samplename,
                      presorted = TRUE)

# vignette recommends smoothing before further analysis
smoothed.CNA.snps <- smooth.CNA(snps.CNA)

# run segmentation
segment.smoothed.CNA.snps <- segment(smoothed.CNA.snps, verbose=2)

# save whole segmentation object for future analysis
saveRDS(segment.smoothed.CNA.snps,
         file = paste0(
         	"/gpfs/data/kline-lab/NIHdata/gladfiles/",
         	samplename,
         	".segment.smoothed.CNA.rds")
         )

# pull out output table
output_table <- segment.smoothed.CNA.snps$output

# save csv of output table
write_csv(output_table, paste0(
	"/gpfs/data/kline-lab/NIHdata/gladfiles/",
	samplename,
	".segment.smoothed.CNA.output.csv")
)
