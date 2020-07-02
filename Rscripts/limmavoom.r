## -----------------------------------------------------

## Limma Voom with Quality Weights on cluster

## Expects cleaned expression matrix of counts filtered for expression
##    Matrix stored in DGElist
##    TMM Normalization factors stored in DGElist
## Expects design matrix specifying batch effect and groups of interest
##    Gene list input must be a list
## Outputs Elist object from voom and MArrayLM object containing results of linear model fits

## M. Leukam
## 7/21/2019

## -----------------------------------------------------

# clear workspace
rm(list = ls())

# load packages
library("tidyverse")
library("limma")

# read in data
dge <- readRDS("/gpfs/data/kline-lab/inputs/dgelist_for_limma.rds")
design <- readRDS("/gpfs/data/kline-lab/inputs/design_for_limma.rds")

# voom transformation is applied to the normalized and filtered DGEList object
v <- voomWithQualityWeights(dge, 
                            design, 
                            plot = FALSE, 
                            normalize = "quantile")

# fit linear model and estimate DGE
fit <- lmFit(v, design)

# write out results
saveRDS <- saveRDS(v, "/gpfs/data/kline-lab/tcga_macs/output/limma_voom_elist_1.rds")
saveRDS <- saveRDS(fit, "/gpfs/data/kline-lab/tcga_macs/output/limma_voom_marraylmfit.rds")