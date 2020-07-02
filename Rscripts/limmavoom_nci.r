## -----------------------------------------------------

## Limma Voom with Quality Weights on cluster

## Expects cleaned expression matrix of counts filtered for expression
##    Matrix stored in DGElist
##    TMM Normalization factors stored in DGElist
## Expects design matrix specifying batch effect and groups of interest
##    Gene list input must be a list
## Outputs Elist object from voom and MArrayLM object containing results of linear model fits

## M. Leukam
## 6/8/2020

## -----------------------------------------------------

# clear workspace
rm(list = ls())

# load packages
library("tidyverse")
library("limma")

# read in data
dge <- readRDS("/gpfs/data/kline-lab/inputs/v5_nci_dgelist_for_limma.rds")
design <- readRDS("/gpfs/data/kline-lab/inputs/v5_nci_design_for_limma.rds")

# voom transformation is applied to the normalized and filtered DGEList object
v <- voomWithQualityWeights(dge, 
                            design, 
                            plot = FALSE, 
                            normalize = "quantile")

# fit linear model and estimate DGE
fit <- lmFit(v, design)

# write out results
saveRDS <- saveRDS(v, "/gpfs/data/kline-lab/inputs/nci_limma_voom_elist_1.rds")
saveRDS <- saveRDS(fit, "/gpfs/data/kline-lab/inputs/nci_limma_voom_marraylmfit.rds")
