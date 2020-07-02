#!/usr/bin/env Rscript

library("GLAD")
library("tidyverse")

DLBCL10116 <- read_csv("/gpfs/data/kline-lab/NIHdata/gladfiles/DLBCL10116.csv") %>% 
	rename(PosBase = Position)

DLBCL10116 <- as.data.frame(DLBCL10116)
DLBCL10116_CGH <- as.profileCGH(DLBCL10116)

DLBCL10116_res <- daglad(DLBCL10116_CGH, mediancenter=FALSE, normalrefcenter=FALSE, 
	genomestep=FALSE, smoothfunc="haarseg", lkern="Exponential", 
    model="Gaussian", qlambda=0.999, bandwidth=1, base=FALSE, round=1.5,
    lambdabreak=8, lambdaclusterGen=40, param=c(d=6), alpha=0.001, msize=2,
    method="centroid", nmin=1, nmax=8, amplicon=1, deletion=-5, deltaN=0.2,
    forceGL=c(-0.3,0.3), nbsigma=3, MinBkpWeight=0.35, CheckBkpPos=TRUE)

saveRDS(DLBCL10116_res, "/gpfs/data/kline-lab/NIHdata/gladfiles/DLBCL10116_res.rds")
