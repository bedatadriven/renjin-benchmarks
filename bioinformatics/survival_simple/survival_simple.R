#
# Copyright (c) 2015 Phil Cheng
# Copyright (c) 2016 BeDataDriven B.V.
# License: ...
#

##### set up session #####
set.seed(1000)
DEBUGGING <- TRUE
##### loading packages #####
library(survival)
##### Set global vars #####
files = list.files(path = ".", pattern = "txt$")

##### Blocks for timing #####
do.load <- function(){
  DATA <- readRDS("pat.gene.rda")
  DATA
}

do.analyse <- function(DATA){
  if (DEBUGGING) cat("> START: do.analyse()\n")
  # Performs calculations and plotting
  pat.gene <- DATA
  m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
  if (DEBUGGING) cat(">>> DONE: Surv()\n")
  sdf <- survdiff(m.surv ~ pat.gene$gene2)
  if (DEBUGGING) cat(">>> DONE: survdiff()\n")
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
  if (DEBUGGING) cat(">>> DONE: survfit()\n")
  half <- summary(survplot)$table[,"median"]
  if (DEBUGGING) cat("> END: do.analyse()\n")
  res <- list(pat.gene, m.surv, sdf, p.val, survplot, half)
}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

DATA <- do.load()


res <- do.analyse(DATA)

print(res)
