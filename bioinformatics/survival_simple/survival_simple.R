#
# Copyright (c) 2015 Phil Cheng
# Copyright (c) 2016 BeDataDriven B.V.
# License: ...
#

##### set up session #####
set.seed(1000)
DEBUGGING <- FALSE
##### loading packages #####
library(survival)
##### Set global vars #####
files = list.files(path = ".", pattern = "txt$")

##### Blocks for timing #####
do.load <- function(){
  return(readRDS("pat.gene.rda"))
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
  print(half)
  if (DEBUGGING) cat("> END: do.analyse()\n")

}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

DATA <- do.load()


res <- do.analyse(DATA)

if (DEBUGGING) cat("print(str(res)):\n")
if (DEBUGGING) print(str(res))
if (DEBUGGING) cat("-------------end file ---------\n")
