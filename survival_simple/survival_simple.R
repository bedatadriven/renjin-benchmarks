# Code by Parham Solaimani
# Test-case workflow for TCGA Browser Shiny app which is developed in Mitch Levesque lab by Phil Cheng
# Analysis code provided by Phil Cheng.

##### set up session #####
set.seed(1000)

##### loading packages #####
library(survival)

##### Set global vars #####
files = list.files(path = ".", pattern = "txt$")

##### Blocks for timing #####
do.load <- function(){
  return(readRDS("pat.gene.rda"))
}

do.analyse <- function(DATA){
  # Performs calculations and plotting
  pat.gene <- DATA
  m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
  sdf <- survdiff(m.surv ~ pat.gene$gene2)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

  survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
  half <- summary(survplot)$table[,"median"]
  print(half)

}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

DATA <- do.load()


do.analyse(DATA)
