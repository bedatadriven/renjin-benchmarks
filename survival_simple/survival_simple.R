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
  cat("> START: do.analyse()\n")
  # Performs calculations and plotting
  pat.gene <- DATA
  m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
  cat(">>> DONE: Surv()\n")
  sdf <- survdiff(m.surv ~ pat.gene$gene2)
  cat(">>> DONE: survdiff()\n")
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

  survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
  cat(">>> DONE: survfit()\n")
  half <- summary(survplot)$table[,"median"]
  print(half)
  cat("> END: do.analyse()\n")

}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

DATA <- do.load()


res <- do.analyse(DATA)

cat("print(str(res)):\n")
print(str(res))
cat("-------------end file ---------\n")