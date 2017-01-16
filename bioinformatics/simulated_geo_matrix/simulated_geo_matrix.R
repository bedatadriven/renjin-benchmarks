START_WORKFLOW <- as.numeric(Sys.time())
# Copyright (c) 2015 MIT DB Group
# based on code from https://github.com/mitdbg/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R
# Copyright (c) 2015 Hannes Muehleisen
# based on code from https://github.com/hannesmuehleisen/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R
# Copyright (c) 2015 Ieuan Clay
# based on code from https://github.com/biolion/genbench
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
#
# # this is a plain-R version of vanilla_R_benchmark.R, without extra data management packages
# only works with NGENES and NPATIENTS >= 1000

## set up session
set.seed(8008)
DEBUGGING <- FALSE

# load packages
library(biclust)
library(s4vd)
library(irlba)

# data collection
load.data <- function() {
  GEO <-      readRDS("GEO-500-500.rds")
  GO <-       readRDS("GO-500-500.rds")
  GENES <-    readRDS("GeneMetaData-500-500.rds")
  PATIENTS <- readRDS("PatientMetaData-500-500.rds")
  data <- list(geo = GEO, go = GO, genes = GENES, patients = PATIENTS)
  data
}

# plain-R q&d replacement for acast(A, list(names(A)[1], names(A)[2]))
df2mxc <- function(df) {
  d1 <- factor(df[ , 1])
  d2 <- factor(df[ , 2])
  m <- matrix(data = NA, nrow = length(levels(d1)),
    ncol = length(levels(d2)), dimnames = list(levels(d1), levels(d2)))
  m[cbind(d1, d2)] <- df[ , 3]
  m
}

# plain-R q&d replacement for sparseMatrix(go[,1], go[,2], x=go[,3])
df2mxs <- function(df) {
  d1 <- df[ , 1]
  d2 <- df[ , 2]
  m <- matrix(data = NA, nrow = max(d1),
    ncol = max(d2))
  m[cbind(d1, d2)] <- df[ , 3]
  m
}


regression <- function(data) {
  if (DEBUGGING) cat("> START: regression()\n")

  ### Data Management ops start ###
  geo      <- data$geo
  genes    <- data$genes
  patients <- data$patients

  # subset
  sub_gmd = genes[genes$func < 250, ]
  colnames(sub_gmd)[1] = "geneid"

  response = patients[ , "drug.response"]

  # join
  A = merge(geo, sub_gmd)[ , c("patientid", "geneid", "expression.value")]

  # matrix cast
  A <- df2mxc(A)

  # run regression
  res <- lm.fit(x = A, y = response)

  # report data.frame of coefficients
  res <- data.frame(coeff = as.character(names(res$coefficients)),
                    p = res$coefficients)
  if (DEBUGGING) cat("> END: regression()\n")
  return(res)
}

covariance <- function(data) {
  if (DEBUGGING) cat("> START: covariance()\n")

  ### Data Management ops start ###

  geo      <- data$geo
  genes    <- data$genes
  patients <- data$patients

  sub_pmd <- patients[patients$disease == 5, ]

  # convert to data tables
  colnames(sub_pmd)[1] = "patientid"

  # join
  A <- merge(geo, sub_pmd)[ , c("patientid", "geneid", "expression.value")]

  # convert to matrix
  A <- df2mxc(A)

  # calculate covariance
  covar <- stats::cov(A)

  covar <- which(covar > 0.75 * (max(covar)), arr.ind = T)

  # return 2 column data.frame
  res <- data.frame(covar)
  res[ , 1] <- as.character(res[ , 1])
  if (DEBUGGING) cat("> END: covariance()\n")
  return(res[complete.cases(res), ])

}

biclustering <- function(data) {
  if (DEBUGGING) cat("> START: biclustering()\n")

  ### Data Management ops start ###
  geo      <- data$geo
  patients <- data$patients

  sub_pmd <- patients[patients$gender == 1 & patients$age <= 40, ]
  colnames(sub_pmd)[1] <- "patientid"
  A <- merge(geo, sub_pmd)[ , c("patientid", "geneid", "expression.value")]
  A <- df2mxc(A)

  # run biclustering
  res <- biclust(A, method = BCssvd, K = 5)
  # report a data.frame
  res <- data.frame(id = as.character(rownames(A)), clust = biclust::writeclust(res))
  if (DEBUGGING) cat("> END: biclustering()\n")
  return(res[complete.cases(res), ])

}

svd_irlba <- function(data) {
  if (DEBUGGING) cat("> START: svd_irlba()\n")
  ### Data Management ops start ###
  geo      <- data$geo
  genes    <- data$genes

  sub_gmd <- genes[genes$func < 250, ]

  # convert to data tables
  colnames(sub_gmd)[1] = "geneid"
  # join
  A <- merge(geo, sub_gmd)[ , c("patientid", "geneid", "expression.value")]

  # store as matrix
  A <- df2mxc(A)

  # run svd
  res <- irlba(A, nu = 50, nv = 50) # compute largest singular values

  # return dataframe
  res <- data.frame(id = as.character(1:length(res$d)), sv = res$d)
  if (DEBUGGING) cat("> END: svd_irlba()\n")
  return(res[complete.cases(res), ])
}

stats <- function(data, percentage = 1) {
  if (DEBUGGING) cat("> START: stats()\n")
  # [percentage] = percentage of rows and columns to runs stats on
  ### Data Management ops start ###
  geo      <- data$geo
  go       <- data$go

  # update code to start all ids at 1
  geo[ , 1] <- geo[ , 1] + 1
  geo[ , 2] <- geo[ , 2] + 1

  # select subset of patients, but breaks if we do. why?? too few
  # geo <- geo[geo$patientid < 0.0025 * max(geo$patientid),]
  A <- df2mxc(geo)

  go[ , 1] <- go[ , 1] + 1
  go[ , 2] <- go[ , 2] + 1
  go <- df2mxs(go)

  # run comparisons
  res <- lapply(1:as.integer( (dim(go)[2] / 100) * percentage ), function(i) {
              s1 <- A[go[ , i] == 1, ]
              s2 <- A[go[ , i] == 0, ]
              data.frame(id = i, p = wilcox.test(s1, s2, alternative = "less")$p.value)
            })

  # combine and return results with low p vals
  res <- sapply(res, rbind)
  res <- as.data.frame(t(res))
  print(head(res))
  res <- subset(res, p < 1e-3)
  if (DEBUGGING) cat("> END: stats()\n")
  return(res)
}

### reporting of timings
# load data
data <- load.data()
# regression
regression_res <- regression(data)

# svd
svd_irlba_res <- svd_irlba(data)

# covariance
covariance_res <- covariance(data)

# biclust
biclustering_res <- biclustering(data)

# stats
stats_res <- stats(data, percentage = 1)

print(svd_irlba_res)
print(regression_res)
print(covariance_res)
print(biclustering_res)
print(stats_res)

END_WORKFLOW <- as.numeric(Sys.time())
TOTAL_TIME <- END_WORKFLOW - START_WORKFLOW
print(TOTAL_TIME)
write(TOTAL_TIME,file="TIMINGS",append=TRUE)

# final clean up
rm(list = ls())
gc()
