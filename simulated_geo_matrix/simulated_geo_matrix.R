# this is a plain-R version of vanilla_R_benchmark.R, without extra data management packages
# only works with NGENES and NPATIENTS >= 1000

## set up session
rm(list=ls())
set.seed(8008)

# packages
library(biclust)
library(s4vd)
library(irlba)

# data collection
NGENES <- "500"
NPATIENTS <- "500"
GEO <-      paste('GEO-', NGENES, '-', NPATIENTS, '.rds', sep="")
GO <-       paste('GO-', NGENES, '-', NPATIENTS, '.rds', sep="")
GENES <-    paste('GeneMetaData-', NGENES, '-', NPATIENTS, '.rds', sep="")
PATIENTS <- paste('PatientMetaData-', NGENES, '-', NPATIENTS, '.rds', sep="")

# plain-R q&d replacement for acast(A, list(names(A)[1], names(A)[2]))
df2mxc <- function(df) {
  d1 <- factor(df[,1])
  d2 <- factor(df[,2])
  m <- matrix(data=NA, nrow=length(levels(d1)),
    ncol=length(levels(d2)), dimnames=list(levels(d1), levels(d2)))
  m[cbind(d1, d2)] <- df[,3]
  m
}

# plain-R q&d replacement for sparseMatrix(go[,1], go[,2], x=go[,3])
df2mxs <- function(df) {
  d1 <- df[,1]
  d2 <- df[,2]
  m <- matrix(data=NA, nrow=max(d1),
    ncol=max(d2))
  m[cbind(d1, d2)] <- df[,3]
  m
}


regression <- function() {
  cat("> START: regression()\n")

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)
  patients <- readRDS(PATIENTS)

  # subset
  sub_gmd = genes[genes$func < 250,]
  colnames(sub_gmd)[1] = "geneid"

  response = patients[,"drug.response"]

  # join
  A = merge(geo, sub_gmd)[,c("patientid", "geneid", "expression.value")]

  # matrix cast
  A <- df2mxc(A)

  # run regression
  res <- lm.fit(x=A, y=response)

  # report data.frame of coefficients
  res <- data.frame(coeff=as.character(names(res$coefficients)),
                    p=res$coefficients)
  cat("> END: regression()\n")
  return(res)
}

covariance <- function() {
  cat("> START: covariance()\n")

  ### Data Management ops start ###

  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)
  patients <- readRDS(PATIENTS)

  sub_pmd <- patients[patients$disease==5,]

  # convert to data tables
  colnames(sub_pmd)[1] = "patientid"

  # join
  A <- merge(geo, sub_pmd)[,c("patientid", "geneid", "expression.value")]

  # convert to matrix
  A <- df2mxc(A)

  # calculate covariance
  covar <- stats::cov(A)

  covar <- which(covar>0.75*(max(covar)), arr.ind=T)

  # return 2 column data.frame
  res <- data.frame(covar)
  res[,1] <- as.character(res[,1])
  cat("> END: covariance()\n")
  return(res[complete.cases(res),])

}

biclustering <- function() {
  cat("> START: biclustering()\n")

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  patients <- readRDS(PATIENTS)

  sub_pmd <- patients[patients$gender == 1 & patients$age <= 40, ]
  colnames(sub_pmd)[1] <- "patientid"
  A <- merge(geo, sub_pmd)[,c("patientid", "geneid", "expression.value")]
  A <- df2mxc(A)

  # run biclustering
  res <- biclust(A, method=BCssvd, K=5)
  # report a data.frame
  res <- data.frame(id=as.character(rownames(A)), clust=biclust::writeclust(res))
  cat("> END: biclustering()\n")
  return(res[complete.cases(res),])

}

svd_irlba <- function() {
  cat("> START: svd_irlba()\n")
  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)

  sub_gmd <- genes[genes$func < 250,]

  # convert to data tables
  colnames(sub_gmd)[1] = "geneid"
  # join
  A <- merge(geo, sub_gmd)[,c("patientid", "geneid", "expression.value")]

  # store as matrix
  A <- df2mxc(A)

  # run svd
  res <- irlba(A, nu=50, nv=50) # compute largest singular values

  # return dataframe
  res <- data.frame(id=as.character(1:length(res$d)), sv=res$d)
  cat("> END: svd_irlba()\n")
  return(res[complete.cases(res),])
}

stats <- function(percentage=1) {
  cat("> START: stats()\n")
  # [percentage] = percentage of rows and columns to runs stats on
  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  go       <- readRDS(GO)

  # update code to start all ids at 1
  geo[,1] <- geo[,1]+1
  geo[,2] <- geo[,2]+1

  # select subset of patients, but breaks if we do. why?? too few
  #geo <- geo[geo$patientid < 0.0025 * max(geo$patientid),]
  A <- df2mxc(geo)

  go[,1] <- go[,1] + 1
  go[,2] <- go[,2] + 1
  go <- df2mxs(go)

  # run comparisons
  res <- lapply(1:as.integer((dim(go)[2]/100)*percentage), function(ii){
    lapply(1:as.integer((dim(A)[1]/100)*percentage), function(jj){
      set1 <- A[jj, go[,ii] == 1]
      set2 <- A[jj, go[,ii] == 0]

      data.frame(id=paste(ii, jj), p=wilcox.test(set1, set2, alternative="less")$p.value)
    })
  })

  # combine and return results with low p vals
  res <- do.call("rbind",unlist(res, recursive = F))
  res <- subset(res, p < 1e-3)
  cat("> END: stats()\n")
  return(res[complete.cases(res),])
}

### reporting of timings
# regression
regression()

# svd
svd_irlba()

# covariance
covariance()

# biclust
biclustering()

# stats
stats(percentage=1)

# final clean up
rm(list=ls())
gc()
