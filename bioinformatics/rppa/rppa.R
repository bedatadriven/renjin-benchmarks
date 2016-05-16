# Copyright (c) 2015 Ieuan Clay
# based on code from https://github.com/biolion/genbench
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
# RPPA classification

### set up session
set.seed(8008)
DEBUGGING <- FALSE
## packages
require(stats)
#require(rjson)
## global vars
rppa <- read.csv("./data_20160126_rppa.csv", header = TRUE, row.names = NULL,
                 sep = ",", blank.lines.skip = TRUE, stringsAsFactors = FALSE)

# drop non-numeric columns
rows <- rppa$TCGA_ID
rppa <- rppa[ , sapply(rppa, is.numeric)]
row.names(rppa) <- rows

# check it
# str(rppa)
dim(rppa) # 3467  131

if (DEBUGGING) cat("Loading Data: complete\n")

do.dist <- function(input_data) {
  ## compute distance matrix
  if (DEBUGGING) cat("Calculating Distance Matrix\n")
  # transpose input data to get distances between samples, not features
  # convert pearson correlation to distance (i.e. bound 0-1, 0 is close)
  dist_mat <- as.dist( (1 - cor(t(input_data), method="pearson")) / 2 )
  if (DEBUGGING) cat("Calculating Distance Matrix: complete\n")
  return(dist_mat)
}

## unsupervised clustering
do.within.ss <- function(d = NULL, clustering) {
  # cut from 'fpc' package function cluster.stats()
  #   a generalisation of the within clusters sum of squares
  #   (k-means objective function), which is obtained if d is a Euclidean
  #   distance matrix. For general distance measures, this is half the sum of the
  #   within cluster squared dissimilarities divided by the cluster size.

  # variables
  cluster.size <- numeric(0)
  dmat <- as.matrix(d)
  within.cluster.ss <- 0
  di <- list()

  # iterate thought distance matrix
  for (i in 1:max(clustering)) {
    cluster.size[i] <- sum(clustering == i)
    di <- as.dist(dmat[clustering == i, clustering == i])
    if (i <= max(clustering)) {
      within.cluster.ss <- within.cluster.ss + sum(di^2) / cluster.size[i]
    }

  }
  return(within.cluster.ss)
}

do.elbow <- function(df) {
  ## return optimal cluster size according to elbow method
  # this is a simple method, not nessecarily recommended,
  # which looks to find a compromise between the minimal number
  # of clusters and within group variance
  # see: http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
  # see: http://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set
  # EXPECTS: names(df) == c("cluster", "within_ss")

  # check ordering of cluster sizes, small to large
  df <- df[order(df$cluster, decreasing = FALSE) , ]

  # calculate delta to next cluster size
  df$delta <- c(df[1:(nrow(df) -1), "within_ss"] - df[2:nrow(df), "within_ss"], 0)

  # fit curve, to remove small local variability and predict k=1 to k=2
  df$k <- df$cluster - 1
  df$smooth <- predict(lm(delta ~ poly((cluster), 4, raw = TRUE), data = df), data.frame(cluster = df$k))

  # return first cluster number (k) where
  # the difference between the change in within_ss for (k) and (k+1)
  # is less than 10% of the starting (i.e. k=1 to k=2)
  return(min(which(df$smooth < df[1, "smooth"] / 10)))

}

do.hc <- function(dist_mat){
  if (DEBUGGING) cat("Calculating Hierarchical clustering\n")

  require(stats)

  # hierarchical clustering, WARD as linkage
  res <- hclust(d = dist_mat, method="ward")

  # determine clustering statistics (within cluster SS), for a range of 'cuts'
  if (DEBUGGING) cat("Calculating Hierarchical clustering: cutting tree\n")
  cuts <- lapply(2:25, FUN = function(idx){ # note: 1 cluster => 'Inf' error
        cat(paste("    >", idx, "clusters\n"))
        return(data.frame(
               cluster = idx,
               within_ss = do.within.ss(dist_mat, cutree(res, k = idx))
               )
        )
    }
  )

  # determine optimal cut and return labels
  # use within cluster SS to be consistent with kmeans
  best_cut <- do.elbow(do.call("rbind", cuts)) # 8 according to paper
  res <- cutree(res, best_cut)
  res <- data.frame(id = attr(dist_mat, which = "Labels"), cluster = res)
  if (DEBUGGING) cat("Calculating Hierarchical clustering: complete\n")
  return(res)
}

# kmeans
do.km <- function(dist_mat){
  if (DEBUGGING) cat("Calculating K-means clustering\n")
  require(stats)

  # kmeans clustering for a range of cluster numbers
  res <- lapply(2:25, FUN = function(i){
      cat(paste("    >", i, "clusters\n"))
      kmeans(dist_mat, algorithm = "Hartigan-Wong", centers = i)
    }
  )

  # determine clustering statistics for a range of 'cuts'
  cuts <- lapply(res, function(x) {
                                    data.frame(cluster = max(x$cluster),
                                               within_ss = sum(x$withinss)
                                              )
                                  } )

  # determine optimal cut and return labels for this
  best_cut <- do.elbow(do.call("rbind", cuts)) # 8 according to paper
  res <- res[[best_cut]]
  res <- data.frame(id = attr(dist_mat, which = "Labels"), cluster = res$cluster)
  if (DEBUGGING) cat("Calculating K-means clustering: complete\n")
  return(res)
}

dist_mat <- do.dist(input_data = rppa)
hc <- do.hc(dist_mat)
km <- do.km(dist_mat = dist_mat)
print(str(dist_mat))
print(str(hc))
print(str(km))
# final clean up
rm(list=ls())
gc()
