# Copyright (c) 2015 Ieuan Clay
# based on code from https://github.com/biolion/genbench
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
#

# reproducibility
set.seed(8008)
DEBUGGING <- FALSE

## packages
library(stats)
library(utils)
## general
do.unpack <- function() {
  if (DEBUGGING) cat("> Start: do.unpack()\n")

  # 'family' data

  fam_files <- data.frame(
		member = c(
      "Daniel MacArthur", "Luke Jostins", "Dan Vorhaus", "Caroline Wright",
			"Kate Morley", "Vincent Plagnol", "Jeff Barrett", "Jan Aerts",
      "Joe Pickrell", "Don Conrad", "Carl Anderson", "Ilana Fisher"),
		dataset_id = c(
      "DGM001", "LXJ001", "DBV001", "CFW001", "KIM001", "VXP001",
      "JCB001", "JXA001", "JKP001", "DFC001", "CAA001", "IPF001"),
		target = c(
      "DGM001_genotypes.zip", "LXJ001_genotypes.zip", "DBV001_genotypes.zip",
      "CFW001_genotypes.zip", "KIM001_genotypes.zip", "VXP001_genotypes.zip",
      "JCB001_genotypes.zip", "JXA001_genotypes.zip", "JKP001_genotypes.zip",
      "DFC001_genotypes.zip", "CAA001_genotypes.zip", "IPF001_genotypes.zip"),
		stringsAsFactors = FALSE
		)

  lapply(
    1:nrow(fam_files),
    function(x) {
          unzip(fam_files[x, "target"])
    })

  if (DEBUGGING) cat("> End: do.unpack()\n")
  return(TRUE)
}

## mutation data

do.pop.load <- function() {
  if (DEBUGGING) cat("> Start: do.pop.load()\n")

  readLines("laml.maf", 3) # no header info other than col names

  maf <- read.delim("laml.maf", blank.lines.skip = TRUE, , stringsAsFactors = FALSE)

  ## sample data
  meta <- read.delim("laml.meta.tsv", blank.lines.skip = TRUE, stringsAsFactors = FALSE)

  stopifnot(length(intersect(unique(maf$TCGA_id), unique(meta$bcr_patient_barcode))) == 197)

  if (DEBUGGING) cat("> End: do.pop.load()\n")
  return(list(maf = maf, meta = meta))
}

do.pop.fig1a <- function(pop) {
  if (DEBUGGING) cat("> Start: do.pop.fig1a()\n")
  ## figure 1A: mutations per sample, split by mutation tier and disease status
  # split-apply-combine
  fig_1a <- do.call(
                "rbind",
                lapply(
                    split(pop$maf, f = pop$maf$TCGA_id),
                    function(df) {
                        tmp <- data.frame(table(df$tier))
                        names(tmp) <- c("tier", "mut_count")
                        tmp$TCGA_id <- df[1, "TCGA_id"]
                        return(tmp)})
                )

  # merge metadata
  fig_1a <- merge(x = fig_1a, y = pop$meta,
                  by.x = "TCGA_id", by.y = "bcr_patient_barcode",
                  all.x = TRUE)
  # subset data
  fig_1a <- subset(fig_1a, tier == "tier1")
  # plot
  plot(y = fig_1a$mut_count,
       x = factor(fig_1a$acute_myeloid_leukemia_calgb_cytogenetics_risk_category))

  # return
  if (DEBUGGING) cat("> End: do.pop.fig1a()\n")
  return(fig_1a[ , c("TCGA_id", "tier", "mut_count")])
}

do.pop.fig1b <- function(pop) {
  if (DEBUGGING) cat("> Start: do.pop.fig1b()\n")
  ## figure 1B: samples per mutated gene
  fig_1b <- do.call(
                "rbind",
                lapply(
                  split(pop$maf, f = pop$maf$gene_name),
                  function(df) {
                      tmp <- data.frame(table(df$tier))
                      names(tmp) <- c("tier", "sample_count")
                      tmp$gene_name <- df[1,"gene_name"]
                      return(tmp)})
                )

  # reorder
  fig_1b <- fig_1b[order(fig_1b$sample_count, decreasing = TRUE), ]
  fig_1b <- head(subset(fig_1b, tier == "tier1"), n = 100) # top 100 tier 1 genes by count

  # return
  if (DEBUGGING) cat("> End: do.pop.fig1a()\n")
  return(fig_1b[ , c("gene_name", "tier", "sample_count")])
}

# "window" apply
# modified from: http://www.r-bloggers.com/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
wapply <- function(x, width, by = NULL, FUN = NULL, ...) {
  if (DEBUGGING) cat("> Start: wapply()\n")
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width

  starts <- seq(1, length(x) - width + 1, by = by)
  indices <- lapply(starts, function(x) x:(x + width - 1))

  windows <- base:::simplify2array(
    lapply(indices, function(a) FUN(x[a], ...)),
    higher = TRUE)
  if (DEBUGGING) cat("> End: wapply()\n")
  return(windows)
}


## familial data
do.fam.load <- function(chromosomes = c(10)) {
  if (DEBUGGING) cat("> Start: do.fam.load()\n")

  DATA_DIR <- normalizePath("./")

  ## load data for individuals from genomesunzipped.org
  # chromosomes : limit chromosomes loaded

  # read in each file, and strip down to chromosome 1 (again, could exapnd later if needed)
  #fam <- lapply(lapply(dir(DATA_DIR, full.names = TRUE, pattern = "genotypes.txt"),
  #           function(x) read.delim(x, skip= 14,
  #                                  blank.lines.skip=TRUE, stringsAsFactors=FALSE)
  #           ),
  #           function(x) subset(x, chromosome %in% chromosomes)
  #           )
  fam <- lapply(
            dir(DATA_DIR, full.names = TRUE, pattern = "genotypes.txt"),
            function(x) {
              read.delim(x, skip= 14, blank.lines.skip=TRUE, stringsAsFactors=FALSE)[1:578320,]
             })

  names(fam) <- dir(DATA_DIR, full.names = FALSE, pattern = "genotypes.txt")

  # return
  if (DEBUGGING) cat("> End: do.fam.load()\n")
  return(fam)
}

do.fam.prepare <- function(fam) {
  if (DEBUGGING) cat("> Start: do.fam.prepare()\n")
  ## scores for individual each SNP when compared between individuals
  # get SNP allele frequencies
  genotypes <- do.call("rbind", lapply(fam, FUN = function(x) table(x$genotype)))
  # create empty score matrix
  scores <- matrix( nrow = ncol(genotypes), ncol = ncol(genotypes) ,
                    dimnames = list(colnames(genotypes), colnames(genotypes)))
  # function for calculating number of alleles
  # allowing for ordering and homozygosity
  do.matching.alleles <- function(i, j) {

    ## returns number of matching alleles for two allele pairs
    # expects i and j to be 2 character strings, i.e. i="AA", j="GT"

    # returns sum of non-zero row sums,
    # when each character in i is compared to each character in j,
    # such that:
    # in case of both pairs being heterozygous:
    # - return 2 if both match,
    # - 1 if one allele matches,
    # - 0 otherwise

    #   > sapply(c("A", "B"), function(x) x == c("A", "B"))
    #           A     B
    #   [A]  TRUE FALSE  -> TRUE
    #   [B] FALSE  TRUE  -> TRUE
    #                        sum = 2
    #   > sapply(c("B", "A"), function(x) x == c("A", "B"))
    #         B     A
    #   [A] FALSE  TRUE  -> TRUE
    #   [B]  TRUE FALSE  -> TRUE
    #                       sum = 2
    #
    #   > sapply(c("A", "B"), function(x) x == c("B", "C"))
    #         A     B
    #   [B] FALSE  TRUE  -> TRUE
    #   [C] FALSE FALSE  -> FALSE
    #                       sum = 1
    #
    # in case of both pairs being homozygous:
    # - return 2 if identical,
    # - 0 otherwise

    #   > sapply(c("B", "B"), function(x) x == c("B", "B"))
    #         B    B
    #   [1,] TRUE TRUE  -> TRUE
    #   [2,] TRUE TRUE  -> TRUE
    #                     sum = 2
    #   > sapply(c("B", "B"), function(x) x == c("C", "C"))
    #         B     B
    #   [C] FALSE FALSE  -> FALSE
    #   [C] FALSE FALSE  -> FALSE
    #                     sum = 0
    #
    # in case of one pair being homozygous, the other heterozygous:
    # - return 1 if a heterozygous allele matches the homozygous pair
    # - 0 otherwise

    #   > sapply(c("A", "A"), function(x) x == c("A", "B"))
    #           A     A
    #   [A]  TRUE  TRUE  -> TRUE
    #   [B] FALSE FALSE  -> FALSE
    #                     sum = 1
    #   > sapply(c("B", "B"), function(x) x == c("A", "B"))
    #         B     B
    #   [A] FALSE FALSE  -> FALSE
    #   [B]  TRUE  TRUE  -> TRUE
    #                       sum = 1
    #   > sapply(c("A", "A"), function(x) x == c("C", "B"))
    #         A     A
    #   [C] FALSE FALSE  -> FALSE
    #   [B] FALSE FALSE  -> FALSE
    #                       sum = 0
    #
    #   > sapply(c("A", "B"), function(x) x == c("B", "B"))
    #         A    B
    #   [B] FALSE TRUE     -> TRUE
    #   [B] FALSE TRUE     -> TRUE
    #                         rowsum = 2
    #                       ***SHOULD EQUAL 1!!***
    #                           take colsum if second match is homozygous
    #                           colsum = 1

    if (nchar(i) == 2 && nchar(j) == 2) {

      return(

        sum(
          apply(
            sapply(strsplit(i, '')[[1]], # for each character in i
                   function(x) x == strsplit(j, '')[[1]]), # compare to each character in j
            ( length(unique(strsplit(j, '')[[1]])) == 1 ) +1, #colsum if homozygous
               function(x) sum(x) > 0) # sum of matches > 0 ?
          )

        )
    }
    else { return(0) }

  }

  # fill score matrix with precalculated match scores
  # for all possible allele pair comparisons
  lapply(rownames(scores), function(i) {
    lapply(colnames(scores), function(j) {
      scores[i, j] <<-              # fill score matrix
        do.matching.alleles(i, j) * # number of matching alleles
        (                          # normalised to "liklihood" of seeing this pair of allele pairs
          (sum(genotypes) - sum(genotypes[ , c(i, j)])) / sum(genotypes)
          )
    }
    )
  } )
  round(scores, 3)

  # return
  return(scores)
  if (DEBUGGING) cat("> End: do.fam.prepare()\n")
}

do.fam.check <- function(fam) {
  if (DEBUGGING) cat("> Start: do.fam.check()\n")
  ## check that all SNPs match
  checks <- c()

  # all sets have the same number of SNPs
  if (length(unique(lapply(fam, nrow))) == 1) {
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }

  # all SNP rsids the same
  if (
    all(
      sapply(fam[1],
             # check that the list of SNPs matches the first list of SNPs
             # i.e. intersection length is equal to unintersected
             function(x) {
               length(intersect(fam[[1]]$X..rsid, x$X..rsid)) == length(fam[[1]]$X..rsid)
             }
             )
        )
    ) {
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }
  if (DEBUGGING) cat(">>>> DONE: SNP rsid same?\n")
  # all SNPs in the same order
  if (
    all(
      sapply(fam[1],
             # check that the list of SNPs matches the first list of SNPs
             # i.e. intersection length is equal to unintersected
             function(x) all(order(fam[[1]]$X..rsid) == order(x$X..rsid))
      )
    )
  ) {
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }
  if (DEBUGGING) cat(">>>> DONE: SNP order same?\n")

  # all checks good?
  if (DEBUGGING) cat("> End: do.fam.check()\n")
  return(all(checks))
}

do.ibd.vector <- function(fam, scores) {
  if (DEBUGGING) cat("> Start: do.ibd.vector()\n")
  do.ibd <- function(maf1, maf2, scores) {
    if (DEBUGGING) cat(">>> Start: do.ibd()\n")

    ### IBD
    # http://en.wikipedia.org/wiki/Identity_by_descent
    # (browning and browning, 2007)[http://www.sciencedirect.com/science/article/pii/S0002929707638828]
    # (rapid calculation of IBD matrices)[http://www.biomedcentral.com/content/pdf/1297-9686-33-5-453.pdf]

    # combine to 'vector' of paired alleles
    genotype.vec <- cbind(maf1$genotype, maf2$genotype)
    if (DEBUGGING) cat(">>> DONE: cbind()\n")
    # score each pair
    genotype.score <- apply(genotype.vec, 1, function(x) scores[x[1], x[2]])
    if (DEBUGGING) cat(">>> DONE: apply(genotype,score)\n")
    if (DEBUGGING) cat(">>> End: do.ibd()\n")

    # return as annotated df for further analysis
    return(

      cbind(maf1[ , names(maf1) != "genotype"], data.frame(score = genotype.score))

      )

  }

  fam.scores <- lapply(combn(1:length(fam), 2, simplify = FALSE), # all pairwise combinations
                       # calculate ibd score vector
                       function(x) return(
                         list(df = do.ibd(fam[[x[1]]], fam[[x[2]]], scores = scores),
                              pair = names(fam)[x]
                              )
                         )
                       )
  if (DEBUGGING) cat("> End: do.ibd.vector()\n")
  return(fam.scores)
}

do.ibd.window <- function(fam.scores, window.sizes = seq(5, 50, 5)) {
  if (DEBUGGING) cat("> Start: do.ibd.window()\n")

  ## calculate some sliding windows summing the scores

  # calculate sliding window over each comparison,
  # calculating average score per SNP in that window

  do.call(
    "rbind",
    unlist(recursive = FALSE,
    # for a range of window sizes
      lapply(
        window.sizes,
        function(win) {
          # run sliding window, summing for each and dividing by number of markers
          lapply(fam.scores,
            function(x) {
              data.frame(
                width = win, pair = x$pair,
                q3 = summary(
                      wapply(
                        x$df$score, width = win,
                        by = as.integer(win / 2),
                        FUN = sum
                      ) / win )["3rd Qu."] ) } )
        } )
    )
  )


  if (DEBUGGING) cat("> End: do.ibd.window()\n")
}


### reporting
## collect data
#do.unpack()

## population data
# load data and compute matrix
pop <- do.pop.load()

# patient level summaries
do.pop.fig1a(pop = pop)

# gene level summaries
do.pop.fig1b(pop = pop)

# todo: something more advanced once packages are implemented in renjin
# maybe something from: http://cran.r-project.org/web/views/Genetics.html
# for example something from package "gap", like gcontrol()
# http://www.inside-r.org/packages/cran/gap/docs/gcontrol


# clean up
rm(pop)
gc()

## "family" data - comparing between (non related) individuals
# load data and check SNP lists
fam <- do.fam.load()

do.fam.check(fam = fam)

# prepare pre-calculated scoring matrix
scores <- do.fam.prepare(fam = fam)

# create vectors comparing pairwise all individuals
fam.scores <- do.ibd.vector(fam = fam, scores = scores)

# score on sliding window
do.ibd.window(fam.scores = fam.scores)

# final clean up
rm(list = ls())
gc()
