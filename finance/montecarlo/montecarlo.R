
# Copyright (c) 2016 BeDataDriven B.V.
# License: http://www.apache.org/licenses/LICENSE-2.0 Apache License version 2.0

cat("Creating dummy data... ")
# Create the set of 124 factors used by Moody's KMV GCORR model:
kmv.factors <- list(
  # country factors:
  country = list(letter = "C", count = 49),
  # industry factors:
  industry = list(letter = "N", count = 61),
  # base factors:
  base = list(letter = "S", count = 14)
  )

factors <- do.call(c, lapply(kmv.factors, function(type) {
  sprintf("%s%02d", type$letter, seq(type$count))
  }))

nfactors <- length(factors)

# Create a large set of dummy obligors:
nobligors <- 10000

nultimates <- sample(100, nobligors, replace = TRUE)
obligors <- data.frame(
  obligor.id = seq(nobligors),
  pod = runif(nobligors),
  number.obligors = nultimates,
  default.distribution = ifelse(nultimates == 1L, "Bernoulli", "Poisson"),
  stringsAsFactors = FALSE
  )

# Create the factor loadings per obligor:
sampleFactors <- function() {
  # Sample 1 country factor, 1 or 2 industry factors and add 14 base factors:
  with(kmv.factors,
       c(sprintf("%s%02d", country$letter, sample(country$count, 1)),
         sprintf("%s%02d", industry$letter, sample(industry$count, ceiling(runif(1, max = 2)))),
         sprintf("%s%02d", base$letter, seq(base$count))
       ))
}

factor.loadings <- do.call(rbind, lapply(seq(nobligors), function(obligor.id) {
  factors <- sampleFactors()
  data.frame(
    obligor.id = rep(obligor.id, length(factors)),
    factor = factors,
    factor.val = rnorm(length(factors)),
    stringsAsFactors = FALSE)
}))
cat("done.\n")

# Create the factor loading matrix:
cat("Creating factor loading matrix... ")
t <- system.time({
  A <- do.call(rbind,
               lapply(split(factor.loadings, factor.loadings$obligor.id), function(obligor) {
                 v <- structure(rep_len(0, length(factors)), .Names = factors)
                 v[obligor$factor] <- obligor$factor.val
                 v
               }))
})
cat("done.\n")
print(t)

nSimulations <- 1e6
nFactorDraws <- 10

cat("Performing Monte Carlo simulations... ")
t <- system.time({
  # Sample the state of the economy (dim(z) = nfactors*nFactorDraws):
  z <- matrix(rnorm(nfactors*nFactorDraws), nrow = nfactors)
  # Calculate the conditional probability (note that the result of qnorm() gets 
  # recycled to match the dimensions of A%*%z; dim(cp) = nobligors*nFactorDraws):
  cp <- pnorm((qnorm(obligors$pod) - A %*% z)/(1 - rowSums(A^2)))
  
  losses <- do.call(rbind, lapply(seq(nSimulations), function(i) {
    # Sample from the default probabilities
    is.bern <- matrix(ifelse(obligors$default.distribution == "Bernouilli", 1, 0), nrow = 1)
    is.pois <- matrix(ifelse(obligors$default.distribution == "Poisson", 1, 0), nrow = 1)
    y.vec <- is.bern %*% matrix(rbinom(length(cp), 1, cp), nrow = nobligors) +
      is.pois %*% matrix(rpois(length(cp), cp*obligors$number.obligors), nrow = nobligors)
    
    y.vec
  }))
})

cat("done.\n")
print(t)

cat("mean losses:\n")
print(rowMeans(losses))


