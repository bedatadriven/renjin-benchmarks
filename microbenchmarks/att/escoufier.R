
# Calculation of Escoufier's equivalent vectors


# Calculate the trace of a matrix (sum of its diagonal elements)
Trace <- function(y) {
    sum( c(y)[1 + 0:(min(dim(y)) - 1) * (dim(y)[1] + 1)], na.rm=FALSE)
}

numVariables <- 100

x <- abs(rnorm(numVariables*numVariables)); 
dim(x) <- c(numVariables, numVariables)

run <- function() {
  p <- ncol(x)
  vt <- 1:p                                  # Variables to test
  vr <- NULL                                 # Result: ordered variables
  RV <- 1:p                                  # Result: correlations
  vrt <- NULL
  for (j in 1:p) {                           # loop on the variable number
    Rvmax <- 0
    for (k in 1:(p-j+1)) {                   # loop on the variables
      x2 <- cbind(x, x[,vr], x[,vt[k]])
      R <- cor(x2)                           # Correlations table
      Ryy <- R[1:p, 1:p]
      Rxx <- R[(p+1):(p+j), (p+1):(p+j)]
      Rxy <- R[(p+1):(p+j), 1:p]
      Ryx <- t(Rxy)
      rvt <- Trace(Ryx %*% Rxy) / sqrt(Trace(Ryy %*% Ryy) * Trace(Rxx %*% Rxx)) # RV calculation
      if (rvt > Rvmax) {
        Rvmax <- rvt                         # test of RV
        vrt <- vt[k]                         # temporary held variable
      }
    }
    vr[j] <- vrt                             # Result: variable
    RV[j] <- Rvmax                           # Result: correlation
    vt <- vt[vt!=vr[j]]                      # reidentify variables to test
  }
  
  # Return the results to ensure they are not optimized away
  list(vr, RV)
}
  
