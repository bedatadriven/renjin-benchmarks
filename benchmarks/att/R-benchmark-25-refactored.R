# R Benchmark 2.5 (06/2008) [Simon Urbanek]
# version 2.5: scaled to get roughly 1s per test, R 2.7.0 @ 2.6GHz Mac Pro
# R Benchmark 2.4 (06/2008) [Simon Urbanek]
# version 2.4 adapted to more recent Matrix package
# R Benchmark 2.3 (21 April 2004)
# Warning: changes are not carefully checked yet!
# version 2.3 adapted to R 1.9.0
# Many thanks to Douglas Bates (bates@stat.wisc.edu) for improvements!
# version 2.2 adapted to R 1.8.0
# version 2.1 adapted to R 1.7.0
# version 2, scaled to get 1 +/- 0.1 sec with R 1.6.2
# using the standard ATLAS library (Rblas.dll)
# on a Pentium IV 1.6 Ghz with 1 Gb Ram on Win XP pro

# revised and optimized for R v. 1.5.x, 8 June 2002
# Requires additionnal libraries: Matrix, SuppDists
# Author : Philippe Grosjean
# eMail  : phgrosjean@sciviews.org
# Web    : http://www.sciviews.org
# License: GPL 2 or above at your convenience (see: http://www.gnu.org)
#
# Several tests are adapted from the Splus Benchmark Test V. 2
# by Stephan Steinhaus (stst@informatik.uni-frankfurt.de) 
# Reference for Escoufier's equivalents vectors (test III.5):
# Escoufier Y., 1970. Echantillonnage dans une population de variables
# aleatoires r�elles. Publ. Inst. Statis. Univ. Paris 19 Fasc 4, 1-47.
#
# type source("c:/<dir>/R2.R") to start the test

runs <- 3			# Number of times the tests are executed
times <- rep(0, 15); dim(times) <- c(5,3)
require(Matrix)		# Optimized matrix operations
require(SuppDists)	# Optimized random number generators
#Runif <- rMWC1019	# The fast uniform number generator
Runif <- runif
# If you don't have SuppDists, you can use: Runif <- runif
#a <- rMWC1019(10, new.start=TRUE, seed=492166)	# Init. the generator
#Rnorm <- rziggurat	# The fast normal number generator
# If you don't have SuppDists, you can use: Rnorm <- rnorm
#b <- rziggurat(10, new.start=TRUE)	# Init. the generator
Rnorm <- rnorm
options(object.size=100000000)


benchmark <- function(id, block) {
  cumulate <- 0;

  expr <- substitute(block)
  
  
  # Only the 'timed' block within the benchmark is measured
  timed <- function(timingBlock) {
      timing <- system.time({
        cat('evaluating timing block\n')
        timingBlock 
      })[3]
    cumulate <<- cumulate + timing
  }

  # Get the expr as unevaluated so we can repeat
  env <- new.env()
  env$timed <- timed
  
  replicate(runs, {   
    cat('committe\n')
    completed <- tryCatch({
      eval(expr, envir =  env)
      TRUE
    }, error = function(cond) {
      cat(sprintf("%s: ERROR: %s\n", id, cond))
      FALSE
    }, warning = function(cond) {
      TRUE
    })
    if(!completed) {
      cat("Error thrown during test execution\n")
      return
    }
  })

  timing <- cumulate/runs

  if(completed) {
    cat(sprintf("%s,%f\n", id, timing))
  }
}

# (1) Creation, transp., deformation of a 2500x2500 matrix
benchmark('matrix-1', {
  timed({
    a <- matrix(Rnorm(2500*2500)/10, ncol=2500, nrow=2500);
    b <- t(a)
    dim(b) <- c(1250, 5000);
    a <- t(b)
  })
})

# (2) 2400x2400 normal distributed random matrix ^1000
benchmark('matrix-2',  {
  a <- abs(matrix(Rnorm(2500*2500)/2, ncol=2500, nrow=2500));
  timed({ 
    b <- a^1000 
  })
})

# (3) Sorting of 7,000,000 random values
benchmark('matrix-3',  {
  a <- Rnorm(7000000)
  timed({
    b <- sort(a, method="quick")	# Sort is modified in v. 1.5.x
    # And there is now a quick method that better competes with other packages!!!
  })
})

# (4) 2800x2800 cross-product matrix (b = a' * a)
benchmark('matrix-4',  {
  a <- Rnorm(2800*2800); dim(a) <- c(2800, 2800)
  timed({
    b <- crossprod(a)		# equivalent to: b <- t(a) %*% a
  })
})

# (5) Linear regr. over a 3000x3000 matrix (c = a \\ b')
benchmark('matrix-5',  {
  a <- new("dgeMatrix", x = Rnorm(2000*2000), Dim = as.integer(c(2000,2000)))
  b <- as.double(1:2000)
  timed({
    c <- solve(crossprod(a), crossprod(a,b))
  })
  
  # This is the old method
  #a <- Rnorm(600*600); dim(a) <- c(600,600)
  #b <- 1:600
  #invisible(gc())
  #timed({
  #  qra <- qr(a, tol = 1e-7);
  #  c <- qr.coef(qra, b)
  #  #Rem: a little faster than c <- lsfit(a, b, inter=F)$coefficients
  #})[3]
  #cumulate <- cumulate + timing
})

# (1) FFT over 2,400,000 random values
benchmark('matrix-6',  {
  a <- Rnorm(2400000)
  timed({
    b <- fft(a)
  })
})

benchmark('matrix-calc-1',  {
  a <- array(Rnorm(600*600), dim = c(600, 600))
  # Only needed if using eigen.Matrix(): Matrix.class(a)
  timed({
  	b <- eigen(a, symmetric=FALSE, only.values=TRUE)$Value
  	# Rem: on my machine, it is faster than:
  #	 b <- La.eigen(a, symmetric=F, only.values=T, method="dsyevr")$Value
  #	 b <- La.eigen(a, symmetric=F, only.values=T, method="dsyev")$Value
  #  b <- eigen.Matrix(a, vectors = F)$Value
  })
})

# (3) Determinant of a 2500x2500 random matrix 
benchmark('matrix-calc-2',  {
  a <- Rnorm(2500*2500); dim(a) <- c(2500, 2500)
  #Matrix.class(a)
  timed({
    #b <- determinant(a, logarithm=F)
    # Rem: the following is slower on my computer!
    # b <- det.default(a)
    b <- det(a)
  })
})

# (4) Cholesky decomposition of a 3000x3000 matrix
benchmark('matrix-calc-3',  {
  a <- crossprod(new("dgeMatrix", x = Rnorm(3000*3000),
                       Dim = as.integer(c(3000, 3000))))
  #a <- Rnorm(900*900); dim(a) <- c(900, 900)
  #a <- crossprod(a, a)
  timed({
    b <- chol(a)
  })
})

# (5) Inverse of a 1600x1600 random matrix
benchmark('matrix-calc-5',  {
  a <- new("dgeMatrix", x = Rnorm(1600*1600), Dim = as.integer(c(1600, 1600)))
  #a <- Rnorm(400*400); dim(a) <- c(400, 400)
  timed({
  #  b <- qr.solve(a)
    # Rem: a little faster than
    b <- solve(a)
  })
})

# III. Programmation
# (1) 3,500,000 Fibonacci numbers calculation (vector calc)

benchmark('progr-1',  {
  phi <- 1.6180339887498949
  a <- floor(Runif(3500000)*1000)
  timed({
    b <- (phi^a - (-phi)^(-a))/sqrt(5)
  })
})


# (2) Creation of a 3000x3000 Hilbert matrix (matrix calc)
benchmark('progr-2',  {
  timed({
    b <- rep(1:a, a); dim(b) <- c(a, a);
    b <- 1 / (t(b) + 0:(a-1))
    # Rem: this is twice as fast as the following code proposed by R programmers
    # a <- 1:a; b <- 1 / outer(a - 1, a, "+")
  })
})

# (3) Grand common divisors of 400,000 pairs (recursion)
gcd2 <- function(x, y) {
  if (sum(y > 1.0E-4) == 0) 
    x 
  else {
    y[y == 0] <- x[y == 0]
    Recall(y, x %% y)
  }
}

benchmark('progr-3',  {
  a <- ceiling(Runif(400000)*1000)
  b <- ceiling(Runif(400000)*1000)
  timed({	  
    c <- gcd2(a, b)                            # gcd2 is a recursive function
  })
})


# (4) Creation of a 500x500 Toeplitz matrix (loops)
benchmark('progr-4',  {
  b <- rep(0, 500*500)
  dim(b) <- c(500, 500)
  timed({
  	# Rem: there are faster ways to do this
  	# but here we want to time loops (220*220 'for' loops)! 
    for (j in 1:500) {
      for (k in 1:500) {
        b[k,j] <- abs(j - k) + 1
      }
    }
  })
})

# (5) Escoufier's method on a 45x45 matrix (mixed)
cumulate <- 0; p <- 0; vt <- 0; vr <- 0; vrt <- 0; rvt <- 0; RV <- 0; j <- 0; k <- 0;
x2 <- 0; R <- 0; Rxx <- 0; Ryy <- 0; Rxy <- 0; Ryx <- 0; Rvmax <- 0
# Calculate the trace of a matrix (sum of its diagonal elements)
Trace <- function(y) {
  sum(c(y)[1 + 0:(min(dim(y)) - 1) * (dim(y)[1] + 1)], na.rm=FALSE)

}
benchmark('p',  {
  x <- abs(Rnorm(45*45)); dim(x) <- c(45, 45)
  timed({
    # Calculation of Escoufier's equivalent vectors
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
  })
})