
# Creation of a 500x500 Toeplitz matrix (loops)

a <- rep(0, 500*500); 
dim(a) <- c(500, 500)

run <- function() {
    
  b <- a
  
  # Rem: there are faster ways to do this
  # but here we want to time loops (220*220 'for' loops)! 
  for (j in 1:500) {
      for (k in 1:500) {
         b[k,j] <- abs(j - k) + 1
      }
  }
  b
}



