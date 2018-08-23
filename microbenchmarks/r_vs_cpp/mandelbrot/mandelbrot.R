# original source: https://www.ibm.com/developerworks/community/blogs/jfp/entry/How_To_Compute_Mandelbrodt_Set_Quickly?lang=en_us
library(Rcpp)
xmin   <- -0.74877
xmax   <- -0.74872
width  <- 1000
ymin   <- 0.06505
ymax   <- 0.06510
height <- 1000

r1 <- seq(xmin, xmax, length.out = width)
r2 <- seq(ymin, ymax, length.out = height)

mandelbrot <- function(z, maxiter) {
  vc <- z
  for (n in seq_len(maxiter)) {
    if(abs(z) > 2) {
      return(n)
    }
    z = z * z + vc
  }
  return(maxiter)
}

mandelbrot_set <- function(r1, r2, maxiter) {
  res <- matrix(0, nrow = length(r1), ncol = length(r2))
  for (i in seq_len(length(r1))) {
    for (j in seq_len(length(r2))) {
      res[i, j] <- mandelbrot(complex(real = r1[i], imaginary= r2[j]), maxiter)
    }
  }
}



mandelbrot_sep <- function(creal, cimag, maxiter) {
  real = creal
  imag = cimag
  for (n in 1:maxiter) {
    real2 <- real * real
    imag2 <- imag * imag
    if ((real2 + imag2) > 4.0) {
      return(n)
    }
    imag <- 2 * real * imag + cimag
    real <- real2 - imag2 + creal
  }
  0
}
mandelbrot_set_sep <- function(r1, r2, maxiter) {
  n3 <- matrix(0, nrow = length(r1), ncol = length(r2))
  for (i in seq_len(length(r1))) {
    for (j in seq_len(length(r2))) {
      n3[i, j] <- mandelbrot_sep(r1[i], r2[j], maxiter)
    }
  }
  n3
}

if(GNUR) {
  library(Renjin)
  
	t1 <- system.time( m1 <- mandelbrot_set(r1, r2, 2048) )
	t2 <- system.time( m2 <- mandelbrot_set_sep(r1, r2, 2048) )
	t3 <- system.time( m3 <- renjin(mandelbrot_set(r1, r2, 2048) ))
	t4 <- system.time( m4 <- renjin(mandelbrot_set_sep(r1, r2, 2048) ))
	t5 <- system.time( m5 <- mandelbrot_set_cpp(r1, r2, 2048) )

	stopifnot(identical(m1, m2))
	stopifnot(identical(m1, m3))
	stopifnot(identical(m1, m4))
	stopifnot(identical(m1, m5))
	
  timings <- rbind(t1, t2, t3, t4, t5)
  print(timings)
} else {
  library(Renjin)
  
	t1 <- system.time( m1 <- mandelbrot_set(r1, r2, 2048) )
	t2 <- system.time( m2 <- mandelbrot_set_sep(r1, r2, 2048) )
	t3 <- system.time( m3 <- mandelbrot_set_cpp(r1, r2, 2048) )
	
	stopifnot(identical(m1, m2))
	stopifnot(identical(m1, m3))
	
  timings <- rbind(t1, t2, t3)
  print(timings)
}


