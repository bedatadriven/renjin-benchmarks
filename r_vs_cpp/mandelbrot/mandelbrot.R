# original source: https://www.ibm.com/developerworks/community/blogs/jfp/entry/How_To_Compute_Mandelbrodt_Set_Quickly?lang=en_us
library(Rcpp)

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

mandelbrot_set <- function(xmin, xmax, ymin, ymax, width, height, maxiter) {
  r1 <- seq(xmin, xmax, length.out = width)
  r2 <- seq(ymin, ymax, length.out = height)
  res <- matrix(0, nrow = width, ncol = height)
  for (i in seq_len(width)) {
    for (j in seq_len(height)) {
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
mandelbrot_set_sep <- function(xmin, xmax, ymin, ymax, width, height, maxiter) {
  r1 <- seq(xmin, xmax, length.out = width)
  r2 <- seq(ymin, ymax, length.out = height)
  n3 <- matrix(0, nrow = width, ncol = height)
  for (i in seq_len(width)) {
    for (j in seq_len(height)) {
      n3[i, j] <- mandelbrot_sep(r1[i], r2[j], maxiter)
    }
  }
  n3
}


# TODO: implement C++ version

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
	t1 <- system.time( m1 <- mandelbrot_set(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) )
	t2 <- system.time( m2 <- mandelbrot_set_sep(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) )
	t3 <- system.time( m3 <- renjin(mandelbrot_set(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) ))
	t4 <- system.time( m4 <- renjin(mandelbrot_set_sep(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) ))
	# TODO: time C++ func

	stopifnot(identical(m1, m2))
	stopifnot(identical(m1, m3))
	stopifnot(identical(m1, m4))
	# TODO: check C++ result
	# TODO: print C++ time
  timings <- rbind(t1, t2, t3, t4)
  print(timings)
} else {
  library(Renjin)
  
	t1 <- system.time( m1 <- mandelbrot_set(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) )
	t2 <- system.time( m2 <- mandelbrot_set_sep(-0.74877,-0.74872,0.06505,0.06510,1000,1000,2048) )
	# TODO: time C++ func

	stopifnot(identical(m1, m2))
	# TODO: check C++ result
	# TODO: print C++ time
  timings <- rbind(t1, t2)
  print(timings)
}


