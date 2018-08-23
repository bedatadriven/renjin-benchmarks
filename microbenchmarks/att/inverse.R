

a <- rnorm(400*400)
dim(a) <- c(400, 400)

run <- function() {
  solve(a)
}
