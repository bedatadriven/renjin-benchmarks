
# original source:http://dirk.eddelbuettel.com/blog/2010/09/07/


# Xian's code, using <- for assignments and passing x down
f <- function(n, x=1) { for (i in 1:n) { x=1/(1+x) }; x }
g <- function(n, x=1) { for (i in 1:n) { x=(1/(1+x)) }; x }
h <- function(n, x=1) { for (i in 1:n) { x=(1+x)^(-1) }; x }
j <- function(n, x=1) { for (i in 1:n) { x={1/{1+x}} }; x }
k <- function(n, x=1) { for (i in 1:n) { x=1/{1+x} }; x }

# According to Dirk, h() is the slowest in R,
# so let's use that as our comparison

# Loading required package: methods
# test replications elapsed relative
# 6 l(N, 1)           10   0.122    1.000
# 5 k(N, 1)           10   9.880   80.984
# 1 f(N, 1)           10   9.978   81.787
# 4 j(N, 1)           10  11.293   92.566
# 2 g(N, 1)           10  12.027   98.582
# 3 h(N, 1)           10  15.372  126.000

run <- function(n = 1e6) {
  h(n)
}



