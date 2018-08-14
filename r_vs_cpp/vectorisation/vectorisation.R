# original source: http://adv-r.had.co.nz/Rcpp.html
library(Rcpp)


GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  vacc1a <- function(age, female, ily) {
    p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
    p <- p * if (female) 1.25 else 0.75
    p <- max(0, p)
    p <- min(1, p)
    p
  }
  
  vacc1 <- function(age, female, ily) {
    n <- length(age)
    out <- numeric(n)
    for (i in seq_len(n)) {
      out[i] <- vacc1a(age[i], female[i], ily[i])
    }
    out
  }
  
  vacc2 <- function(age, female, ily) {
    p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
    p <- p * ifelse(female, 1.25, 0.75)
    p <- pmax(0, p)
    p <- pmin(1, p)
    p
  }
  
  vacc3 <- cppFunction("
                       double vacc3a(double age, bool female, bool ily)
                       {
                       double p = 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily;
                       p = p * (female ? 1.25 : 0.75);
                       p = std::max(p, 0.0);
                       p = std::min(p, 1.0);
                       return p;
                       }
                       
                       NumericVector vacc3(NumericVector age, LogicalVector female, 
                       LogicalVector ily) 
                       {
                       int n = age.size();
                       NumericVector out(n);
                       
                       for(int i = 0; i < n; ++i) {
                       out[i] = vacc3a(age[i], female[i], ily[i]);
                       }
                       
                       return out;
                       }
                       ")

n <- 1000000
age <- rnorm(n, mean = 50, sd = 10)
female <- sample(c(T, F), n, rep = TRUE)
ily <- sample(c(T, F), n, prob = c(0.8, 0.2), rep = TRUE)

t1  = system.time(x1 <- vacc1(age, female, ily = ily))
t1r = system.time(x1r <- renjin(vacc1(age, female, ily = ily)))
t2  = system.time(x2 <- vacc2(age, female, ily = ily))
t2r = system.time(x2r <- renjin(vacc2(age, female, ily = ily)))
t3  = system.time(x3 <- vacc3(age, female, ily = ily))


stopifnot(
  all.equal(x1, x2),
  all.equal(x2, x3),
  all.equal(x1, x1r),
  all.equal(x2, x2r)
)

timings <- rbind(t1, t1r, t2, t2r, t3)

print(timings)

} else {
  
  vacc1a <- function(age, female, ily) {
    p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
    p <- p * if (female) 1.25 else 0.75
    p <- max(0, p)
    p <- min(1, p)
    p
  }
  
  vacc1 <- function(age, female, ily) {
    n <- length(age)
    out <- numeric(n)
    for (i in seq_len(n)) {
      out[i] <- vacc1a(age[i], female[i], ily[i])
    }
    out
  }
  
  vacc2 <- function(age, female, ily) {
    p <- 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily
    p <- p * ifelse(female, 1.25, 0.75)
    p <- pmax(0, p)
    p <- pmin(1, p)
    p
  }
  
  vacc3 <- cppFunction("
                       double vacc3a(double age, bool female, bool ily)
                       {
                       double p = 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily;
                       p = p * (female ? 1.25 : 0.75);
                       p = std::max(p, 0.0);
                       p = std::min(p, 1.0);
                       return p;
                       }
                       
                       NumericVector vacc3(NumericVector age, LogicalVector female, 
                       LogicalVector ily) 
                       {
                       int n = age.size();
                       NumericVector out(n);
                       
                       for(int i = 0; i < n; ++i) {
                       out[i] = vacc3a(age[i], female[i], ily[i]);
                       }
                       
                       return out;
                       }
                       ")

  n <- 1000000
  age <- rnorm(n, mean = 50, sd = 10)
  female <- sample(c(TRUE, FALSE), n, rep = TRUE)
  ily <- sample(c(TRUE, FALSE), n, prob = c(0.8, 0.2), rep = TRUE)
  
  t1 = system.time(x1 <- vacc1(age, female, ily = ily))
  t2 = system.time(x2 <- vacc2(age, female, ily = ily))
  t3 = system.time(x3 <- vacc3(age, female, ily = ily))


  stopifnot(
    all.equal(x1, x2),
    all.equal(x2, x3)
  )

  timings <- rbind(t1, t2, t3)
  
  print(timings)
}

