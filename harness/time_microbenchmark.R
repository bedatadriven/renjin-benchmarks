
# This is a script that runs in a fresh process
# and executes 1000 iterations of the benchmark, timing all the runs
#
# The following command line arguments are expected:
# [benchmark script] [results file] 
# 
# This script must run in the benchmark's directory

args <- commandArgs(trailingOnly = TRUE)
benchmarkScript <- args[1]
resultsFile <- args[2]
iterations <- 1000

# Source the benchmark script, loading any libraries and running
# any set up code

env <- new.env()
source(file = benchmarkScript, local = env, chdir = TRUE)

# The benchmark script must also define a run() function that we will time
if(is.null(env$run)) {
  stop(sprintf("microbenchmark %s is missing run() function", benchmarkScript))
}

result <- microbenchmark::microbenchmark(env$run(), times = iterations)

writeLines(text = as.character(result$time), con = resultsFile)