
# This is a script that runs in a fresh process
# and executes a single iteration of the microbenchmark in order to verify it's output
#
# The following command line arguments are expected:
# [benchmark script] [rds output file]
# 
# This script must run in the benchmark's directory

args <- commandArgs(trailingOnly = TRUE)
benchmarkScript <- args[1]
outputFile <- args[2]

# Source the benchmark script, loading any libraries and running
# any set up code

env <- new.env()
source(file = benchmarkScript, local = env, chdir = TRUE)

# The benchmark script must also define a run() function that we will time
if(is.null(env$run)) {
  stop(sprintf("microbenchmark %s is missing run() function", benchmarkScript))
}

result <- env$run()

saveRDS(result, file = outputFile)