
# Benchmarking harness for R benchmarks written as simple scripts
# Rscript SOURCE.R warmupCount runCount

args <- commandArgs(trailingOnly = TRUE)
sourceFile <- args[1]
warmupRuns <- as.integer(args[2])
runs <- as.integer(args[3])

cat(sprintf("Running %s with %d warm up runs and %d timing runs\n", sourceFile, warmupRuns, runs))

benchmark <- parse(file = sourceFile)

# Run a few times to get things warmed up
for(i in 1:warmupRuns) {
  cat(sprintf("Warmup run #%d\n", i))
  eval(benchmark, new.env())
}

# Now run the counts for real
# The ceremonhy part here borrowed from benchmark-25 - not
# sure if really needed

cumulate <- 0
for (i in 1:runs) {
  invisible(gc())
  timing <- system.time({
    eval(benchmark, new.env())
  })
  cumulate <- cumulate + timing['elapsed']
}
timing <- cumulate/runs

# Write the result to standard out
cat(sprintf("RESULT = %f\n", timing))

